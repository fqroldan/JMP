function extend_state_space!(h::Hank, qʰ_mat, qᵍ_mat, T_mat, R_mat)

	Npn = length(h.pngrid)

	ϕa_ext = SharedArray{Float64}(h.Nω, h.Nϵ, h.Nb, h.Nμ, h.Nσ, h.Nw, h.Nζ, h.Nz, Npn)
	ϕb_ext = SharedArray{Float64}(h.Nω, h.Nϵ, h.Nb, h.Nμ, h.Nσ, h.Nw, h.Nζ, h.Nz, Npn)
	ϕc_ext = SharedArray{Float64}(h.Nω, h.Nϵ, h.Nb, h.Nμ, h.Nσ, h.Nw, h.Nζ, h.Nz, Npn)

	all_knots = (h.ωgrid, h.ϵgrid, h.bgrid, h.μgrid, h.σgrid, h.wgrid, h.ζgrid, h.zgrid)
	agg_knots = (h.bgrid, h.μgrid, h.σgrid, h.wgrid, h.ζgrid, h.zgrid)
	itp_vf = interpolate(all_knots, h.vf, Gridded(Linear()))
	itp_R  = interpolate(agg_knots, R_mat, Gridded(Linear()))

	print_save("\nExtending the state space ($(Npn) iterations needed)")

	@sync @parallel for jpn in 1:Npn

		pnv = h.pngrid[jpn]
		
		N = size(h.Jgrid, 1)

		wage_pn, labor_pn, profits_pn = Array{Float64, 1}(N), Array{Float64, 1}(N), Array{Float64, 1}(N)
		for js in 1:N
			jw = h.Jgrid[js, 4]
			jζ = h.Jgrid[js, 5]
			jz = h.Jgrid[js, 6]

			wv = h.wgrid[jw]
			ζv = h.ζgrid[jζ]
			zv = h.zgrid[jz]

			jdef = (ζv != 1.0)

			labor_pn[js], wage_pn[js], profits_pn[js] = labor_market(h, jdef, zv, wv, pnv)
		end

		pC = price_index(h, pnv)
		pC_mat = ones(h.Nb, h.Nμ, h.Nσ, h.Nw, h.Nζ, h.Nz) * pC

		T_mat = govt_bc(h, wage_pn.*labor_pn) - reshape(profits_pn, h.Nb, h.Nμ, h.Nσ, h.Nw, h.Nζ, h.Nz)

		wL_mat  = reshape(wage_pn.*labor_pn, h.Nb, h.Nμ, h.Nσ, h.Nw, h.Nζ, h.Nz) * (1.0 - h.τ)

		# Re-solve for these values of wn and pn
		_, ϕa, ϕb, ϕc = opt_value(h, qʰ_mat, qᵍ_mat, wL_mat, T_mat, pC_mat, itp_R, itp_vf)
			
		ϕa_ext[:,:,:,:,:,:,:,:,jpn] = reshape(ϕa, h.Nω, h.Nϵ, h.Nb, h.Nμ, h.Nσ, h.Nw, h.Nζ, h.Nz)
		ϕb_ext[:,:,:,:,:,:,:,:,jpn] = reshape(ϕb, h.Nω, h.Nϵ, h.Nb, h.Nμ, h.Nσ, h.Nw, h.Nζ, h.Nz)
		ϕc_ext[:,:,:,:,:,:,:,:,jpn] = reshape(ϕc, h.Nω, h.Nϵ, h.Nb, h.Nμ, h.Nσ, h.Nw, h.Nζ, h.Nz)
	end

	h.ϕa_ext = ϕa_ext
	h.ϕb_ext = ϕb_ext
	h.ϕc_ext = ϕc_ext

	Void
end

transform_vars(m::Float64, cmin, cmax) = cmax - (cmax-cmin)/(1+exp(m))

function _unpack_origvars(x, xmin, xmax)
	y = zeros(x)
	for (jx, xv) in enumerate(x)
		y[jx] = transform_vars(xv, xmax[jx], xmin[jx])
	end
	return y
end

function labor_demand(h::Hank, w, z, pN; get_both::Bool = false)

	Ld_nontradables = (h.α_N * pN / w).^(1.0/(1.0-h.α_N))
	Ld_tradables    = (h.α_T * z  / w).^(1.0/(1.0-h.α_T))

	if get_both
		return Ld_nontradables, Ld_tradables
	else
		return Ld_nontradables + Ld_tradables
	end
end

function labor_market(h::Hank, jdef, zv, wv, pNv)
	""" Finds w and Lᵈ at the current state given a guess of pNv """
	TFP = ifelse(jdef, (1.0 - h.Δ) * zv, zv)
	w_constraint = h.γw * wv

	# Step 1: Assume w_t is at the constraint, find labor demand, and check whether the eq'm wage is above or below
	Ld_N, Ld_T = labor_demand(h, w_constraint, TFP, pNv; get_both=true)
	Ld = Ld_N + Ld_T

	# Step 2: If labor demand is lower than supply, find the wage above γw * wv that clears the labor mkt
	Ls = 1.0
	w_max = maximum(h.wgrid)
	if w_max < w_constraint
		print_save("\nSomething wrong with wages")
		w_max = w_constraint
	end

	w_new = w_constraint

	if Ld - Ls > 1e-4
		res = Optim.optimize(
			w -> (labor_demand(h, w, TFP, pNv) - Ls)^2,
				w_constraint, w_max, GoldenSection()
			)
		w_new = res.minimizer
		minf = Ls - labor_demand(h, w_new, TFP, pNv)
		abs(minf) > 1e-4? print_save("\nWARNING: Labor exc supply = $(@sprintf("%0.3g",minf)) at (w, ̄w) = ($(@sprintf("%0.3g",w_new)), $(@sprintf("%0.3g",w_max)))"): Void
		Ld_N, Ld_T = labor_demand(h, w_new, TFP, pNv; get_both=true)
		Ld = Ld_N + Ld_T
	end

	profits = pNv*zv .* Ld_N.^h.α_N + zv .* Ld_T.^h.α_T - w_new * (Ld_N + Ld_T)

	return Ld, w_new, profits
end


function mkt_clearing(h::Hank, itp_ϕc, Bpv, pNv, pNmin, pNmax, bv, μv, σv, wv, ζv, zv, jdefault; orig_vars::Bool = true, get_others::Bool = false)
	pN = pNv[1]
	if orig_vars == false
		pN = transform_vars(pN, pNmin, pNmax)
	end

	Ld, w_new, profits = labor_market(h, jdefault, zv, wv, pN)

	# Step 3: Get the household's policies at these prices

	val_A, val_B, val_C, sum_prob = 0., 0., 0., 0.
	for (jϵ, ϵv) in enumerate(h.ϵgrid)
		for jω = 1:length(h.ωgrid_fine)-1
			ωv  = h.ωgrid_fine[jω]
			ω1v = h.ωgrid_fine[jω+1]
			ωmv = 0.5*(ωv+ω1v)

			prob = pdf(LogNormal(μv, σv), ωmv-h.ωmin) * h.λϵ[jϵ] * (ω1v - ωv)

			ϕc = itp_ϕc[ωmv, ϵv, bv, μv, σv, wv, ζv, zv, pN]

			val_C += prob * ϕc
			
			sum_prob += prob
		end
	end

	# Step 4: Check market clearing for nontradables
	TFP = ifelse(jdefault, (1.0 - h.Δ) * zv, zv)
	Ld_N, _ = labor_demand(h, w_new, TFP, pN; get_both=true)
	supply_N = TFP * Ld_N^(h.α_N)

	demand = val_C / sum_prob

	""" Recover nontraded demand from total consumption """
	pC = price_index(h, pN)
	demand_N = demand * h.ϖ * (pN/pC)^(-h.η)

	F = supply_N - demand_N

	# """ Acá elegir qᵍ′ de manera consistente con μ′ y σ′ que a su vez son funciones de qᵍ′ """
	# ω′ = val_A + val_B * (h.κ + (1.0 - h.ρ) * qᵍ′)

	# μ′, σ′ = tomorrow_dist(h, val_A, val_B)

	# # Step 4: Compute NFI



	if get_others
		return w_new, Ld
	else
		return F
	end
end

function find_prices(h::Hank, itp_ϕc, Bpv, pNg, pNmin, pNmax, bv, μv, σv, wv, ζv, zv, jdefault)

	function wrap_mktclear!(pN::Vector, fvec=similar(x))

		out = mkt_clearing(h, itp_ϕc, Bpv, pN, pNmin, pNmax, bv, μv, σv, wv, ζv, zv, jdefault; orig_vars = false)

		fvec[:] = out
	end
	
	res = fsolve(wrap_mktclear!, [pNg])
	if res.:converged == false
		res2 = fsolve(wrap_mktclear!, [pNg], method=:lmdif)

		if res2.:converged || sum(res2.:f.^2) < sum(res.:f.^2)
			res = res2
		end
	end

	minf = res.:f[1]

	pN = transform_vars(res.:x[1], pNmin, pNmax)

	w, Ld = mkt_clearing(h, itp_ϕc, Bpv, pN, pNmin, pNmax, bv, μv, σv, wv, ζv, zv, jdefault; get_others=true)

	results = [w; pN; Ld]

	if abs(minf) > 1e-4
		print_save("\nNontradables exc supply = $(@sprintf("%0.3g",minf)) at pN = $(@sprintf("%0.3g",pN))")
	end

	return results, minf
end

function find_all_prices(h::Hank, itp_ϕc, B′_vec)

	N = size(h.Jgrid, 1)

	results = SharedArray{Float64}(N, 3)
	minf	= SharedArray{Float64}(N, 1)

	pN_guess = h.pN

	@sync @parallel for js in 1:N
		Bpv = B′_vec[js]
		pNg = pN_guess[js]

		jb = h.Jgrid[js, 1]
		jμ = h.Jgrid[js, 2]
		jσ = h.Jgrid[js, 3]
		jw = h.Jgrid[js, 4]
		jζ = h.Jgrid[js, 5]
		jz = h.Jgrid[js, 6]

		bv = h.bgrid[jb]
		μv = h.μgrid[jμ]
		σv = h.σgrid[jσ]
		wv = h.wgrid[jw]
		ζv = h.ζgrid[jζ]
		zv = h.zgrid[jz]		

		jdefault = (ζv != 1.0)

		pNmin, pNmax = minimum(h.pngrid), maximum(h.pngrid)

		results[js, :], minf[js, :] = find_prices(h, itp_ϕc, Bpv, pNg, pNmin, pNmax, bv, μv, σv, wv, ζv, zv, jdefault)
	end
		
	
	return results, minf
end

function update_state_functions!(h::Hank, upd_η::Float64)
	all_knots = (h.ωgrid, h.ϵgrid, h.bgrid, h.μgrid, h.σgrid, h.wgrid, h.ζgrid, h.zgrid, h.pngrid)
	# agg_knots = (h.bgrid, h.μgrid, h.σgrid, h.wgrid, h.ζgrid, h.zgrid)

	itp_ϕa  = interpolate(all_knots, h.ϕa_ext, Gridded(Linear()))
	itp_ϕb  = interpolate(all_knots, h.ϕb_ext, Gridded(Linear()))
	itp_ϕc  = interpolate(all_knots, h.ϕc_ext, Gridded(Linear()))

	results, minf = find_all_prices(h, itp_ϕc, h.issuance)

	dist = Array{Float64,1}(3)
	dist[1] = sqrt.(sum( (results[:, 1] - h.wage).^2 )) / sqrt.(sum(h.wage.^2))
	dist[2] = sqrt.(sum( (results[:, 2] - h.pN).^2 ))   / sqrt.(sum(h.pN.^2))
	dist[3] = sqrt.(sum( (results[:, 3] - h.Ld).^2 ))   / sqrt.(sum(h.Ld.^2))

	h.wage 	= upd_η * results[:, 1] + (1.0-upd_η) * h.wage
	h.pN 	= upd_η * results[:, 2] + (1.0-upd_η) * h.pN
	h.Ld 	= upd_η * results[:, 3] + (1.0-upd_η) * h.Ld

	h.w′	= h.wage

	meanf = zeros(size(minf)[end])
	for jf in 1:size(minf)[end]
		meanf[jf] = mean(minf[:,:,:,:,jf])
	end
	return meanf, dist
end

function find_q(h::Hank, q, a, b, var_a, var_b, cov_ab, Bpv, wpv, thres, zpv, jdef, itp_qᵍ, reentry; get_μσ::Bool=false)

	ζpv = 1
	if jdef && reentry==false
		ζpv = 3
	end
	if jdef == false && zpv <= thres
		ζpv = 2
	end
	
	R = (ζpv==1) * h.κ + (1.0 - h.ℏ * (ζpv==2)) .* ((1.0-h.ρ)*q)

	Eω   = a + R*b
	varω = var_a + R^2 * var_b + 2*R * cov_ab

	# print_save("\nEω, varω = $Eω, $varω")
	Eσ2 = 1.0 + varω / ( (Eω - h.ωmin)^2 )
	if Eσ2 < 0
		print_save("\n1 + vω / (Eω-ωmin)² = $( Eσ2)")
	end

	σ2 = log( Eσ2 )

	μpv = log(Eω - h.ωmin) - 0.5 * σ2
	σpv = sqrt(σ2)

	new_q = itp_qᵍ[Bpv, μpv, σpv, wpv, ζpv, zpv]

	if get_μσ
		return μpv, σpv
	else
		return new_q
	end
end


function compute_stats_logN(h::Hank, js, a, b, var_a, var_b, cov_ab, itp_qᵍ, Bpv, wpv, thres)

	ζv = h.ζgrid[h.Jgrid[js, 5]]
	jdef = (ζv != 1.0)

	μ, σ = Array{Float64, 2}(h.Nz, 2), Array{Float64, 2}(h.Nz, 2)

	for (jzp, zpv) in enumerate(h.zgrid)
		reentry = true
		qmin, qmax = minimum(h.qᵍ), maximum(h.qᵍ)

		res = Optim.optimize(
			q -> (find_q(h, q, a, b, var_a, var_b, cov_ab, Bpv, wpv, thres, zpv, jdef, itp_qᵍ, reentry) - q)^2,
			qmin, qmax, GoldenSection()
			)
		qᵍ = res.minimizer
		res.minimum > 1e-4? print_save("WARNING: Can't find consistent qᵍ"): Void

		μ[jzp, 1], σ[jzp, 1] = find_q(h, qᵍ, a, b, var_a, var_b, cov_ab, Bpv, wpv, thres, zpv, jdef, itp_qᵍ, reentry; get_μσ = true)

		if jdef
			reentry = false
			res = Optim.optimize(
				q -> (find_q(h, q, a, b, var_a, var_b, cov_ab, Bpv, wpv, thres, zpv, jdef, itp_qᵍ, reentry) - q)^2,
				qmin, qmax, GoldenSection()
				)
			qᵍ = res.minimizer
			res.minimum > 1e-4? print_save("WARNING: Can't find consistent qᵍ"): Void

			μ[jzp, 2], σ[jzp, 2] = find_q(h, qᵍ, a, b, var_a, var_b, cov_ab, Bpv, wpv, thres, zpv, jdef, itp_qᵍ, reentry; get_μσ = true)
		else
			μ[jzp, 2], σ[jzp, 2] = μ[jzp, 1], σ[jzp, 1]
		end
	end
	return μ, σ
end

function new_expectations(h::Hank, itp_ϕa, itp_ϕb, itp_qᵍ, Bpv, wpv, thres, js, jdef)

	bv = h.bgrid[h.Jgrid[js, 1]]
	μv = h.μgrid[h.Jgrid[js, 2]]
	σv = h.σgrid[h.Jgrid[js, 3]]
	wv = h.wgrid[h.Jgrid[js, 4]]
	ζv = h.ζgrid[h.Jgrid[js, 5]]
	zv = h.zgrid[h.Jgrid[js, 6]]
	
	val_a, val_b, val_a2, val_b2, val_ab, sum_prob = 0., 0., 0., 0., 0., 0.
	for (jϵ, ϵv) in enumerate(h.ϵgrid)
		for jω = 1:length(h.ωgrid_fine)-1
			ωv  = h.ωgrid_fine[jω]
			ω1v = h.ωgrid_fine[jω+1]
			ωmv = 0.5*(ωv+ω1v)

			prob = pdf(LogNormal(μv, σv), ωmv-h.ωmin) * h.λϵ[jϵ] * (ω1v - ωv)

			ϕa = itp_ϕa[ωmv, ϵv, bv, μv, σv, wv, ζv, zv]
			ϕb = itp_ϕb[ωmv, ϵv, bv, μv, σv, wv, ζv, zv]

			val_a  += prob * ϕa
			val_a2 += prob * ϕa^2
			val_b  += prob * ϕb
			val_b2 += prob * ϕb^2
			val_ab += prob * ϕa * ϕb
			
			sum_prob += prob
		end
	end

	a  = val_a  / sum_prob
	a2 = val_a2 / sum_prob
	b  = val_b  / sum_prob
	b2 = val_b2 / sum_prob
	ab = val_ab / sum_prob

	var_a  = a2 - a^2
	var_b  = b2 - b^2
	cov_ab = ab - a*b

	# print_save("\nVa, Vb, cov = $var_a, $var_b, $cov_ab")

	μ′, σ′ = compute_stats_logN(h, js, a, b, var_a, var_b, cov_ab, itp_qᵍ, Bpv, wpv, thres)

	return μ′, σ′
end


function find_all_expectations(h::Hank, itp_ϕa, itp_ϕb, itp_qᵍ, B′_vec, w′_vec, thres_vec)
	N = size(h.Jgrid, 1)

	μ′ = SharedArray{Float64}(N, h.Nz, 2)
	σ′ = SharedArray{Float64}(N, h.Nz, 2)

	@sync @parallel for js in 1:N
		Bpv = B′_vec[js]
		wpv = w′_vec[js]
		thres = thres_vec[js]

		jζ = h.Jgrid[js, 5]
		jdefault = (jζ != 1.0)

		μ′[js, :, :], σ′[js, :, :] = new_expectations(h, itp_ϕa, itp_ϕb, itp_qᵍ, Bpv, wpv, thres, js, jdefault)
	end
		
	
	return μ′, σ′
end

function update_expectations!(h::Hank, upd_η::Float64)
	""" 
	Computes mean and variance of tomorrow's distribution and deduces parameters for logN
	"""

	qᵍmt = reshape(h.qᵍ, h.Nb, h.Nμ, h.Nσ, h.Nw, h.Nζ, h.Nz)

	all_knots = (h.ωgrid, h.ϵgrid, h.bgrid, h.μgrid, h.σgrid, h.wgrid, h.ζgrid, h.zgrid)
	agg_knots = (h.bgrid, h.μgrid, h.σgrid, h.wgrid, h.ζgrid, h.zgrid)

	itp_ϕa = interpolate(all_knots, h.ϕa, Gridded(Linear()))
	itp_ϕb = interpolate(all_knots, h.ϕb, Gridded(Linear()))
	itp_qᵍ = interpolate(agg_knots, qᵍmt, Gridded(Linear()))

	μ′_new, σ′_new = find_all_expectations(h, itp_ϕa, itp_ϕb, itp_qᵍ, h.issuance, h.w′, h.def_thres)
	# μ′_new, σ′_new = h.μ′, h.σ′

	print_save("\nOld μ min, mean, max = $(minimum(h.μ′)), $(mean(h.μ′)), $(maximum(h.μ′))")
	print_save("\nOld σ min, mean, max = $(minimum(h.σ′)), $(mean(h.σ′)), $(maximum(h.σ′))")
	h.μ′ = upd_η * μ′_new + (1.0 - upd_η) * h.μ′
	h.σ′ = upd_η * σ′_new + (1.0 - upd_η) * h.σ′

	h.μ′ = max(min(h.μ′, maximum(h.μgrid)), minimum(h.μgrid))
	h.σ′ = max(min(h.σ′, maximum(h.σgrid)), minimum(h.σgrid))

	print_save("\nNew μ min, mean, max = $(minimum(h.μ′)), $(mean(h.μ′)), $(maximum(h.μ′))")
	print_save("\nNew σ min, mean, max = $(minimum(h.σ′)), $(mean(h.σ′)), $(maximum(h.σ′))")

	Void
end
