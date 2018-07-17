TFP_N(z, Δ, ζ) = 1.0   # * (1.0 - Δ*(ζ==2))
TFP_T(z, Δ, ζ) = exp(z) * (1.0 - Δ*(ζ==2))

function extend_state_space!(h::Hank, qʰ_mat, qᵍ_mat, T_mat)

	Npn = length(h.pngrid)

	ϕa_ext = Array{Float64}(h.Nω, h.Nϵ, h.Nb, h.Nμ, h.Nσ, h.Nw, h.Nζ, h.Nz, Npn)
	ϕb_ext = Array{Float64}(h.Nω, h.Nϵ, h.Nb, h.Nμ, h.Nσ, h.Nw, h.Nζ, h.Nz, Npn)
	ϕc_ext = Array{Float64}(h.Nω, h.Nϵ, h.Nb, h.Nμ, h.Nσ, h.Nw, h.Nζ, h.Nz, Npn)

	itp_vf = make_itp(h, h.vf; agg=false)
	itp_qᵍ = make_itp(h, h.qᵍ; agg=true)

	print_save("\nExtending the state space ($(Npn) iterations needed)")

	for jpn in 1:Npn

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

			labor_pn[js], wage_pn[js], profits_pn[js], _ = labor_market(h, ζv, zv, wv, pnv)
		end

		pC = price_index(h, pnv)
		pC_mat = ones(h.Nb, h.Nμ, h.Nσ, h.Nw, h.Nζ, h.Nz) * pC

		T_mat = govt_bc(h, wage_pn.*labor_pn)# - reshape(profits_pn - h.profits, h.Nb, h.Nμ, h.Nσ, h.Nw, h.Nζ, h.Nz)
		Π_mat = reshape(profits_pn, h.Nb, h.Nμ, h.Nσ, h.Nw, h.Nζ, h.Nz)

		wL_mat  = reshape(wage_pn.*labor_pn, h.Nb, h.Nμ, h.Nσ, h.Nw, h.Nζ, h.Nz) * (1.0 - h.τ)

		# Re-solve for these values of wn and pn
		_, ϕa, ϕb, ϕe, ϕc = opt_value(h, qʰ_mat, qᵍ_mat, wL_mat, T_mat, pC_mat, Π_mat, itp_qᵍ, itp_vf; resolve = true, verbose = false)

		isapprox(sum(abs.(ϕc)), 0)? print_save("\nWARNING: ϕc(pN = $(round(pnv, 2))) ≡ 0 when extending state space"): Void

		for jz in 1:h.Nz, jζ in 1:h.Nζ, jw in 1:h.Nw, jσ in 1:h.Nσ, jμ in 1:h.Nμ, jb in 1:h.Nb, jϵ in 1:h.Nϵ, jω in 1:h.Nω
			ϕa_ext[jω,jϵ,jb,jμ,jσ,jw,jζ,jz,jpn] = ϕa[jω,jϵ,jb,jμ,jσ,jw,jζ,jz]
			ϕb_ext[jω,jϵ,jb,jμ,jσ,jw,jζ,jz,jpn] = ϕb[jω,jϵ,jb,jμ,jσ,jw,jζ,jz]
			ϕc_ext[jω,jϵ,jb,jμ,jσ,jw,jζ,jz,jpn] = ϕc[jω,jϵ,jb,jμ,jσ,jw,jζ,jz]
		end
	end

	!isnan(sum(ϕa_ext)) || print_save("ERROR: $(isnan(sum(ϕa_ext))) NaN counts in ϕa_ext")
	!isnan(sum(ϕa_ext)) || throw(error("$(isnan(sum(ϕa_ext))) NaN counts in ϕa_ext"))

	isapprox(sum(abs.(ϕc_ext)), 0)? print_save("\nWARNING: ϕc ≡ 0 when extending state space"): Void

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

function labor_demand(h::Hank, w, z, ζ, pN; get_both::Bool = false)

	Δ = h.Δ

	Ld_nontradables = (h.α_N * pN * TFP_N(z,Δ,ζ) / w).^(1.0/(1.0-h.α_N))
	Ld_tradables    = (h.α_T * 1. * TFP_T(z,Δ,ζ) / w).^(1.0/(1.0-h.α_T))

	if get_both
		return Ld_nontradables, Ld_tradables
	else
		return Ld_nontradables + Ld_tradables
	end
end

function labor_market(h::Hank, ζv, zv, wv, pNv)
	""" Finds w and Lᵈ at the current state given a guess of pNv """
	w_constraint = h.γw * wv

	# Step 1: Assume wₜ is at the constraint, find labor demand, and check whether the eq'm wage is above or below
	Ld_N, Ld_T = labor_demand(h, w_constraint, zv, ζv, pNv; get_both=true)
	Ld = Ld_N + Ld_T

	# Step 2: If labor demand is lower than supply, find the wage above γw * wv that clears the labor mkt
	Ls = 1.0
	w_max = maximum(h.wgrid)
	if w_max < w_constraint
		# print_save("\nSomething wrong with wages")
		w_max = w_constraint
	end

	w_new = w_constraint

	if Ld - Ls > 1e-4
		res = Optim.optimize(
			w -> (labor_demand(h, w, zv, ζv, pNv) - Ls)^2,
				w_constraint, 1.*w_max, GoldenSection()
			)
		w_new = res.minimizer
		minf = Ls - labor_demand(h, w_new, zv, ζv, pNv)
		abs(minf) < 1e-4 || print_save("\nWARNING: Labor exc supply = $(@sprintf("%0.3g",minf)) at (w, w_max) = ($(@sprintf("%0.3g",w_new)), $(@sprintf("%0.3g",w_max))) at (z,ζ,pN) = $([zv, ζv, pNv])")
		Ld_N, Ld_T = labor_demand(h, w_new, zv, ζv, pNv; get_both=true)
		Ld = Ld_N + Ld_T
	end

	output 	= pNv * TFP_N(zv,h.Δ,ζv) * Ld_N^h.α_N + TFP_T(zv,h.Δ,ζv) * Ld_T^h.α_T
	profits = output - w_new * (Ld_N + Ld_T)

	return Ld, w_new, profits, output
end


function mkt_clearing(h::Hank, itp_ϕc, G, Bpv, pNv, pNmin, pNmax, bv, μv, σv, wv, jζ, jz, jdefault; orig_vars::Bool = true, get_others::Bool = false, get_both::Bool=false)
	typeof(pNv) == Vector{Float64}?	pN = pNv[1]: pN = pNv
	if orig_vars == false
		pN = transform_vars(pN, pNmin, pNmax)
	end

	isnan(pN)? print_save("\nWARNING: pNv, pN = $(pNv[1]), $(pN)"): Void

	ζv, zv = h.ζgrid[jζ], h.zgrid[jz]

	Ld, w_new, profits, output = labor_market(h, ζv, zv, wv, pN)

	# Check market clearing for nontradables
	Ld_N, _  = labor_demand(h, w_new, zv, ζv, pN; get_both=true)
	supply_N = TFP_N(zv, h.Δ, ζv) * Ld_N^(h.α_N)

	# Get the household's policies at these prices
	ωmin_int, ωmax_int = quantile.(LogNormal(μv, σv), [.005; .995]) + h.ωmin
	val_int_C = 0.
	for (jϵ, ϵv) in enumerate(h.ϵgrid)
		f(ω) = pdf(LogNormal(μv, σv), ω-h.ωmin) * h.λϵ[jϵ] * itp_ϕc[ω, jϵ, bv, μv, σv, wv, jζ, jz, pN]
	
		(val, err) = hquadrature(f, ωmin_int, ωmax_int, reltol=1e-12, abstol=0, maxevals=0)
	
		val_int_C += val
	end
	# val_C, sum_prob = 0., 0.
	# for (jϵ, ϵv) in enumerate(h.ϵgrid)
	# 	for jω = 1:length(h.ωgrid_fine)-1
	# 		ωv  = h.ωgrid_fine[jω]
	# 		ω1v = h.ωgrid_fine[jω+1]
	# 		ωmv = 0.5*(ωv+ω1v)

	# 		prob = pdf(LogNormal(μv, σv), ωmv-h.ωmin) * h.λϵ[jϵ] * (ω1v - ωv)

	# 		ϕc = itp_ϕc[ωmv, jϵ, bv, μv, σv, wv, jζ, jz, pN]

	# 		val_C  += prob * ϕc
	# 		sum_prob += prob
	# 	end
	# end
	# val_int_C = val_C / sum_prob
	# isapprox(sum_prob, 1) || print_save("\nWARNING: Something wrong computing aggregate consumption")

 	# Recover nontraded demand from total consumption
	pC = price_index(h, pN)
	demand_N_cons = val_int_C * h.ϖ * (pN/pC)^(-h.η)
	demand_N_govt = G / pN

	demand_N = demand_N_cons + demand_N_govt

	F = supply_N - demand_N

	if get_others
		return w_new, Ld, output
	elseif get_both
		return supply_N, demand_N
	else
		return F
	end
end

function find_prices(h::Hank, itp_ϕc, G, Bpv, pNg, pNmin, pNmax, bv, μv, σv, wv, jζ, jz, jdefault)

	# First find eq'm assuming the constraint does not bind
	w_slack = 0.5 * h.wgrid[1]
	res = Optim.optimize(
		pN -> mkt_clearing(h, itp_ϕc, G, Bpv, pN, pNmin, pNmax, bv, μv, σv, w_slack, jζ, jz, jdefault; orig_vars = true)^2,
		pNmin, pNmax, GoldenSection()
		)
	pN = res.minimizer
	w, Ld, output = mkt_clearing(h, itp_ϕc, G, Bpv, pN, pNmin, pNmax, bv, μv, σv, w_slack, jζ, jz, jdefault; orig_vars=true, get_others=true)

	if w >= h.γw * wv && res.minimum < 1e-6
		pN > pNmax? exc_dem = 1: exc_dem = 0
		pN < pNmin? exc_sup = 1: exc_sup = 0

		minf = mkt_clearing(h, itp_ϕc, G, Bpv, pN, pNmin, pNmax, bv, μv, σv, wv, jζ, jz, jdefault; orig_vars = true)

		results = [w; pN; Ld; output]

		return results, minf, exc_dem, exc_sup
	end

	res = Optim.optimize(
		pN -> mkt_clearing(h, itp_ϕc, G, Bpv, pN, pNmin, pNmax, bv, μv, σv, wv, jζ, jz, jdefault; orig_vars = true)^2,
		0.9*pNmin, 1.1*pNmax, GoldenSection()
		)
	pN = res.minimizer
	minf = mkt_clearing(h, itp_ϕc, G, Bpv, pN, pNmin, pNmax, bv, μv, σv, wv, jζ, jz, jdefault; orig_vars = true)

	pN > pNmax? exc_dem = 1: exc_dem = 0
	pN < pNmin? exc_sup = 1: exc_sup = 0
	# if res.minimum > 1e-6
	# 	exc_dem, exc_sup = 1, 1
	# end

	if false # Deprecated fsolve-based method
		function wrap_mktclear!(pN::Vector, fvec=similar(x))

			out = mkt_clearing(h, itp_ϕc, G, Bpv, pN, pNmin, pNmax, bv, μv, σv, wv, jζ, jz, jdefault; orig_vars = false)

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
	end

	w, Ld, output = mkt_clearing(h, itp_ϕc, G, Bpv, pN, pNmin, pNmax, bv, μv, σv, wv, jζ, jz, jdefault; get_others=true)

	results = [w; pN; Ld; output]

	if abs(minf) > 1e-4
		# print_save("\nNontradables exc supply = $(@sprintf("%0.4g",minf)) at pN = $(@sprintf("%0.4g",pN))")
	end

	return results, minf, exc_dem, exc_sup
end

function find_all_prices(h::Hank, itp_ϕc, B′_vec, G_vec)

	N = size(h.Jgrid, 1)

	results = SharedArray{Float64}(N, 4)
	minf	= SharedArray{Float64}(N, 1)
	exc_dem = SharedArray{Float64}(N)
	exc_sup = SharedArray{Float64}(N)

	pN_guess = h.pN

	@sync @parallel for js in 1:N
		Bpv = B′_vec[js]
		G = G_vec[js]
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

		results[js, :], minf[js, :], exc_dem[js], exc_sup[js] = find_prices(h, itp_ϕc, G, Bpv, pNg, pNmin, pNmax, bv, μv, σv, wv, jζ, jz, jdefault)
	end

	return results, minf, exc_dem, exc_sup
end

function update_state_functions!(h::Hank, upd_η::Float64)
	itp_ϕc = make_itp(h, h.ϕc_ext; agg=false)

	results, minf, exc_dem, exc_sup = find_all_prices(h, itp_ϕc, h.issuance, h.spending)

	dist = Array{Float64,1}(3)
	dist[1] = sqrt.(sum( (results[:, 1] - h.wage).^2 )) / sqrt.(sum(h.wage.^2))
	dist[2] = sqrt.(sum( (results[:, 2] - h.pN).^2 ))   / sqrt.(sum(h.pN.^2))
	dist[3] = sqrt.(sum( (results[:, 3] - h.Ld).^2 ))   / sqrt.(sum(h.Ld.^2))

	h.pN = upd_η * results[:, 2] + (1.0-upd_η) * h.pN

	consistent_others = true

	if consistent_others
		N = size(h.Jgrid,1)
		wage, Ld, output = zeros(N), zeros(N), zeros(N)
		for js in 1:N
			Bpv = h.issuance[js]
			G = h.spending[js]

			pN = h.pN[js]
			pNmin, pNmax = minimum(h.pngrid), maximum(h.pngrid)

			bv = h.bgrid[h.Jgrid[js, 1]]
			μv = h.μgrid[h.Jgrid[js, 2]]
			σv = h.σgrid[h.Jgrid[js, 3]]
			wv = h.wgrid[h.Jgrid[js, 4]]
			jζ = h.Jgrid[js, 5]
			jz = h.Jgrid[js, 6]

			wage[js], Ld[js], output[js] = mkt_clearing(h, itp_ϕc, G, Bpv, pN, pNmin, pNmax, bv, μv, σv, wv, jζ, jz, (jζ!=1); get_others=true)
		end
		h.wage, h.Ld, h.output = wage, Ld, output
	else
		h.wage 	 = upd_η * results[:, 1] + (1.0-upd_η) * h.wage
		h.Ld 	 = upd_η * results[:, 3] + (1.0-upd_η) * h.Ld
		h.output = upd_η * results[:, 4] + (1.0-upd_η) * h.output
	end

	h.profits = h.output - h.wage .* h.Ld
	h.w′ = h.wage
	mean_f = mean(minf)
	max_f = minf[indmax(abs.(minf))]

	exc_dem_prop = sum(exc_dem) / length(exc_dem)
	exc_sup_prop = sum(exc_sup) / length(exc_sup)

	return exc_dem_prop, exc_sup_prop, mean_f, max_f, dist
end

function update_grids_pw!(h::Hank, exc_dem_prop, exc_sup_prop)

	pN_down = minimum(h.pngrid)
	if exc_sup_prop > 0.01
		pN_down = pN_down * 0.95
	elseif exc_sup_prop == 0.
		pN_down = pN_down * 1.01
	end
	pN_up = maximum(h.pngrid)
	if exc_dem_prop > 0.01
		pN_up = pN_up * 1.05
	elseif exc_dem_prop == 0.
		pN_up = pN_up * 0.99
	end

	Ls = 1.0
	res1 = Optim.optimize(
			w -> (labor_demand(h, w, exp(h.zgrid[end]), pN_up) - Ls).^2,
			h.wgrid[1], h.wgrid[end] * 2.0, GoldenSection()
			)
	res2 = Optim.optimize(
			w -> (labor_demand(h, w, exp(h.zgrid[end]), pN_down) - Ls).^2,
			h.wgrid[1], h.wgrid[end] * 2.0, GoldenSection()
			)
	w_up = max(res1.minimizer, res2.minimizer) * 1.25
	res1 = Optim.optimize(
			w -> (labor_demand(h, w, (1.0-h.Δ) * exp(h.zgrid[1]), pN_up) - Ls).^2,
			0.5 * h.wgrid[1], h.wgrid[end], GoldenSection()
			)
	res2 = Optim.optimize(
			w -> (labor_demand(h, w, (1.0-h.Δ) * exp(h.zgrid[1]), pN_down) - Ls).^2,
			0.5 * h.wgrid[1], h.wgrid[end], GoldenSection()
			)
	w_down = min(res1.minimizer, res2.minimizer) * 0.75

	h.pngrid = collect(linspace(pN_down, pN_up, length(h.pngrid)))
	new_wgrid = collect(linspace(w_down, w_up, h.Nw))

	return new_wgrid
end


function find_q(h::Hank, q, a, b, var_a, var_b, cov_ab, Bpv, wpv, exp_rep, jzp, jdef, itp_qᵍ, reentry; get_μσ::Bool=false)

	zpv = h.zgrid[jzp]

	ζpv = 1
	haircut = 0.0
	if jdef && reentry==false
		ζpv = 2
		haircut = 0.0
	end
	if jdef == false && exp_rep[jzp] < 0.5
		ζpv = 2
		haircut = h.ℏ
	end

	R = (ζpv==1) * h.κ + (1.0 - haircut) .* ((1.0-h.ρ)*q)

	Eω   = a + R*b
	varω = var_a + R^2 * var_b + 2*R * cov_ab

	varω > 0. || print_save("\nvar_a, var_b, cov_ab = $(var_a), $(var_b), $(cov_ab)")

	Eσ2 = 1.0 + varω / ( (Eω - h.ωmin)^2 )
	if isapprox(Eσ2, 1.)
		Eσ2 = min(Eσ2, 1.)
	end

	Eσ2 >= 1. || print_save("\n1 + vω / (Eω-ωmin)² = $(Eσ2)")

	σ2 = log( Eσ2 )

	μpv = log(Eω - h.ωmin) - 0.5 * σ2
	σpv = sqrt(σ2)

	new_q = itp_qᵍ[(1.0 - haircut) .* Bpv, μpv, σpv, wpv, ζpv, jzp]

	if get_μσ
		return μpv, σpv
	else
		return new_q
	end
end


function compute_stats_logN(h::Hank, js, a, b, var_a, var_b, cov_ab, itp_qᵍ, Bpv, wpv, exp_rep)

	ζv = h.ζgrid[h.Jgrid[js, 5]]
	jdef = (ζv != 1.0)


	μ, σ, q = Array{Float64, 2}(h.Nz, 2), Array{Float64, 2}(h.Nz, 2), Array{Float64, 2}(h.Nz, 2)
	alarm_mat = Array{Float64, 2}(h.Nz, 2)

	for (jzp, zpv) in enumerate(h.zgrid)
		reentry = true
		qmin, qmax = minimum(h.qᵍ), maximum(h.qᵍ)

		res = Optim.optimize(
			q -> (find_q(h, q, a, b, var_a, var_b, cov_ab, Bpv, wpv, exp_rep, jzp, jdef, itp_qᵍ, reentry) - q)^2,
			qmin, qmax, GoldenSection()
			)
		q[jzp, 1] = res.minimizer
		res.minimum > 1e-4? alarm_mat[jzp, 1] = 1: alarm_mat[jzp, 1] = 0

		μ[jzp, 1], σ[jzp, 1] = find_q(h, q[jzp, 1], a, b, var_a, var_b, cov_ab, Bpv, wpv, exp_rep, jzp, jdef, itp_qᵍ, reentry; get_μσ = true)

		if jdef
			reentry = false
			res = Optim.optimize(
				q -> (find_q(h, q, a, b, var_a, var_b, cov_ab, Bpv, wpv, exp_rep, jzp, jdef, itp_qᵍ, reentry) - q)^2,
				qmin, qmax, GoldenSection()
				)
			q[jzp, 2] = res.minimizer
			res.minimum > 1e-4? alarm_mat[jzp, 2] = 1: alarm_mat[jzp, 2] = 0

			μ[jzp, 2], σ[jzp, 2] = find_q(h, q[jzp, 2], a, b, var_a, var_b, cov_ab, Bpv, wpv, exp_rep, jzp, jdef, itp_qᵍ, reentry; get_μσ = true)
		else
			μ[jzp, 2], σ[jzp, 2], q[jzp, 2], alarm_mat[jzp, 2] = μ[jzp, 1], σ[jzp, 1], q[jzp, 1], alarm_mat[jzp, 1]
		end
	end
	return μ, σ, q, alarm_mat
end

function new_expectations(h::Hank, itp_ϕa, itp_ϕb, itp_qᵍ, Bpv, wpv, exp_rep, js, jdef)

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

	val_a, val_b, val_a2, val_b2, val_ab, sum_prob = 0., 0., 0., 0., 0., 0.

	ωmin_int, ωmax_int = quantile.(LogNormal(μv, σv), [.005; .995]) + h.ωmin
	ωmax_int = max(ωmax_int, h.ωmax)
	for (jϵ, ϵv) in enumerate(h.ϵgrid)
		fA(ω) = pdf(LogNormal(μv, σv), ω-h.ωmin) * max(h.ωmin, itp_ϕa[ω, jϵ, bv, μv, σv, wv, jζ, jz])
		(valA, err) = hquadrature(fA, ωmin_int, ωmax_int, reltol=1e-10, abstol=1e-8, maxevals=0)
		val_a += valA * h.λϵ[jϵ]
		fA2(ω) = pdf(LogNormal(μv, σv), ω-h.ωmin) * max(h.ωmin, itp_ϕa[ω, jϵ, bv, μv, σv, wv, jζ, jz])^2
		(valA2, err) = hquadrature(fA2, ωmin_int, ωmax_int, reltol=1e-10, abstol=1e-8, maxevals=0)
		val_a2 += valA2 * h.λϵ[jϵ]
		fB(ω) = pdf(LogNormal(μv, σv), ω-h.ωmin) * max(0., itp_ϕb[ω, jϵ, bv, μv, σv, wv, jζ, jz])
		(valB, err) = hquadrature(fB, ωmin_int, ωmax_int, reltol=1e-10, abstol=1e-8, maxevals=0)
		val_b += valB * h.λϵ[jϵ]
		fB2(ω) = pdf(LogNormal(μv, σv), ω-h.ωmin) * max(0., itp_ϕb[ω, jϵ, bv, μv, σv, wv, jζ, jz])^2
		(valB2, err) = hquadrature(fB2, ωmin_int, ωmax_int, reltol=1e-10, abstol=1e-8, maxevals=0)
		val_b2 += valB2 * h.λϵ[jϵ]
		fAB(ω) = pdf(LogNormal(μv, σv), ω-h.ωmin) * max(h.ωmin, itp_ϕa[ω, jϵ, bv, μv, σv, wv, jζ, jz]) * max(0., itp_ϕb[ω, jϵ, bv, μv, σv, wv, jζ, jz])
		(valAB, err) = hquadrature(fAB, ωmin_int, ωmax_int, reltol=1e-10, abstol=1e-8, maxevals=0)
		val_ab += valAB * h.λϵ[jϵ]
	end
	sum_prob = 1.
	# for (jϵ, ϵv) in enumerate(h.ϵgrid)
	# 	for jω = 1:length(h.ωgrid_fine)-1
	# 		ωv  = h.ωgrid_fine[jω]
	# 		ω1v = h.ωgrid_fine[jω+1]
	# 		ωmv = 0.5*(ωv+ω1v)

	# 		prob = pdf(LogNormal(μv, σv), ωmv-h.ωmin) * h.λϵ[jϵ] * (ω1v - ωv)

	# 		ϕa = max(itp_ϕa[ωmv, jϵ, bv, μv, σv, wv, jζ, jz], h.ωmin)
	# 		ϕb = max(itp_ϕb[ωmv, jϵ, bv, μv, σv, wv, jζ, jz], 0.)

	# 		val_a  += prob * ϕa
	# 		val_a2 += prob * ϕa^2
	# 		val_b  += prob * ϕb
	# 		val_b2 += prob * ϕb^2
	# 		val_ab += prob * ϕa * ϕb

	# 		sum_prob += prob
	# 	end
	# end

	!isnan(val_a+val_a2+val_b+val_b2+val_ab) || print_save("\na,a2,b,b2,ab = $([val_a,val_a2,val_b,val_b2,val_ab])")
	!isapprox(sum_prob, 0.) || print_save("\nsum_prob = $(sum_prob)")

	a  = val_a  / sum_prob
	a2 = val_a2 / sum_prob
	b  = val_b  / sum_prob
	b2 = val_b2 / sum_prob
	ab = val_ab / sum_prob

	var_a  = a2 - a^2
	var_b  = b2 - b^2
	cov_ab = ab - a*b

	!isnan(var_a+var_b+cov_ab) || print_save("\nVa, Vb, cov = $var_a, $var_b, $cov_ab at $([jb, jμ, jσ, jw, jζ, jz])")

	μ′, σ′, q, alarm_vec = compute_stats_logN(h, js, a, b, var_a, var_b, cov_ab, itp_qᵍ, Bpv, wpv, exp_rep)

	return μ′, σ′, alarm_vec
end


function find_all_expectations(h::Hank, itp_ϕa, itp_ϕb, itp_qᵍ, B′_vec, w′_vec, thres_vec)
	N = size(h.Jgrid, 1)

	μ′ = SharedArray{Float64}(N, h.Nz, 2)
	σ′ = SharedArray{Float64}(N, h.Nz, 2)
	alarm_vec = SharedArray{Float64}(N, h.Nz, 2)

	@sync @parallel for js in 1:N
		Bpv = B′_vec[js]
		wpv = w′_vec[js]
		thres = thres_vec[js]

		# js % 20 != 0 || print_save("\n$js")

		jb = h.Jgrid[js, 1]
		jμ = h.Jgrid[js, 2]
		jσ = h.Jgrid[js, 3]
		jw = h.Jgrid[js, 4]
		jζ = h.Jgrid[js, 5]
		jz = h.Jgrid[js, 6]

		rep_mat = reshape(h.repay, h.Nb, h.Nμ, h.Nσ, h.Nw, h.Nζ, h.Nz, h.Nz)
		exp_rep = rep_mat[jb, jμ, jσ, jw, jζ, jz, :]


		jζ = h.Jgrid[js, 5]
		jdefault = (jζ != 1.0)

		μ′[js, :, :], σ′[js, :, :], alarm_vec[js, :, :] = new_expectations(h, itp_ϕa, itp_ϕb, itp_qᵍ, Bpv, wpv, exp_rep, js, jdefault)
	end

	print_save("\n")
	if sum(alarm_vec[:]) >= 1
		print_save("WARNING: ")
	end
	print_save("Couldn't find qᵍ $(round(100*sum(alarm_vec[:])/N,0))% of the time")

	return μ′, σ′
end

function update_expectations!(h::Hank, upd_η::Float64)
	"""
	Computes mean and variance of tomorrow's distribution and deduces parameters for logN
	"""

	μ′_old = copy(h.μ′)
	σ′_old = copy(h.σ′)

	dist_exp = Array{Float64,1}(2)

	itp_ϕa = make_itp(h, h.ϕa; agg=false)
	itp_ϕb = make_itp(h, h.ϕb; agg=false)
	itp_qᵍ = make_itp(h, h.qᵍ; agg=true)

	μ′_new, σ′_new = find_all_expectations(h, itp_ϕa, itp_ϕb, itp_qᵍ, h.issuance, h.w′, h.def_thres)

	function new_grid(x′, xgrid)
		xmax = maximum(x′)
		xmin = minimum(x′)

		Nx = length(xgrid)

		# Expand grids if x′ goes beyond the bounds
		xmax > maximum(xgrid)? Xmax = maximum(xgrid) + 0.25: Xmax = maximum(xgrid)# * (maximum(x′) - xmax): Void
		xmin < minimum(xgrid)? Xmin = minimum(xgrid) - 0.25: Xmin = minimum(xgrid)# * (xmin - minimum(x′)): Void

		# Retract grids if x′ doesn't reach the bounds
		xmax < maximum(xgrid)? Xmax = maximum(xgrid) - 0.05: Xmax = maximum(xgrid)# * (xmax - maximum(x′)): Void
		xmin > minimum(xgrid)? Xmin = minimum(xgrid) + 0.05: Xmin = minimum(xgrid)# * (minimum(x′) - xmin): Void

		return collect(linspace(Xmin, Xmax, Nx))
	end


	new_μgrid = new_grid(μ′_new, h.μgrid)
	new_σgrid = h.σgrid

	μ′_new = max.(min.(μ′_new, maximum(h.μgrid)), minimum(h.μgrid))
	σ′_new = max.(min.(σ′_new, maximum(h.σgrid)), minimum(h.σgrid))

	dist_exp[1] = sqrt.(sum( (μ′_new - μ′_old).^2 )) / sqrt.(sum(μ′_old.^2))
	dist_exp[2] = sqrt.(sum( (σ′_new - σ′_old).^2 )) / sqrt.(sum(σ′_old.^2))

	μ′_new = upd_η * μ′_new + (1.0 - upd_η) * μ′_old
	σ′_new = upd_η * σ′_new + (1.0 - upd_η) * σ′_old

	h.μ′ = μ′_new
	h.σ′ = σ′_new

	return dist_exp, new_μgrid, new_σgrid
end

function update_grids!(h::Hank; new_μgrid::Vector=[], new_σgrid::Vector=[], new_wgrid::Vector=[])

	if new_μgrid==[]
		new_μgrid = h.μgrid
	end
	if new_σgrid==[]
		new_σgrid = h.σgrid
	end
	if new_wgrid==[]
		new_wgrid = h.wgrid
	end

	function reinterp(h::Hank, y; agg::Bool=false)
		knots = (h.ωgrid, h.ϵgrid, h.bgrid, h.μgrid, h.σgrid, h.wgrid, h.ζgrid, h.zgrid)
		if agg
			knots = (h.bgrid, h.μgrid, h.σgrid, h.wgrid, h.ζgrid, h.zgrid)
			y = reshape(y, h.Nb, h.Nμ, h.Nσ, h.Nw, h.Nζ, h.Nz)
		end

		itp_obj_y = interpolate(knots, y, Gridded(Linear()))
		itp_y = extrapolate(itp_obj_y, Linear())

		if agg
			y_new = itp_y[h.bgrid, new_μgrid, new_σgrid, new_wgrid, h.ζgrid, h.zgrid]
			return reshape(y_new, size(h.Jgrid,1))
		else
			y_new = itp_y[h.ωgrid, h.ϵgrid, h.bgrid, new_μgrid, new_σgrid, new_wgrid, h.ζgrid, h.zgrid]
			return y_new
		end
	end

	h.ϕa = reinterp(h, h.ϕa, agg=false)
	h.ϕb = reinterp(h, h.ϕb, agg=false)
	h.ϕc = reinterp(h, h.ϕc, agg=false)
	h.vf = reinterp(h, h.vf, agg=false)

	h.Ld 		= reinterp(h, h.Ld, agg=true)
	h.wage 		= reinterp(h, h.wage, agg=true)
	h.issuance 	= reinterp(h, h.issuance, agg=true)
	h.issuance 	= min.(max.(h.issuance, minimum(h.bgrid)), maximum(h.bgrid))
	h.spending 	= reinterp(h, h.spending, agg=true)
	h.pN 		= reinterp(h, h.pN, agg=true)
	h.w′ 		= reinterp(h, h.w′, agg=true)

	knots 		= (h.bgrid, h.μgrid, h.σgrid, h.wgrid, h.ζgrid, h.zgrid, h.zgrid)
	repay_mat 	= reshape(h.repay, h.Nb, h.Nμ, h.Nσ, h.Nw, h.Nζ, h.Nz, h.Nz)
	itp_repay 	= interpolate(knots, repay_mat, Gridded(Linear()))
	rep_new 	= itp_repay[h.bgrid, h.μgrid, h.σgrid, h.wgrid, h.ζgrid, h.zgrid, h.zgrid]
	h.repay 	= reshape(rep_new, h.Nb*h.Nμ*h.Nσ*h.Nw*h.Nζ*h.Nz*h.Nz)

	for jzp in 1:h.Nz
		for jreent in 1:2
			h.μ′[:,jzp,jreent] = reinterp(h, h.μ′[:,jzp,jreent], agg=true)
			h.σ′[:,jzp,jreent] = reinterp(h, h.σ′[:,jzp,jreent], agg=true)
		end
	end

	h.μgrid = new_μgrid
	h.σgrid = new_σgrid
	h.wgrid = new_wgrid

	h.μ′ = max.(min.(h.μ′, maximum(h.μgrid)), minimum(h.μgrid))
	h.σ′ = max.(min.(h.σ′, maximum(h.σgrid)), minimum(h.σgrid))
	h.w′ = max.(min.(h.w′, maximum(h.wgrid)), minimum(h.wgrid))

	Void
end
