TFP_N(z, Δ, ζ) = 1.0    * max(0, (1.0 - Δ*1*(ζ==2) - Δ*(1.0-1)*exp(.25*z)^3 ))
TFP_T(z, Δ, ζ) = exp(z) * max(0, (1.0 - Δ*1*(ζ==2) - Δ*(1.0-1)*exp(z)^3))

function extend_state_space!(h::Hank, qʰ_mat, qᵍ_mat, T_mat)

	Npn = length(h.pngrid)

	ϕa_ext = Array{Float64}(undef, h.Nω, h.Nϵ, h.Nb, h.Nμ, h.Nσ, h.Nξ, h.Nζ, h.Nz, Npn)
	ϕb_ext = Array{Float64}(undef, h.Nω, h.Nϵ, h.Nb, h.Nμ, h.Nσ, h.Nξ, h.Nζ, h.Nz, Npn)
	ϕc_ext = Array{Float64}(undef, h.Nω, h.Nϵ, h.Nb, h.Nμ, h.Nσ, h.Nξ, h.Nζ, h.Nz, Npn)

	itp_vf = make_itp(h, h.vf; agg=false)
	itp_qᵍ = make_itp(h, h.qᵍ; agg=true)

	print_save("\nExtending the state space ($(Npn) iterations needed)")

	for jpn in 1:Npn

		pnv = h.pngrid[jpn]

		N = size(h.Jgrid, 1)

		wage_pn, labor_pn, profits_pn = Array{Float64, 1}(undef, N), Array{Float64, 1}(undef, N), Array{Float64, 1}(undef, N)
		for js in 1:N
			jζ = h.Jgrid[js, 5]
			jz = h.Jgrid[js, 6]

			wv = h.wbar
			ζv = h.ζgrid[jζ]
			zv = h.zgrid[jz]

			# wv = 0.9

			labor_pn[js], wage_pn[js], profits_pn[js], _ = labor_market(h, ζv, zv, pnv)
		end

		pC = price_index(h, pnv)
		pC_mat = ones(h.Nb, h.Nμ, h.Nσ, h.Nξ, h.Nζ, h.Nz) * pC

		T_mat = govt_bc(h, wage_pn.*labor_pn)# - reshape(profits_pn - h.profits, h.Nb, h.Nμ, h.Nσ, h.Nξ, h.Nζ, h.Nz)
		Π_mat = reshape(profits_pn, h.Nb, h.Nμ, h.Nσ, h.Nξ, h.Nζ, h.Nz)

		wL_mat  = reshape(wage_pn.*labor_pn, h.Nb, h.Nμ, h.Nσ, h.Nξ, h.Nζ, h.Nz) * (1.0 - h.τ)

		if h.ϕa_ext[1,1,1,1,1,1,1,1,jpn] == 0
			guess_a = h.ϕa
			guess_b = h.ϕb
		else
			guess_a = h.ϕa_ext[:,:,:,:,:,:,:,:,jpn]
			guess_b = h.ϕb_ext[:,:,:,:,:,:,:,:,jpn]
		end

		# Re-solve for these values of wn and pn
		_, ϕa, ϕb, ϕe, ϕc, _ = opt_value(h, qʰ_mat, qᵍ_mat, wL_mat, T_mat, pC_mat, Π_mat, itp_qᵍ, itp_vf; resolve = true, verbose = false, guess_a=guess_a, guess_b=guess_b)

		isapprox(sum(abs.(ϕc)), 0) ? print_save("\nWARNING: ϕc(pN = $(round(pnv, 2))) ≡ 0 when extending state space") : nothing

		for jz in 1:h.Nz, jζ in 1:h.Nζ, jξ in 1:h.Nξ, jσ in 1:h.Nσ, jμ in 1:h.Nμ, jb in 1:h.Nb, jϵ in 1:h.Nϵ, jω in 1:h.Nω
			ϕa_ext[jω,jϵ,jb,jμ,jσ,jξ,jζ,jz,jpn] = ϕa[jω,jϵ,jb,jμ,jσ,jξ,jζ,jz]
			ϕb_ext[jω,jϵ,jb,jμ,jσ,jξ,jζ,jz,jpn] = ϕb[jω,jϵ,jb,jμ,jσ,jξ,jζ,jz]
			ϕc_ext[jω,jϵ,jb,jμ,jσ,jξ,jζ,jz,jpn] = ϕc[jω,jϵ,jb,jμ,jσ,jξ,jζ,jz]
		end
	end

	!isnan(sum(ϕa_ext)) || print_save("ERROR: $(isnan(sum(ϕa_ext))) NaN counts in ϕa_ext")
	!isnan(sum(ϕa_ext)) || throw(error("$(isnan(sum(ϕa_ext))) NaN counts in ϕa_ext"))

	isapprox(sum(abs.(ϕc_ext)), 0) ? print_save("\nWARNING: ϕc ≡ 0 when extending state space") : nothing

	h.ϕa_ext = ϕa_ext
	h.ϕb_ext = ϕb_ext
	h.ϕc_ext = ϕc_ext

	nothing
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

function labor_market(h::Hank, ζv, zv, pNv; w_slack::Bool=false)
	""" Finds w and Lᵈ at the current state given a guess of pNv """
	w_constraint = h.wbar
	if w_slack
		w_constraint = 0.5
	end

	# Step 1: Assume wₜ is at the constraint, find labor demand, and check whether the eq'm wage is above or below
	Ld_N, Ld_T = labor_demand(h, w_constraint, zv, ζv, pNv; get_both=true)
	Ld = Ld_N + Ld_T

	# Step 2: If labor demand is lower than supply, find the wage above γw * wv that clears the labor mkt
	Ls = 1.0
	αmax = max(h.α_N, h.α_T)
	w_max = 1.1 * αmax * ( TFP_T(zv,h.Δ,ζv)^(1.0/(1.0-αmax)) + (pNv*TFP_N(zv,h.Δ,ζv))^(1.0/(1.0-αmax)) )^(1.0-αmax) 

	if w_max < w_constraint
		# print_save("\nSomething wrong with wages")
		w_max = w_constraint
	end

	w_new = w_constraint

	if Ld - Ls > 1e-4
		res = Optim.optimize(
			w -> (labor_demand(h, w, zv, ζv, pNv) - Ls)^2,
				w_constraint, 1.0.*w_max, GoldenSection()
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


function mkt_clearing(h::Hank, itp_ϕc, G, Bpv, pNv, pNmin, pNmax, bv, μv, σv, ξv, jζ, jz, jdefault; w_slack::Bool=false, orig_vars::Bool = true, get_others::Bool = false, get_both::Bool=false)
	typeof(pNv) == Vector{Float64} ? pN = pNv[1] : pN = pNv
	if orig_vars == false
		pN = transform_vars(pN, pNmin, pNmax)
	end

	isnan(pN) ? print_save("\nWARNING: pNv, pN = $(pNv[1]), $(pN)") : nothing

	ζv, zv = h.ζgrid[jζ], h.zgrid[jz]

	Ld, w_new, profits, output = labor_market(h, ζv, zv, pN; w_slack=w_slack)

	# Check market clearing for nontradables
	Ld_N, _  = labor_demand(h, w_new, zv, ζv, pN; get_both=true)
	supply_N = TFP_N(zv, h.Δ, ζv) * Ld_N^(h.α_N)

	# Get the household's policies at these prices
	ωmin_int, ωmax_int = quantile.(LogNormal(μv, σv), [1e-6; 1-1e-6]) .+ h.ωmin
	ωmax_int = min(ωmax_int, maximum(h.ωgrid))
	ωmin_int = h.ωmin
	val_C, sum_prob = 0., 0.
	itp_ϕc = extrapolate(itp_ϕc, Interpolations.Flat())
	for (jϵ, ϵv) in enumerate(h.ϵgrid)
		f_pdf(ω) = pdf(LogNormal(μv, σv), ω-h.ωmin)
		(val_pdf, err) = hquadrature(f_pdf, ωmin_int, ωmax_int, rtol=1e-10, atol=1e-12, maxevals=500)
		sum_prob += val_pdf * h.λϵ[jϵ]

		f(ω) = f_pdf(ω) * itp_ϕc(ω, jϵ, bv, μv, σv, ξv, jζ, jz, pN)
		(val, err) = hquadrature(f, ωmin_int, ωmax_int, rtol=1e-12, atol=0, maxevals=500)
	
		val_C += val * h.λϵ[jϵ]
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
	# isapprox(sum_prob, 1) || print_save("\nWARNING: Something wrong computing aggregate consumption")

	val_int_C = val_C / sum_prob
 	
 	# Recover nontraded demand from total consumption
	pC = price_index(h, pN)
	demand_N_cons = val_int_C * h.ϖ * (pN/pC)^(-h.η)
	demand_N_govt = G * h.ϑ / pN

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

function find_prices(h::Hank, itp_ϕc, G, Bpv, pNg, pNmin, pNmax, bv, μv, σv, ξv, jζ, jz, jdefault)

	# First find eq'm assuming the constraint does not bind
	res = Optim.optimize(
		pN -> mkt_clearing(h, itp_ϕc, G, Bpv, pN, pNmin, pNmax, bv, μv, σv, ξv, jζ, jz, jdefault; orig_vars = true, w_slack=true)^2,
		pNmin, pNmax, GoldenSection()
		)
	pN = res.minimizer
	wage, Ld, output = mkt_clearing(h, itp_ϕc, G, Bpv, pN, pNmin, pNmax, bv, μv, σv, ξv, jζ, jz, jdefault; orig_vars=true, get_others=true, w_slack=true)

	if wage >= h.wbar && res.minimum < 1e-6
		pN >= pNmax - 0.05*(pNmax-pNmin) ? exc_dem = 1 : exc_dem = 0
		pN <= pNmin + 0.05*(pNmax-pNmin) ? exc_sup = 1 : exc_sup = 0

		minf = mkt_clearing(h, itp_ϕc, G, Bpv, pN, pNmin, pNmax, bv, μv, σv, ξv, jζ, jz, jdefault; orig_vars = true)

		results = [wage; pN; Ld; output]

		return results, [minf], exc_dem, exc_sup
	end

	res = Optim.optimize(
		pN -> mkt_clearing(h, itp_ϕc, G, Bpv, pN, pNmin, pNmax, bv, μv, σv, ξv, jζ, jz, jdefault; orig_vars = true)^2,
		# 0.9*pNmin, 1.1*pNmax, GoldenSection()
		pNmin, pNmax, GoldenSection()
		)
	pN = res.minimizer
	minf = mkt_clearing(h, itp_ϕc, G, Bpv, pN, pNmin, pNmax, bv, μv, σv, ξv, jζ, jz, jdefault; orig_vars = true)

	pN >= pNmax - 0.05*(pNmax-pNmin) ? exc_dem = 1 : exc_dem = 0
	pN <= pNmin + 0.05*(pNmax-pNmin) ? exc_sup = 1 : exc_sup = 0
	# if res.minimum > 1e-6
	# 	exc_dem, exc_sup = 1, 1
	# end
	wage, Ld, output = mkt_clearing(h, itp_ϕc, G, Bpv, pN, pNmin, pNmax, bv, μv, σv, ξv, jζ, jz, jdefault; get_others=true)

	results = [wage; pN; Ld; output]

	if abs(minf) > 1e-4
		# print_save("\nNontradables exc supply = $(@sprintf("%0.4g",minf)) at pN = $(@sprintf("%0.4g",pN))")
	end

	return results, [minf], exc_dem, exc_sup
end

function eval_prices_direct(h::Hank, itp_ϕc, G, pN, bv, μv, σv, ξv, jζ, jz, jdefault; get_others::Bool=false)

	ζv, zv = h.ζgrid[jζ], h.zgrid[jz]

	Ls = 1.0
	
	# Step 1: Find unconstrained wage
	res = Optim.optimize(
		w -> (Ls - labor_demand(h, w, zv, ζv, pN))^2,
		0.5, 4*pN, GoldenSection()
		)

	# Step 2: Apply constraint
	wv = max(res.minimizer, h.wbar)

	# Step 3: Compute labor demand in nontradables and supply
	Ld_N, Ld_T  = labor_demand(h, wv, zv, ζv, pN; get_both=true)
	Ld = Ld_N + Ld_T
	supply_N = TFP_N(zv,h.Δ,ζv) * Ld_N^h.α_N
	supply_T = TFP_T(zv,h.Δ,ζv) * Ld_T^h.α_T

	# Step 4: Get traded absorption
	# Get the household's policies
	ωmin_int, ωmax_int = quantile.(LogNormal(μv, σv), [1e-6; 1-1e-6]) .+ h.ωmin
	ωmax_int = min(ωmax_int, maximum(h.ωgrid))
	ωmin_int = h.ωmin
	val_C, sum_prob = 0., 0.
	itp_ϕc = extrapolate(itp_ϕc, Interpolations.Line())
	for (jϵ, ϵv) in enumerate(h.ϵgrid)
		f_pdf(ω) = pdf(LogNormal(μv, σv), ω-h.ωmin)
		(val_pdf, err) = hquadrature(f_pdf, ωmin_int, ωmax_int, rtol=1e-10, atol=1e-12, maxevals=500)
		sum_prob += val_pdf * h.λϵ[jϵ]

		f(ω) = f_pdf(ω) * itp_ϕc(ω, jϵ, bv, μv, σv, ξv, jζ, jz)
		(val, err) = hquadrature(f, ωmin_int, ωmax_int, rtol=1e-12, atol=0, maxevals=500)
	
		val_C += val * h.λϵ[jϵ]
	end

	val_int_C = val_C / sum_prob

 	# Recover traded demand from total consumption
	pC = price_index(h, pN)
	cT = val_int_C * (1.0-h.ϖ) * (1.0/pC)^(-h.η)
	cN = val_int_C * (h.ϖ) * (pN/pC)^(-h.η)
	gN = G * h.ϑ / pN

	yN = max(0.0, supply_N - gN)

	# Get new implied pN
	pN_new = h.ϖ^(1/h.η) / (1-h.ϖ^(1/h.η)) * (cT/yN)^(1/h.η)

	pN_new = min(100.0, pN_new)

	F = pN - pN_new

	if get_others
		output = pN * supply_N + supply_T
		return wv, Ld, output, supply_N, cN + gN, pN_new
	else
		return F
	end
end

function find_prices_direct(h::Hank, itp_ϕc, G, Bpv, pNg, pNmin, pNmax, bv, μv, σv, ξv, jζ, jz, jdefault)

	do_solver_prices = true

	if do_solver_prices
		res = Optim.optimize(
			pNv -> (eval_prices_direct(h, itp_ϕc, G, pNv, bv, μv, σv, ξv, jζ, jz, jdefault))^2,
			pNmin, pNmax, GoldenSection()
			)

		pN = res.minimizer
	else
		if false
			pNg = max(min(pNg, pNmax - 1e-6), pNmin + 1e-6)
			obj_f(x) = (eval_prices_direct(h, itp_ϕc, G, x, bv, μv, σv, ξv, jζ, jz, jdefault))^2
			res = Optim.optimize(
				pN -> obj_f(first(pN)),
				[pNmin], [pNmax], [pNg], Fminbox(LBFGS())#, autodiff=:forward#, Optim.Options(f_tol=1e-12)
				)
			pN = first(res.minimizer)
		else
			pN = pNg
		end
	end

	wage, Ld, output, supply_N, demand_N, _ = eval_prices_direct(h, itp_ϕc, G, pN, bv, μv, σv, ξv, jζ, jz, jdefault; get_others=true)

	pN >= pNmax - 0.05*(pNmax-pNmin) ? exc_dem = 1 : exc_dem = 0
	pN <= pNmin + 0.05*(pNmax-pNmin) ? exc_sup = 1 : exc_sup = 0

	minf = supply_N - demand_N

	results = [wage; pN; Ld; output]

	return results, [minf], exc_dem, exc_sup
end

function find_all_prices(h::Hank, itp_ϕc, B′_vec, G_vec)

	N = size(h.Jgrid, 1)

	results = Array{Float64}(undef, N, 4)
	minf	= Array{Float64}(undef, N, 1)
	exc_dem = Array{Float64}(undef, N)
	exc_sup = Array{Float64}(undef, N)

	pN_guess = h.pN
	
	# for js in 1:N
	Threads.@threads for js in 1:N
		Bpv = B′_vec[js]
		G = G_vec[js]
		
		pNmin, pNmax = minimum(h.pngrid), maximum(h.pngrid)
		isnan(pN_guess[js]) ? pNg = (pNmin+pNmax) / 2 : pNg = 1 * pN_guess[js]

		jb = h.Jgrid[js, 1]
		jμ = h.Jgrid[js, 2]
		jσ = h.Jgrid[js, 3]
		jξ = h.Jgrid[js, 4]
		jζ = h.Jgrid[js, 5]
		jz = h.Jgrid[js, 6]

		bv = h.bgrid[jb]
		μv = h.μgrid[jμ]
		σv = h.σgrid[jσ]
		ξv = h.ξgrid[jξ]
		ζv = h.ζgrid[jζ]
		zv = h.zgrid[jz]

		jdefault = (ζv != 1.0)

		# r, mf, ed, es = find_prices(h, itp_ϕc, G, Bpv, pNg, pNmin, pNmax, bv, μv, σv, ξv, jζ, jz, jdefault)
		r, mf, ed, es = find_prices_direct(h, itp_ϕc, G, Bpv, pNg, pNmin, pNmax, bv, μv, σv, ξv, jζ, jz, jdefault)

		results[js, :] = r
		minf[js, :], exc_dem[js], exc_sup[js] = mf, ed, es
		if isnan(results[js, 2])
			results[js, 2] = pNg
		end
	end

	return results, minf, exc_dem, exc_sup
end

function update_state_functions!(h::Hank, upd_η::Float64)
	# itp_ϕc = make_itp(h, h.ϕc_ext; agg=false)
	itp_ϕc = make_itp(h, h.ϕc; agg=false)
	# itp_ϕc = extrapolate(itp_ϕc, Interpolations.Flat())

	t1 = time()
	results, minf, exc_dem, exc_sup = find_all_prices(h, itp_ϕc, h.issuance, h.spending)
	print_save(" (new prices in $(time_print(time()-t1)))")

	dist = Array{Float64,1}(undef, 3)
	dist[1] = sqrt.(sum( (results[:, 1] - h.wage).^2 )) / sqrt.(sum(h.wage.^2))
	dist[2] = sqrt.(sum( (results[:, 2] - h.pN).^2 ))   / sqrt.(sum(h.pN.^2))
	dist[3] = sqrt.(sum( (results[:, 3] - h.Ld).^2 ))   / sqrt.(sum(h.Ld.^2))

	h.pN = upd_η * results[:, 2] + (1.0-upd_η) * h.pN

	consistent_others = true
	if consistent_others
		N = size(h.Jgrid,1)
		wage, Ld, output = zeros(N), zeros(N), zeros(N)
		Threads.@threads for js in 1:N
		# for js in 1:N
			Bpv = h.issuance[js]
			G = h.spending[js]

			pN = h.pN[js]
			pNmin, pNmax = minimum(h.pngrid), maximum(h.pngrid)

			bv = h.bgrid[h.Jgrid[js, 1]]
			μv = h.μgrid[h.Jgrid[js, 2]]
			σv = h.σgrid[h.Jgrid[js, 3]]
			ξv = h.ξgrid[h.Jgrid[js, 4]]
			jζ = h.Jgrid[js, 5]
			jz = h.Jgrid[js, 6]

			pN = max(min(pN, pNmax), pNmin)

			# wage[js], Ld[js], output[js] = mkt_clearing(h, itp_ϕc, G, Bpv, pN, pNmin, pNmax, bv, μv, σv, ξv, jζ, jz, (jζ!=1); get_others=true)

			wage[js], Ld[js], output[js], _, _, _ = eval_prices_direct(h, itp_ϕc, G, pN, bv, μv, σv, ξv, jζ, jz, (jζ!=1); get_others=true)

			jzmean = h.Nz # ceil(Int, h.Nz/2)
		end
		h.wage, h.Ld, h.output = wage, Ld, output
		update_fiscalrules!(h)
	else
		h.wage 	 = upd_η * results[:, 1] + (1.0-upd_η) * h.wage
		h.Ld 	 = upd_η * results[:, 3] + (1.0-upd_η) * h.Ld
		h.output = upd_η * results[:, 4] + (1.0-upd_η) * h.output
		h.w′ 	 = h.w′
	end

	h.profits = h.output - h.wage .* h.Ld
	mean_f = mean(minf)
	max_f = maximum(abs.(minf))

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
	ζ = 1
	wlow, whigh = minimum(h.wage), maximum(h.wage)

	res1 = Optim.optimize(
			w -> (labor_demand(h, w, h.zgrid[end], ζ, pN_up) - Ls).^2,
			wlow, whigh * 2.0, GoldenSection()
			)
	res2 = Optim.optimize(
			w -> (labor_demand(h, w, h.zgrid[end], ζ, pN_down) - Ls).^2,
			wlow, whigh * 2.0, GoldenSection()
			)
	w_up = max(res1.minimizer, res2.minimizer) * 1.25
	ζ = 2
	res1 = Optim.optimize(
			w -> (labor_demand(h, w, h.zgrid[1], ζ, pN_up) - Ls).^2,
			0.5 * wlow, whigh, GoldenSection()
			)
	res2 = Optim.optimize(
			w -> (labor_demand(h, w, h.zgrid[1], ζ, pN_down) - Ls).^2,
			0.5 * wlow, whigh, GoldenSection()
			)
	w_down = min(res1.minimizer, res2.minimizer) * 0.75

	h.pngrid = collect(range(pN_down, pN_up, length=length(h.pngrid)))

	nothing
end


function find_q(h::Hank, q, a, b, var_a, var_b, cov_ab, Bpv, ξpv, ζpv, jzp, jdef, itp_qᵍ, reentry; get_μσ::Bool=false)

	zpv = h.zgrid[jzp]

	haircut = 0.0
	if jdef && reentry==false
		ζpv = 2
		haircut = 0.0
	end
	if jdef == false && ζpv == 2
		ζpv = 2
		haircut = h.ℏ
	end

	R = (ζpv==1) * h.κ + (1.0 - haircut) .* ((1.0-h.ρ)*q)

	Eω   = a + R*b
	varω = var_a + R^2 * var_b + 2*R * cov_ab
	varω = max(varω, 0.)

	# varω >= 0. || print_save("\nvar_a, var_b, cov_ab, R, q = $(var_a), $(var_b), $(cov_ab), $(R), $(q)")

	μpv, σpv = make_logN(max(0.0, Eω - h.ωmin), varω)
	# μpv = min(max(μpv, minimum(h.μgrid)), maximum(h.μgrid))
	# σpv = min(max(σpv, minimum(h.σgrid)), maximum(h.σgrid))

	itp_qᵍ = extrapolate(itp_qᵍ, Interpolations.Flat())
	new_q = itp_qᵍ((1.0 - haircut) .* Bpv, μpv, σpv, ξpv, ζpv, jzp)

	if get_μσ
		return μpv, σpv
	else
		return new_q
	end
end


function compute_stats_logN(h::Hank, js, a, b, var_a, var_b, cov_ab, itp_qᵍ, Bpv, exp_rep)

	ζv = h.ζgrid[h.Jgrid[js, 5]]
	jdef = (ζv != 1.0)


	μ, σ, q = zeros(h.Nξ, h.Nz, 2), zeros(h.Nξ, h.Nz, 2), zeros(h.Nξ, h.Nz, 2)
	alarm_mat = Array{Float64, 3}(undef, h.Nξ, h.Nz, 2)

	for (jξp, ξpv) in enumerate(h.ξgrid), (jzp, zpv) in enumerate(h.zgrid)
		# First any case where ζ′ = 1
		reentry = true
		ζpv = 1
		qmin, qmax = minimum(h.qᵍ), maximum(h.qᵍ)

		res = Optim.optimize(
			q -> (find_q(h, q, a, b, var_a, var_b, cov_ab, Bpv, ξpv, ζpv, jzp, jdef, itp_qᵍ, reentry) - q)^2,
			qmin, qmax, Brent()
			)
		q[jξp, jzp, 1] = res.minimizer
		res.minimum > 1e-4 ? alarm_mat[jξp, jzp, 1] = 1 : alarm_mat[jξp, jzp, 1] = 0

		μ[jξp, jzp, 1], σ[jξp, jzp, 1] = find_q(h, q[jξp, jzp, 1], a, b, var_a, var_b, cov_ab, Bpv, ξpv, ζpv, jzp, jdef, itp_qᵍ, reentry; get_μσ = true)

		if jdef
			# If default continues
			reentry = false
			ζpv = 2 # Irrelevant but to stress that the default state continues
			res = Optim.optimize(
				q -> (find_q(h, q, a, b, var_a, var_b, cov_ab, Bpv, ξpv, ζpv, jzp, jdef, itp_qᵍ, reentry) - q)^2,
				qmin, qmax, Brent()
				)
			q[jξp, jzp, 2] = res.minimizer
			res.minimum > 1e-4 ? alarm_mat[jξp, jzp, 2] = 1 : alarm_mat[jξp, jzp, 2] = 0

			μ[jξp, jzp, 2], σ[jξp, jzp, 2] = find_q(h, q[jξp, jzp, 2], a, b, var_a, var_b, cov_ab, Bpv, ξpv, ζpv, jzp, jdef, itp_qᵍ, reentry; get_μσ = true)
		else
			# Entering default
			reentry = true # Irrelevant
			ζpv = 2
			res = Optim.optimize(
				q -> (find_q(h, q, a, b, var_a, var_b, cov_ab, Bpv, ξpv, ζpv, jzp, jdef, itp_qᵍ, reentry) - q)^2,
				qmin, qmax, Brent()
				)
			q[jξp, jzp, 2] = res.minimizer
			res.minimum > 1e-4 ? alarm_mat[jξp, jzp, 2] = 1 : alarm_mat[jξp, jzp, 2] = 0

			μ[jξp, jzp, 2], σ[jξp, jzp, 2] = find_q(h, q[jξp, jzp, 2], a, b, var_a, var_b, cov_ab, Bpv, ξpv, ζpv, jzp, jdef, itp_qᵍ, reentry; get_μσ = true)
		end
	end
	return μ, σ, q, alarm_mat
end

function new_expectations(h::Hank, itp_ϕa, itp_ϕb, itp_qᵍ, Bpv, exp_rep, js, jdef)

	itp_ϕa = extrapolate(itp_ϕa, Interpolations.Line())
	itp_ϕb = extrapolate(itp_ϕb, Interpolations.Line())

	jb = h.Jgrid[js, 1]
	jμ = h.Jgrid[js, 2]
	jσ = h.Jgrid[js, 3]
	jξ = h.Jgrid[js, 4]
	jζ = h.Jgrid[js, 5]
	jz = h.Jgrid[js, 6]

	bv = h.bgrid[jb]
	μv = h.μgrid[jμ]
	σv = h.σgrid[jσ]
	ξv = h.ξgrid[jξ]

	val_a, val_b, val_a2, val_b2, val_ab, sum_prob = 0., 0., 0., 0., 0., 0.

	ωmin_int, ωmax_int = quantile.(LogNormal(μv, σv), [1e-6; 1-1e-6]) .+ h.ωmin
	ωmax_int = min(ωmax_int, maximum(h.ωgrid))
	ωmin_int = h.ωmin
	for (jϵ, ϵv) in enumerate(h.ϵgrid)
		f_pdf(ω) = pdf(LogNormal(μv, σv), ω-h.ωmin)
		(val_pdf, err) = hquadrature(f_pdf, ωmin_int, ωmax_int, rtol=1e-10, atol=1e-12, maxevals=0)
		sum_prob += val_pdf * h.λϵ[jϵ]

		fA(ω) = f_pdf(ω) * max(h.ωmin, itp_ϕa(ω, jϵ, bv, μv, σv, ξv, jζ, jz))
		(valA, err) = hquadrature(fA, ωmin_int, ωmax_int, rtol=1e-10, atol=1e-12, maxevals=0)
		val_a += valA * h.λϵ[jϵ]

		fA2(ω) = f_pdf(ω) * max(h.ωmin, itp_ϕa(ω, jϵ, bv, μv, σv, ξv, jζ, jz))^2
		(valA2, err) = hquadrature(fA2, ωmin_int, ωmax_int, rtol=1e-10, atol=1e-12, maxevals=0)
		val_a2 += valA2 * h.λϵ[jϵ]

		fB(ω) = f_pdf(ω) * max(0., itp_ϕb(ω, jϵ, bv, μv, σv, ξv, jζ, jz))
		(valB, err) = hquadrature(fB, ωmin_int, ωmax_int, rtol=1e-10, atol=1e-12, maxevals=0)
		val_b += valB * h.λϵ[jϵ]

		fB2(ω) = f_pdf(ω) * max(0., itp_ϕb(ω, jϵ, bv, μv, σv, ξv, jζ, jz))^2
		(valB2, err) = hquadrature(fB2, ωmin_int, ωmax_int, rtol=1e-10, atol=1e-12, maxevals=0)
		val_b2 += valB2 * h.λϵ[jϵ]

		fAB(ω) = f_pdf(ω) * max(h.ωmin, itp_ϕa(ω, jϵ, bv, μv, σv, ξv, jζ, jz)) * max(0., itp_ϕb(ω, jϵ, bv, μv, σv, ξv, jζ, jz))
		(valAB, err) = hquadrature(fAB, ωmin_int, ωmax_int, rtol=1e-10, atol=1e-12, maxevals=0)
		val_ab += valAB * h.λϵ[jϵ]
	end

	# for (jϵ, ϵv) in enumerate(h.ϵgrid)
	# 	for jω = 1:length(h.ωgrid_fine)-1
	# 		ωv  = h.ωgrid_fine[jω]
	# 		ω1v = h.ωgrid_fine[jω+1]
	# 		ωmv = 0.5*(ωv+ω1v)

	# 		prob = pdf(LogNormal(μv, σv), ωmv-h.ωmin) * h.λϵ[jϵ] * (ω1v - ωv)

	# 		ϕa = max(itp_ϕa[ωmv, jϵ, bv, μv, σv, ξv, jζ, jz], h.ωmin)
	# 		ϕb = max(itp_ϕb[ωmv, jϵ, bv, μv, σv, ξv, jζ, jz], 0.)

	# 		val_a  += prob * ϕa
	# 		val_a2 += prob * ϕa^2
	# 		val_b  += prob * ϕb
	# 		val_b2 += prob * ϕb^2
	# 		val_ab += prob * ϕa * ϕb

	# 		sum_prob += prob
	# 	end
	# end

	!isnan(sum_prob) || throw(error("NaN in sum_prob"))
	!isapprox(sum_prob, 0.) || throw(error("\nsum_prob = $(sum_prob) at $([jb, jμ, jσ, jξ, jζ, jz])"))

	!isnan(val_a+val_a2+val_b+val_b2+val_ab) || print_save("\na,a2,b,b2,ab = $([val_a,val_a2,val_b,val_b2,val_ab])")

	a  = val_a  / sum_prob
	a2 = val_a2 / sum_prob
	b  = val_b  / sum_prob
	b2 = val_b2 / sum_prob
	ab = val_ab / sum_prob

	var_a  = a2 - a^2
	var_b  = b2 - b^2
	cov_ab = ab - a*b

	if max(var_a, var_b) < 1e-8
		var_a, var_b, cov_ab = 0., 0., 0.
	end

	!isnan(var_a+var_b+cov_ab) || print_save("\nVa, Vb, cov = $var_a, $var_b, $cov_ab at $([jb, jμ, jσ, jξ, jζ, jz])")

	μ′, σ′, q, alarm_vec = compute_stats_logN(h, js, a, b, var_a, var_b, cov_ab, itp_qᵍ, Bpv, exp_rep)

	return μ′, σ′, alarm_vec
end


function find_all_expectations(h::Hank, itp_ϕa, itp_ϕb, itp_qᵍ, B′_vec)
	N = size(h.Jgrid, 1)

	μ′ = Array{Float64}(undef, N, h.Nξ, h.Nz, 2)
	σ′ = Array{Float64}(undef, N, h.Nξ, h.Nz, 2)
	alarm_vec = Array{Float64}(undef, N, h.Nξ, h.Nz, 2)

	# Threads.@threads for js in 1:N
	for js in 1:N
		Bpv = B′_vec[js]

		jb = h.Jgrid[js, 1]
		jμ = h.Jgrid[js, 2]
		jσ = h.Jgrid[js, 3]
		jξ = h.Jgrid[js, 4]
		jζ = h.Jgrid[js, 5]
		jz = h.Jgrid[js, 6]

		rep_mat = reshape(h.repay, h.Nb, h.Nμ, h.Nσ, h.Nξ, h.Nζ, h.Nz, h.Nξ, h.Nz)
		exp_rep = rep_mat[jb, jμ, jσ, jξ, jζ, jz, :, :]

		jdefault = (jζ != 1.0)

		μ′[js, :, :, :], σ′[js, :, :, :], alarm_vec[js, :, :, :] = new_expectations(h, itp_ϕa, itp_ϕb, itp_qᵍ, Bpv, exp_rep, js, jdefault)
	end

	print_save("\n")
	if sum(alarm_vec[:]) >= 1
		print_save("WARNING: ")
	end
	print_save("Couldn't find qᵍ $(round(100*sum(alarm_vec[:])/N,digits=0))% of the time")

	return μ′, σ′
end

function update_expectations!(h::Hank, upd_η::Float64)
	"""
	Computes mean and variance of tomorrow's distribution and deduces parameters for logN
	"""

	μ′_old = copy(h.μ′)
	σ′_old = copy(h.σ′)

	dist_exp = Array{Float64,1}(undef, 2)

	itp_ϕa = make_itp(h, h.ϕa; agg=false)
	itp_ϕb = make_itp(h, h.ϕb; agg=false)
	itp_qᵍ = make_itp(h, h.qᵍ; agg=true)

	μ′_new, σ′_new = find_all_expectations(h, itp_ϕa, itp_ϕb, itp_qᵍ, h.issuance)

	function new_grid(x′, xgrid; ub::Float64=Inf, lb::Float64=-Inf)
		Nx = length(xgrid)
		xdist = maximum(xgrid) - minimum(xgrid)

		Xmax = maximum(xgrid)
		Xmin = minimum(xgrid)

		# Expand grids if x′ goes beyond the bounds
		quantile(x′[:], 0.99) > maximum(xgrid) ? Xmax = maximum(xgrid) + 0.01*xdist : nothing
		quantile(x′[:], 0.01) < minimum(xgrid) ? Xmin = minimum(xgrid) - 0.01*xdist : nothing

		# Retract grids if x′ doesn't reach the bounds
		maximum(x′) < maximum(xgrid) ? Xmax = maximum(xgrid) - 0.01*xdist : nothing
		minimum(x′) > minimum(xgrid) ? Xmin = minimum(xgrid) + 0.01*xdist : nothing

		Xmax = min(Xmax, ub)
		Xmin = max(Xmin, lb)

		return collect(range(Xmin, Xmax, length=Nx))
	end


	new_μgrid = new_grid(μ′_new, h.μgrid, lb = -3.0, ub = 3.5)
	new_σgrid = new_grid(σ′_new, h.σgrid, lb = 1e-2)

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

function update_grids!(h::Hank; new_μgrid::Vector=[], new_σgrid::Vector=[], new_zgrid::Vector=[])

	if new_μgrid==[]
		new_μgrid = h.μgrid
	end
	if new_σgrid==[]
		new_σgrid = h.σgrid
	end
	if new_zgrid==[]
		new_zgrid = h.zgrid
	end

	function reinterp(h::Hank, y; agg::Bool=false, ext::Bool=false)
		knots = (h.ωgrid, h.ϵgrid, h.bgrid, h.μgrid, h.σgrid, h.ξgrid, h.ζgrid, h.zgrid)
		if agg
			knots = (h.bgrid, h.μgrid, h.σgrid, h.ξgrid, h.ζgrid, h.zgrid)
			y = reshape(y, h.Nb, h.Nμ, h.Nσ, h.Nξ, h.Nζ, h.Nz)
		end
		if ext
			knots = (h.ωgrid, h.ϵgrid, h.bgrid, h.μgrid, h.σgrid, h.ξgrid, h.ζgrid, h.zgrid, h.pngrid)
		end

		itp_obj_y = interpolate(knots, y, Gridded(Linear()))
		itp_y = extrapolate(itp_obj_y, Line())

		if agg
			y_new = itp_y(h.bgrid, new_μgrid, new_σgrid, h.ξgrid, h.ζgrid, new_zgrid)
			return reshape(y_new, length(y_new))
		elseif ext
			y_new = itp_y(h.ωgrid, h.ϵgrid, h.bgrid, new_μgrid, new_σgrid, h.ξgrid, h.ζgrid, new_zgrid, h.pngrid)
			return y_new
		else
			y_new = itp_y(h.ωgrid, h.ϵgrid, h.bgrid, new_μgrid, new_σgrid, h.ξgrid, h.ζgrid, new_zgrid)
			return y_new
		end
	end

	h.ϕa = reinterp(h, h.ϕa, agg=false)
	h.ϕb = reinterp(h, h.ϕb, agg=false)
	h.ϕa_ext = reinterp(h, h.ϕa_ext, agg=false, ext=true)
	h.ϕb_ext = reinterp(h, h.ϕb_ext, agg=false, ext=true)
	h.ϕc = reinterp(h, h.ϕc, agg=false)
	h.ϕc = max.(1e-6, h.ϕc)
	h.vf = reinterp(h, h.vf, agg=false)
	h.vf = max.(1e-20, h.vf)

	h.Ld 		= reinterp(h, h.Ld, agg=true)
	h.output 	= reinterp(h, h.output, agg=true)
	h.wage 		= reinterp(h, h.wage, agg=true)
	h.issuance 	= reinterp(h, h.issuance, agg=true)
	h.issuance 	= min.(max.(h.issuance, minimum(h.bgrid)), maximum(h.bgrid))
	h.spending 	= reinterp(h, h.spending, agg=true)
	h.pN 		= reinterp(h, h.pN, agg=true)

	knots 		= (h.bgrid, h.μgrid, h.σgrid, h.ξgrid, h.ζgrid, h.zgrid, h.ξgrid, h.zgrid)
	repay_mat 	= reshape(h.repay, h.Nb, h.Nμ, h.Nσ, h.Nξ, h.Nζ, h.Nz, h.Nξ, h.Nz)
	itp_repay 	= extrapolate(interpolate(knots, repay_mat, Gridded(Linear())), Line())
	rep_new 	= itp_repay(h.bgrid, new_μgrid, new_σgrid, h.ξgrid, h.ζgrid, new_zgrid, h.ξgrid, new_zgrid)
	rep_new 	= max.(0, min.(1, rep_new))
	h.repay 	= reshape(rep_new, length(rep_new))

	h.welfare   = reinterp(h, h.welfare, agg=true)

	μ′_new = zeros(h.Nb * length(new_μgrid) * length(new_σgrid) * h.Nξ * h.Nζ * length(new_zgrid), h.Nξ, length(new_zgrid), 2)
	σ′_new = zeros(h.Nb * length(new_μgrid) * length(new_σgrid) * h.Nξ * h.Nζ * length(new_zgrid), h.Nξ, length(new_zgrid), 2)
	for jξp in 1:h.Nξ
		for jzp in 1:length(new_zgrid)
			for jreent in 1:2
				try
					μ′_new[:,jξp,jzp,jreent] = reinterp(h, h.μ′[:,jξp,jzp,jreent], agg=true)
					σ′_new[:,jξp,jzp,jreent] = reinterp(h, h.σ′[:,jξp,jzp,jreent], agg=true)
				catch
					μ′_new[:,jξp,jzp,jreent] = reinterp(h, h.μ′[:,jξp,h.Nz,jreent], agg=true)
					σ′_new[:,jξp,jzp,jreent] = reinterp(h, h.σ′[:,jξp,h.Nz,jreent], agg=true)
				end	
			end
		end
	end

	h.μgrid = new_μgrid
	h.σgrid = new_σgrid

	h.μ′ = max.(min.(μ′_new, maximum(h.μgrid)), minimum(h.μgrid))
	h.σ′ = max.(min.(σ′_new, maximum(h.σgrid)), minimum(h.σgrid))

	nothing
end
