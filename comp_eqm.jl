TFP_N(z, Δ, ζv) = 1.0    * max(0, (1.0 - Δ*(ζv==1)))
TFP_T(z, Δ, ζv) = exp(z) * max(0, (1.0 - Δ*(ζv==1)))

function labor_demand(sd::SOEdef, w, z, ζ, pN)
	pars = sd.pars

	Ld_nontradables = (pars[:α_N] * pN * TFP_N(z,pars[:Δ],ζ) / w).^(1.0/(1.0-pars[:α_N]))
	Ld_tradables    = (pars[:α_T] * 1  * TFP_T(z,pars[:Δ],ζ) / w).^(1.0/(1.0-pars[:α_T]))

	return Ld_nontradables, Ld_tradables
end


function eval_prices_direct(sd::SOEdef, val_int_C, G, pN, bv, μv, σv, ξv, ζv, zv)
	pars = sd.pars
	Ls = 1.0
	
	# Step 1: Find unconstrained wage
	res = Optim.optimize(
		w -> (Ls - sum(labor_demand(sd, w, zv, ζv, pN)))^2,
		min(0.5*pN,0.5), max(1,4*pN), GoldenSection()
		)

	# Step 2: Apply constraint
	wv = max(res.minimizer, pars[:wbar])

	# Step 3: Compute labor demand in nontradables and supply
	Ld_N, Ld_T = labor_demand(sd, wv, zv, ζv, pN)
	Ld = Ld_N + Ld_T
	supply_N = TFP_N(zv,pars[:Δ],ζv) * Ld_N^pars[:α_N]
	supply_T = TFP_T(zv,pars[:Δ],ζv) * Ld_T^pars[:α_T]

	# Step 4: Get nontraded demand
 	# Recover traded demand from total consumption
 	ϖ, η, ϑ = [pars[key] for key in [:ϖ, :η, :ϑ]]
	pC = price_index(sd, pN)
	cT = val_int_C * (1-ϖ) * (1.0/pC)^(-η)
	cN = val_int_C * ( ϖ ) * (pN/pC)^(-η)
	gN = G * ϑ / pN

	yN = max(0.0, supply_N - gN)

	# Get new implied pN
	pN_new = ϖ^(1/η) / (1-ϖ^(1/η)) * (cT/yN)^(1/η)

	pN_new = min(100.0, pN_new)

	output = pN * supply_N + supply_T
	F = pN - pN_new

	return Dict(:output => output, :Ld => Ld, :pN_new => pN_new, :C => val_int_C, :wage => wv, :F => F, :supply_N => supply_N, :demand_N => cN+gN)
end

function find_prices_direct(sd::SOEdef, val_int_C, G, Bpv, pNg, pNmin, pNmax, bv, μv, σv, ξv, ζv, zv)

	do_solver_prices = true

	if do_solver_prices
		res = Optim.optimize(
			pNv -> (eval_prices_direct(sd, val_int_C, G, pNv, bv, μv, σv, ξv, ζv, zv)[:F])^2,
			pNmin, pNmax, GoldenSection()
			)

		pN = res.minimizer
	else
		if false
			pNg = max(min(pNg, pNmax - 1e-6), pNmin + 1e-6)
			obj_f(x) = (eval_prices_direct(sd, val_int_C, G, x, bv, μv, σv, ξv, ζv, zv)[:F])^2
			res = Optim.optimize(
				pN -> obj_f(first(pN)),
				[pNmin], [pNmax], [pNg], Fminbox(LBFGS())#, autodiff=:forward#, Optim.Options(f_tol=1e-12)
				)
			pN = first(res.minimizer)
		else
			pN = pNg
		end
	end

	vars = eval_prices_direct(sd, val_int_C, G, pN, bv, μv, σv, ξv, ζv, zv)

	pN >= pNmax - 0.05*(pNmax-pNmin) ? exc_dem = 1 : exc_dem = 0
	pN <= pNmin + 0.05*(pNmax-pNmin) ? exc_sup = 1 : exc_sup = 0

	minf = vars[:supply_N] - vars[:demand_N]

	results = [vars[:wage]; pN; vars[:Ld]; vars[:output]]

	return results, [minf], exc_dem, exc_sup
end

function find_all_prices(sd::SOEdef, itp_ϕc)
	eq, gr = sd.eq, sd.gr
	B′_vec = eq[:issuance]
	G_vec  = eq[:spending]

	Jgrid = agg_grid(sd)
	Nj = size(Jgrid, 1)

	results = Array{Float64}(undef, Nj, 4)
	minf	= Array{Float64}(undef, Nj, 1)
	exc_dem = Array{Float64}(undef, Nj)
	exc_sup = Array{Float64}(undef, Nj)

	pN_guess = eq[:pN]
	
	# for js in 1:N
	Threads.@threads for js in 1:Nj
		Bpv = B′_vec[js]
		G = G_vec[js]
		
		pNmin, pNmax = minimum(gr[:pN]), maximum(gr[:pN])
		isnan(pN_guess[js]) ? pNg = (pNmin+pNmax) / 2 : pNg = 1 * pN_guess[js]

		bv = gr[:b][Jgrid[js, 1]]
		μv = gr[:μ][Jgrid[js, 2]]
		σv = gr[:σ][Jgrid[js, 3]]
		ξv = gr[:ξ][Jgrid[js, 4]]
		ζv = gr[:ζ][Jgrid[js, 5]]
		zv = gr[:z][Jgrid[js, 6]]

		val_int_C = get_agg_C(sd, itp_ϕc, bv, μv, σv, ξv, ζv, zv)
		r, mf, ed, es = find_prices_direct(sd, val_int_C, G, Bpv, pNg, pNmin, pNmax, bv, μv, σv, ξv, ζv, zv)

		results[js, :] = r
		minf[js, :], exc_dem[js], exc_sup[js] = mf, ed, es
		if isnan(results[js, 2])
			results[js, 2] = pNg
		end
	end

	return results, minf, exc_dem, exc_sup
end

function get_agg_C(sd::SOEdef, itp_ϕc, bv, μv, σv, ξv, ζv, zv)
	pars, gr = sd.pars, sd.gr
	ωmin_int, ωmax_int = quantile.(LogNormal(μv, σv), [1e-6; 1-1e-6]) .+ pars[:ωmin]
	ωmax_int = min(ωmax_int, maximum(gr[:ω]))
	ωmax = maximum(gr[:ω])
	ωmin_int = pars[:ωmin]
	val_C, sum_prob = 0., 0.
	λϵ = ergodic_ϵ(sd)
	for (jϵ, ϵv) in enumerate(gr[:ϵ])
		f_pdf(ω) = pdf(LogNormal(μv, σv), ω-pars[:ωmin])
		(val_pdf, err) = hquadrature(f_pdf, ωmin_int, ωmax_int, rtol=1e-10, atol=1e-12, maxevals=500)
		sum_prob += val_pdf * λϵ[jϵ]

		f(ω) = f_pdf(ω) * itp_ϕc(ω, ϵv, bv, μv, σv, ξv, ζv, zv)
		(val, err) = hquadrature(f, ωmin_int, ωmax_int, rtol=1e-12, atol=0, maxevals=500)
	
		val_C += val * λϵ[jϵ]
	end

	val_int_C = val_C / sum_prob
end

function update_state_functions!(sd::SOEdef, upd_η; verbose::Bool=false)
	eq, gr = sd.eq, sd.gr
	itp_ϕc = make_itp(sd, sd.ϕ[:c]; agg=false);

	t1 = time()
	results, minf, exc_dem, exc_sup = find_all_prices(sd, itp_ϕc);
	!verbose || print_save(" (new prices in $(time_print(time()-t1)))")

	dist = Array{Float64,1}(undef, 3)
	dist[1] = sqrt.(sum( (results[:, 1] - eq[:wage]).^2 )) / sqrt.(sum(eq[:wage].^2))
	dist[2] = sqrt.(sum( (results[:, 2] - eq[:pN]).^2 ))   / sqrt.(sum(eq[:pN].^2))
	dist[3] = sqrt.(sum( (results[:, 3] - eq[:Ld]).^2 ))   / sqrt.(sum(eq[:Ld].^2))

	eq[:pN] = upd_η * results[:, 2] + (1.0-upd_η) * eq[:pN]

	Jgrid = agg_grid(sd)
	Threads.@threads for js in 1:size(Jgrid,1)
	# for js in 1:N
		Bpv = eq[:issuance][js]
		G = eq[:spending][js]

		pN = eq[:pN][js]
		pNmin, pNmax = minimum(gr[:pN]), maximum(gr[:pN])

		bv = gr[:b][Jgrid[js, 1]]
		μv = gr[:μ][Jgrid[js, 2]]
		σv = gr[:σ][Jgrid[js, 3]]
		ξv = gr[:ξ][Jgrid[js, 4]]
		ζv = gr[:ζ][Jgrid[js, 5]]
		zv = gr[:z][Jgrid[js, 6]]

		pN = max(min(pN, pNmax), pNmin)

		val_int_C = get_agg_C(sd, itp_ϕc, bv, μv, σv, ξv, ζv, zv)

		others = eval_prices_direct(sd, val_int_C, G, pN, bv, μv, σv, ξv, ζv, zv)

		for (key, val) in others # should contain :wage, :Ld, :output, :C
			if haskey(eq, key)
				eq[key][js] = val
			end
		end
	end
	update_fiscalrules!(sd)

	eq[:profits] = eq[:output] - eq[:wage] .* eq[:Ld]
	mean_f = mean(minf)
	max_f = maximum(abs.(minf))

	exc_dem_prop = sum(exc_dem) / length(exc_dem)
	exc_sup_prop = sum(exc_sup) / length(exc_sup)

	return exc_dem_prop, exc_sup_prop, mean_f, max_f, dist
end

function update_grid_p!(sd::SOEdef, exc_dem_prop, exc_sup_prop)
	gr = sd.gr

	pN_down = minimum(gr[:pN])
	if exc_sup_prop > 0.01
		pN_down = pN_down * 0.95
	elseif exc_sup_prop == 0.
		pN_down = pN_down * 1.01
	end
	pN_up = maximum(gr[:pN])
	if exc_dem_prop > 0.01
		pN_up = pN_up * 1.05
	elseif exc_dem_prop == 0.
		pN_up = pN_up * 0.99
	end

	gr[:pN] = collect(range(pN_down, pN_up, length=length(gr[:pN])))

	nothing
end

function find_q(sd::SOEdef, q, a, b, var_a, var_b, cov_ab, Bpv, ξpv, ζpv, zpv, jdef, itp_qᵍ, reentry)
	pars = sd.pars
	haircut = 0.0
	if jdef && reentry==false
		haircut = 0.0
	end
	if jdef == false && ζpv == 0
		haircut = pars[:ℏ]
	end

	R = (ζpv==1) * pars[:κ] + (1.0 - haircut) .* ((1.0-pars[:ρ])*q)

	Eω   = a + R*b
	varω = var_a + R^2 * var_b + 2*R * cov_ab
	varω = max(varω, 0.)

	# varω >= 0. || print_save("\nvar_a, var_b, cov_ab, R, q = $(var_a), $(var_b), $(cov_ab), $(R), $(q)")

	μpv, σpv = make_logN(max(0.0, Eω - pars[:ωmin]), varω)

	itp_qᵍ = extrapolate(itp_qᵍ, Interpolations.Flat())
	new_q = itp_qᵍ((1.0 - haircut) .* Bpv, μpv, σpv, ξpv, ζpv, zpv)

	return new_q, μpv, σpv
end

function compute_stats_logN(sd::SOEdef, a, b, var_a, var_b, cov_ab, itp_qᵍ, Bpv, exp_rep, jdef)
	gr = sd.gr

	q = zeros(N(sd,:ξ), N(sd,:z), 2)
	μ = [zeros(2) for jξp in 1:N(sd,:ξ), jzp in 1:N(sd,:z)]
	σ = [zeros(2) for jξp in 1:N(sd,:ξ), jzp in 1:N(sd,:z)]

	alarm_mat = Array{Float64, 3}(undef, N(sd,:ξ), N(sd,:z), 2)

	for (jξp, ξpv) in enumerate(gr[:ξ]), (jzp, zpv) in enumerate(gr[:z])
		# First any case where ζ′ = 1 -- no default
		reentry = true
		jζp = 2
		ζpv = gr[:ζ][jζp]
		qmin, qmax = minimum(sd.eq[:qᵍ]), maximum(sd.eq[:qᵍ])

		res = Optim.optimize(
			q -> (find_q(sd, q, a, b, var_a, var_b, cov_ab, Bpv, ξpv, ζpv, zpv, jdef, itp_qᵍ, reentry)[1] - q)^2,
			qmin, qmax, Brent()
			)
		q[jξp, jzp, jζp] = res.minimizer
		res.minimum > 1e-4 ? alarm_mat[jξp, jzp, jζp] = 1 : alarm_mat[jξp, jzp, jζp] = 0

		_, μ[jξp, jzp][jζp], σ[jξp, jzp][jζp] = find_q(sd, q[jξp, jzp, jζp], a, b, var_a, var_b, cov_ab, Bpv, ξpv, ζpv, zpv, jdef, itp_qᵍ, reentry)

		if jdef
			# If default continues (both ζv and ζvp are 0)
			reentry = false
			jζp = 1
			ζpv = gr[:ζ][jζp] # Irrelevant but to stress that the default state continues
			res = Optim.optimize(
				q -> (find_q(sd, q, a, b, var_a, var_b, cov_ab, Bpv, ξpv, ζpv, zpv, jdef, itp_qᵍ, reentry)[1] - q)^2,
				qmin, qmax, Brent()
				)
			q[jξp, jzp, jζp] = res.minimizer
			res.minimum > 1e-4 ? alarm_mat[jξp, jzp, jζp] = 1 : alarm_mat[jξp, jzp, jζp] = 0

			_, μ[jξp, jzp][jζp], σ[jξp, jzp][jζp] = find_q(sd, q[jξp, jzp, jζp], a, b, var_a, var_b, cov_ab, Bpv, ξpv, ζpv, zpv, jdef, itp_qᵍ, reentry)
		else
			# Entering default
			reentry = true # Irrelevant
			jζp = 1
			ζpv = gr[:ζ][jζp]
			res = Optim.optimize(
				q -> (find_q(sd, q, a, b, var_a, var_b, cov_ab, Bpv, ξpv, ζpv, zpv, jdef, itp_qᵍ, reentry)[1] - q)^2,
				qmin, qmax, Brent()
				)
			q[jξp, jzp, jζp] = res.minimizer
			res.minimum > 1e-4 ? alarm_mat[jξp, jzp, jζp] = 1 : alarm_mat[jξp, jzp, jζp] = 0

			_, μ[jξp, jzp][jζp], σ[jξp, jzp][jζp] = find_q(sd, q[jξp, jzp, 2], a, b, var_a, var_b, cov_ab, Bpv, ξpv, ζpv, zpv, jdef, itp_qᵍ, reentry)
		end
	end
	return μ, σ, q, alarm_mat
end

function new_expectations(sd::SOEdef, itp_ϕa, itp_ϕb, itp_qᵍ, Bpv, exp_rep, js, jdef)
	pars, gr = sd.pars, sd.gr

	itp_ϕa = extrapolate(itp_ϕa, Interpolations.Line())
	itp_ϕb = extrapolate(itp_ϕb, Interpolations.Line())

	Jgrid = agg_grid(sd)

	bv = gr[:b][Jgrid[js, 1]]
	μv = gr[:μ][Jgrid[js, 2]]
	σv = gr[:σ][Jgrid[js, 3]]
	ξv = gr[:ξ][Jgrid[js, 4]]
	ζv = gr[:ζ][Jgrid[js, 5]]
	zv = gr[:z][Jgrid[js, 6]]

	λϵ = ergodic_ϵ(sd)

	val_a, val_b, val_a2, val_b2, val_ab, sum_prob = 0., 0., 0., 0., 0., 0.

	ωmin_int, ωmax_int = quantile.(LogNormal(μv, σv), [1e-6; 1-1e-6]) .+ pars[:ωmin]
	ωmax_int = min(ωmax_int, maximum(gr[:ω]))
	ωmin_int = pars[:ωmin]
	for (jϵ, ϵv) in enumerate(gr[:ϵ])
		f_pdf(ω) = pdf(LogNormal(μv, σv), ω-pars[:ωmin])
		(val_pdf, err) = hquadrature(f_pdf, ωmin_int, ωmax_int, rtol=1e-10, atol=1e-12, maxevals=0)
		sum_prob += val_pdf * λϵ[jϵ]

		fA(ω) = f_pdf(ω) * max(pars[:ωmin], itp_ϕa(ω, ϵv, bv, μv, σv, ξv, ζv, zv))
		(valA, err) = hquadrature(fA, ωmin_int, ωmax_int, rtol=1e-10, atol=1e-12, maxevals=0)
		val_a += valA * λϵ[jϵ]

		fA2(ω) = f_pdf(ω) * max(pars[:ωmin], itp_ϕa(ω, ϵv, bv, μv, σv, ξv, ζv, zv))^2
		(valA2, err) = hquadrature(fA2, ωmin_int, ωmax_int, rtol=1e-10, atol=1e-12, maxevals=0)
		val_a2 += valA2 * λϵ[jϵ]

		fB(ω) = f_pdf(ω) * max(0., itp_ϕb(ω, ϵv, bv, μv, σv, ξv, ζv, zv))
		(valB, err) = hquadrature(fB, ωmin_int, ωmax_int, rtol=1e-10, atol=1e-12, maxevals=0)
		val_b += valB * λϵ[jϵ]

		fB2(ω) = f_pdf(ω) * max(0., itp_ϕb(ω, ϵv, bv, μv, σv, ξv, ζv, zv))^2
		(valB2, err) = hquadrature(fB2, ωmin_int, ωmax_int, rtol=1e-10, atol=1e-12, maxevals=0)
		val_b2 += valB2 * λϵ[jϵ]

		fAB(ω) = f_pdf(ω) * max(pars[:ωmin], itp_ϕa(ω, ϵv, bv, μv, σv, ξv, ζv, zv)) * max(0., itp_ϕb(ω, ϵv, bv, μv, σv, ξv, ζv, zv))
		(valAB, err) = hquadrature(fAB, ωmin_int, ωmax_int, rtol=1e-10, atol=1e-12, maxevals=0)
		val_ab += valAB * λϵ[jϵ]
	end

	!isnan(sum_prob) || throw(error("NaN in sum_prob"))
	!isapprox(sum_prob, 0.) || throw(error("\nsum_prob = $(sum_prob)"))

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

	!isnan(var_a+var_b+cov_ab) || print_save("\nVa, Vb, cov = $var_a, $var_b, $cov_ab")

	μ′, σ′, q, alarm_vec = compute_stats_logN(sd, a, b, var_a, var_b, cov_ab, itp_qᵍ, Bpv, exp_rep, jdef)
	return μ′, σ′, alarm_vec
end

function find_all_expectations(sd::SOEdef, itp_ϕa, itp_ϕb, itp_qᵍ)
	B′_vec = sd.eq[:issuance]

	Jgrid = agg_grid(sd)
	Nj = size(Jgrid, 1)

	# μ′ = Array{Float64}(undef, Nj, N(sd,:ξ), N(sd,:z), 2)
	# σ′ = Array{Float64}(undef, Nj, N(sd,:ξ), N(sd,:z), 2)
	μ′ = similar(sd.LoM[:μ])
	σ′ = similar(sd.LoM[:σ])
	alarm_vec = Array{Float64}(undef, Nj, N(sd,:ξ), N(sd,:z), 2)

	# Threads.@threads for js in 1:N
	for js in 1:Nj
		Bpv = B′_vec[js]

		jb = Jgrid[js, 1]
		jμ = Jgrid[js, 2]
		jσ = Jgrid[js, 3]
		jξ = Jgrid[js, 4]
		jζ = Jgrid[js, 5]
		jz = Jgrid[js, 6]

		rep_mat = reshape(sd.gov[:repay], N(sd,:b), N(sd,:μ), N(sd,:σ), N(sd,:ξ), N(sd,:ζ), N(sd,:z), N(sd,:ξ), N(sd,:z))
		exp_rep = rep_mat[jb, jμ, jσ, jξ, jζ, jz, :, :]

		jdefault = (jζ == 1)

		μ′[js, :, :], σ′[js, :, :], alarm_vec[js, :, :, :] = new_expectations(sd, itp_ϕa, itp_ϕb, itp_qᵍ, Bpv, exp_rep, js, jdefault)
	end

	if sum(alarm_vec[:]) >= 1
		print_save("\nWARNING: Couldn't find qᵍ $(round(100*sum(alarm_vec[:])/Nj,digits=0))% of the time")
	end

	return μ′, σ′
end

fl(x) = [x[js][jj] for js in 1:length(x), jj in 1:2][:]
function new_grid(x′, xgrid; ub=Inf, lb=-Inf)
	Nx = length(xgrid)
	xdist = maximum(xgrid) - minimum(xgrid)

	Xmax = maximum(xgrid)
	Xmin = minimum(xgrid)


	# Expand grids if x′ goes beyond the bounds
	quantile(fl(x′)[:], 0.99) > maximum(xgrid) ? Xmax = maximum(xgrid) + 0.01*xdist : nothing
	quantile(fl(x′)[:], 0.01) < minimum(xgrid) ? Xmin = minimum(xgrid) - 0.01*xdist : nothing

	# Retract grids if x′ doesn't reach the bounds
	maximum(fl(x′)) < maximum(xgrid) ? Xmax = maximum(xgrid) - 0.01*xdist : nothing
	minimum(fl(x′)) > minimum(xgrid) ? Xmin = minimum(xgrid) + 0.01*xdist : nothing

	Xmax = min(Xmax, ub)
	Xmin = max(Xmin, lb)

	return collect(range(Xmin, Xmax, length=Nx))
end

function update_expectations!(sd::SOEdef, upd_η::Float64)
	"""
	Computes mean and variance of tomorrow's distribution and deduces parameters for logN
	"""
	LoM, gr = sd.LoM, sd.gr

	μ_old = copy(LoM[:μ])
	σ_old = copy(LoM[:σ])

	dist_exp = Array{Float64,1}(undef, 2)

	itp_ϕa = make_itp(sd, sd.ϕ[:a]; agg=false)
	itp_ϕb = make_itp(sd, sd.ϕ[:b]; agg=false)
	itp_qᵍ = make_itp(sd, sd.eq[:qᵍ]; agg=true)

	μ_new, σ_new = find_all_expectations(sd, itp_ϕa, itp_ϕb, itp_qᵍ)

	new_μgrid = new_grid(μ_new, gr[:μ], lb = -3.0, ub = 3.0)
	new_σgrid = new_grid(σ_new, gr[:σ], lb = 1e-2)

	for js in 1:length(μ_new)
		μ_new[js] = max.(min.(μ_new[js], maximum(gr[:μ])), minimum(gr[:μ]))
		σ_new[js] = max.(min.(σ_new[js], maximum(gr[:σ])), minimum(gr[:σ]))
	end

	dist_exp[1] = sqrt.(sum( (fl(μ_new) - fl(μ_old)).^2 )) / sqrt.(sum(fl(μ_old).^2))
	dist_exp[2] = sqrt.(sum( (fl(σ_new) - fl(σ_old)).^2 )) / sqrt.(sum(fl(σ_old).^2))

	μ_new = upd_η * μ_new + (1.0 - upd_η) * μ_old
	σ_new = upd_η * σ_new + (1.0 - upd_η) * σ_old

	LoM[:μ] = μ_new
	LoM[:σ] = σ_new

	return dist_exp, new_μgrid, new_σgrid
end

function update_grids!(sd::SOEdef; new_μgrid::Vector=[], new_σgrid::Vector=[], new_zgrid::Vector=[])
	gr = sd.gr

	if new_μgrid==[]
		new_μgrid = gr[:μ]
	end
	if new_σgrid==[]
		new_σgrid = gr[:σ]
	end
	if new_zgrid==[]
		new_zgrid = gr[:z]
	end

	function reinterp(sd::SOEdef, y; agg::Bool=false, ext::Bool=false)
		gr = sd.gr
		knots = (gr[:ω], gr[:ϵ], gr[:b], gr[:μ], gr[:σ], gr[:ξ], gr[:ζ], gr[:z])
		if agg
			knots = (gr[:b], gr[:μ], gr[:σ], gr[:ξ], gr[:ζ], gr[:z])
			y = reshape(y, N(sd,:b), N(sd,:μ), N(sd,:σ), N(sd,:ξ), N(sd,:ζ), N(sd,:z))
		end
		if ext
			knots = (gr[:ω], gr[:ϵ], gr[:b], gr[:μ], gr[:σ], gr[:ξ], gr[:ζ], gr[:z], gr[:pN])
		end

		itp_obj_y = interpolate(knots, y, Gridded(Linear()))
		itp_y = extrapolate(itp_obj_y, Line())

		if agg
			y_new = itp_y(gr[:b], new_μgrid, new_σgrid, gr[:ξ], gr[:ζ], new_zgrid)
			return reshape(y_new, length(y_new))
		elseif ext
			y_new = itp_y(gr[:ω], gr[:ϵ], gr[:b], new_μgrid, new_σgrid, gr[:ξ], gr[:ζ], new_zgrid, gr[:pN])
			return y_new
		else
			y_new = itp_y(gr[:ω], gr[:ϵ], gr[:b], new_μgrid, new_σgrid, gr[:ξ], gr[:ζ], new_zgrid)
			return y_new
		end
	end

	for key in [:a, :b, :c, :s, :θ]
		sd.ϕ[key] = reinterp(sd, sd.ϕ[key], agg=false)
	end
	sd.v[:v] = max.(1e-20, reinterp(sd, sd.v[:v], agg=false))
	sd.v[:w] = max.(1e-20, reinterp(sd, sd.v[:w], agg=false))

	for key in [:Ld, :output, :wage, :spending, :pN, :issuance, :welfare]
		sd.eq[key] = reinterp(sd, sd.eq[key], agg=true)
	end
	sd.eq[:issuance] = min.(max.(sd.eq[:issuance], minimum(gr[:b])), maximum(gr[:b]))

	knots 		= (gr[:b], gr[:μ], gr[:σ], gr[:ξ], gr[:ζ], gr[:z], gr[:ξ], gr[:z])
	repay_mat 	= reshape(sd.gov[:repay], N(sd,:b), N(sd,:μ), N(sd,:σ), N(sd,:ξ), N(sd,:ζ), N(sd,:z), N(sd,:ξ), N(sd,:z))
	itp_repay 	= extrapolate(interpolate(knots, repay_mat, Gridded(Linear())), Line())
	rep_new 	= itp_repay(gr[:b], new_μgrid, new_σgrid, gr[:ξ], gr[:ζ], new_zgrid, gr[:ξ], new_zgrid)
	rep_new 	= max.(0, min.(1, rep_new))
	sd.gov[:repay] 	= reshape(rep_new, length(rep_new))

	μmax, μmin = maximum(new_μgrid) - 1e-6, minimum(new_μgrid) + 1e-6
	σmax, σmin = maximum(new_σgrid) - 1e-6, minimum(new_σgrid) + 1e-6

	μ′_new = similar(sd.LoM[:μ])
	σ′_new = similar(sd.LoM[:σ])

	for (jξp, ξpv) in enumerate(gr[:ξ])
		for jzp in 1:length(new_zgrid)
			try
				mvec = [reinterp(sd, [sd.LoM[:μ][js,jξp,jzp][jre] for js in 1:size(sd.LoM[:μ],1)], agg=true) for jre in 1:2]
				svec = [reinterp(sd, [sd.LoM[:σ][js,jξp,jzp][jre] for js in 1:size(sd.LoM[:σ],1)], agg=true) for jre in 1:2]
				for js in 1:size(μ′_new,1)
					μ′_new[js,jξp,jzp] = [max(min(mvec[jre][js], μmax), μmin) for jre in 1:2]
					σ′_new[js,jξp,jzp] = [max(min(svec[jre][js], σmax), σmin) for jre in 1:2]
				end
			catch
				mvec = [reinterp(sd, [sd.LoM[:μ][js,jξp,N(sd,:z)][jre] for js in 1:size(sd.LoM[:μ],1)], agg=true) for jre in 1:2]
				svec = [reinterp(sd, [sd.LoM[:σ][js,jξp,N(sd,:z)][jre] for js in 1:size(sd.LoM[:σ],1)], agg=true) for jre in 1:2]
				for js in 1:size(μ′_new,1)
					μ′_new[js,jξp,jzp] = [max(min(mvec[jre][js], μmax), μmin) for jre in 1:2]
					σ′_new[js,jξp,jzp] = [max(min(svec[jre][js], σmax), σmin) for jre in 1:2]
				end
			end
		end
	end

	gr[:μ] = new_μgrid
	gr[:σ] = new_σgrid

	sd.LoM[:μ] = μ′_new
	sd.LoM[:σ] = σ′_new

	nothing
end

function comp_eqm!(sd::SOEdef; tol::Float64=5e-3, maxiter::Int64=2500, verbose::Bool=false, iter_show::Int64 = 50)
	dist = 1+tol
	iter = 0

	upd_η = 0.5

	tol_vfi = 1e-2
	t0 = time()
	while dist > tol && iter < maxiter
		iter += 1
		if verbose
			print_save("\nIteration $iter")
			print_save("(vfi update tolerance = $(@sprintf("%0.3g",tol_vfi)))")
			print_save(". (upd_η = $(@sprintf("%0.3g", upd_η)))")
			print_save(Dates.format(now(), "HH:MM"))
		end

		iterate_qᵍ!(sd, verbose = verbose)
		update_fiscalrules!(sd)

		Jgrid = agg_grid(sd);
		var(sd.eq[:qʰ]) .< 1e-16 || print_save("\nWARNING: qʰ is not constant. $(var(sd.eq[:qʰ]))")
		!verbose || print_save("\nqᵍ between $(round(minimum(sd.eq[:qᵍ][Jgrid[:,5].==1]),digits=4)) and $(round(maximum(sd.eq[:qᵍ]),digits=4)). risk-free is $(round(mean(sd.eq[:qʰ]),digits=4))")
		!verbose || print_save(" (spread between $(floor(Int,10000*minimum(sd.eq[:spread]))) bps and $(floor(Int,10000*maximum(sd.eq[:spread][Jgrid[:,5].==1]))) bps)")

		""" SOLVE INCOME FLUCTUATIONS PROBLEM """
		consw, dist_v = vfi!(sd, tol = tol_vfi, verbose = false);
		flag = (dist_v < tol_vfi)
		if flag && iter % iter_show == 0
			print_save(" ✓")
		end

		consw = 100*floor(Int,mean(consw))
		consw > 0 ? msg = "\nWARNING: " : msg="\n"
		msg *= "Can't affort consumption $(consw)% of the time"
		!verbose || print_save(msg)

		""" UPDATE STATE FUNCTIONS """
		t1 = time()
		!verbose || print_save("\nUpdating functions of the state")
		exc_dem_prop, exc_sup_prop, mean_excS, max_excS, dists = update_state_functions!(sd, upd_η)
		!verbose || print_save(": done in $(time_print(time()-t1))")
		!verbose || print_save("\nStates with exc supply, demand = $(round(100*exc_sup_prop,digits=2))%, $(round(100*exc_dem_prop,digits=2))%")
		!verbose || print_save("\nAverage, max exc supply = $(@sprintf("%0.3g",mean_excS)), $(@sprintf("%0.3g",max_excS))")

		""" UPDATE GRID FOR PN """
		t1 = time()
		update_grid_p!(sd, exc_dem_prop, exc_sup_prop)
		!verbose || print_save("\nNew pN_grid = [$(@sprintf("%0.3g",minimum(sd.gr[:pN]))), $(@sprintf("%0.3g",maximum(sd.gr[:pN])))]")
		!verbose || print_save("\nDistance in state functions: (dw,dpN,dLd) = ($(@sprintf("%0.3g",mean(dists[1]))),$(@sprintf("%0.3g",mean(dists[2]))),$(@sprintf("%0.3g",mean(dists[3]))))")
		
		dist_s = maximum(dists)

		""" UPDATE EXPECTATIONS AND GRIDS FOR LOMS """
		dist_exp, new_μgrid, new_σgrid = update_expectations!(sd, 0.5 * upd_η)
		update_grids!(sd, new_μgrid = new_μgrid, new_σgrid = new_σgrid)

		!verbose || print_save("\nDistance in expectations: (dμ,dσ) = ($(@sprintf("%0.3g",mean(dist_exp[1]))),$(@sprintf("%0.3g",mean(dist_exp[2]))))")
		!verbose || print_save("\nNew μ_grid = [$(@sprintf("%0.3g",minimum(sd.gr[:μ]))), $(@sprintf("%0.3g",maximum(sd.gr[:μ])))]")
		!verbose || print_save("\nNew σ_grid = [$(@sprintf("%0.3g",minimum(sd.gr[:σ]))), $(@sprintf("%0.3g",maximum(sd.gr[:σ])))]")
		# !verbose || print_save("\nNew σ′ ∈ [$(@sprintf("%0.3g",minimum(fl(sd.LoM[:σ])))), $(@sprintf("%0.3g",maximum(fl(sd.LoM[:σ]))))]")
		!verbose || print_save("\nGrids and expectations updated in $(time_print(time()-t1))")

		dist_exp = maximum(dist_exp)
		if iter % iter_show == 0
			print_save("\nDistances: (10dv, dLoM, dp) = ($(@sprintf("%0.3g",minimum(10*dist_v))), $(@sprintf("%0.3g",minimum(dist_exp))), $(@sprintf("%0.3g",minimum(dist_s))))")
		end
		dist_s = max(dist_s, dist_exp)
		dist = max(10*dist_v, dist_s)
		tol_vfi = max(exp(0.95*log(1+tol_vfi))-1, tol/10)

		upd_η = max(0.9*upd_η, 5e-2)
	end

	if dist <= tol
		print_save("\nConverged in $iter iterations. ")
	else
		print_save("\nStopping at distance $(@sprintf("%0.3g",dist)), target was $(@sprintf("%0.3g",tol)).")
	end

	print_save("\nTime: $(time_print(time()-t0))")
	nothing
end
