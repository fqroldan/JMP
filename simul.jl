function initial_dist(sd::SOEdef, μ0=mean(sd.gr[:μ]), σ0=mean(sd.gr[:σ]))
	λϵ = ergodic_ϵ(sd)	
	λ = zeros(N(sd,:ωf) * N(sd,:ϵ))
	js = 0
	for (jϵ, ϵv) in enumerate(sd.gr[:ϵ])
		for jω = 1:length(sd.gr[:ωf])-1
			js += 1
			ωv  = sd.gr[:ωf][jω]
			ω1v = sd.gr[:ωf][jω+1]
			ωmv = 0.5*(ωv+ω1v)

			λ[js] = pdf(LogNormal(μ0, σ0), ωmv-sd.pars[:ωmin]) * λϵ[jϵ] * (ω1v - ωv)
		end
	end
	λ = λ / sum(λ)
	return λ
end

function integrate_C(sd::SOEdef, Bt, μt, σt, ξt, ζt, zt, λt, itp_ϕc, itp_C)
	C_from_interp = itp_C(Bt, μt, σt, ξt, ζt, zt)

	C_from_λ = 0.0
	λ_mat  = reshape(λt, N(sd,:ωf), N(sd,:ϵ))
	for (jω, ωv) in enumerate(sd.gr[:ωf]), (jϵ, ϵv) in enumerate(sd.gr[:ϵ])
		C_from_λ += itp_ϕc(ωv, ϵv, Bt, μt, σt, ξt, ζt, zt) * λ_mat[jω, jϵ]
	end

	return C_from_λ, C_from_interp - C_from_λ
end


function iter_simul!(sd::SOEdef, p::Path, t, jz_series, itp_ϕa, itp_ϕb, itp_ϕc, itp_ϕs, itp_ϕθ, itp_vf, itp_C, itp_B′, itp_G, itp_pN, itp_qᵍ, itp_repay, itp_W, λt, Qϵ)


	# Enter with a state B, μ, σ, w0, ζ, z.
	# h.zgrid[jz] must equal getfrompath(p, t, :z)
	# B, ζ, and z are decided at the end of the last period
	jz = jz_series[t]

	Bt = getfrompath(p, t, :B)
	μt = getfrompath(p, t, :μ)
	σt = getfrompath(p, t, :σ)
	ξt = getfrompath(p, t, :ξ)
	jξ = findfirst(sd.gr[:ξ] .== ξt)
	ζt = getfrompath(p, t, :ζ)
	zt = getfrompath(p, t, :z)

	zt == sd.gr[:z][jz] || print("something wrong with ztv and jzt")
	abs(zt - sd.gr[:z][jz]) < 1e-8 || throw(error("something wrong with ztv and jzt"))

	# print("\n$(t), μ = $μt")
	if t % 100 == 0
		print_save("\n[Bt, μt, σt, ξt, ζt, zt]\n")
		print_save("$(@sprintf("%0.3g",Bt)), ")
		print_save("$(@sprintf("%0.3g",μt)), ")
		print_save("$(@sprintf("%0.3g",σt)), ")
		print_save("$(@sprintf("%0.3g",ξt)), ")
		print_save("$(@sprintf("%0.3g",ζt)), ")
		print_save("$(@sprintf("%0.3g",zt))")
	end

	Bpv 	= itp_B′(Bt, μt, σt, ξt, ζt, zt)
	lb = sd.gr[:b][end] - sd.gr[:b][1]
	Bpv = max.(min.(Bpv, sd.gr[:b][end]-0.01*lb), sd.gr[:b][1]+0.01*lb)

	Gt 		= itp_G(Bt, μt, σt, ξt, ζt, zt)
	pNg 	= itp_pN(Bt, μt, σt, ξt, ζt, zt)

	if ζt == 0 # Default when ζ == 0., jζ == 1
		Bpv = Bt
	end

	exp_rep = zeros(N(sd,:ξ), N(sd,:z))
	# itp_repay = extrapolate(itp_repay, Interpolations.Flat())
	for (jξp, ξpv) in enumerate(sd.gr[:ξ]), (jzp, zpv) in enumerate(sd.gr[:z])
		exp_rep[jξp, jzp] = max(0, min(1, itp_repay(Bt, μt, σt, ξt, ζt, zt, ξpv, zpv)))
	end

	# Find pN at the current state. Deduce w, L, Π, T.
	pNmin, pNmax = minimum(sd.gr[:pN]), maximum(sd.gr[:pN])
	jdef = (ζt == 0)

	val_int_C, discrepancy = integrate_C(sd, Bt, μt, σt, ξt, ζt, zt, λt, itp_ϕc, itp_C)

	results, _ = find_prices_direct(sd, val_int_C, Gt, Bpv, pNg, pNmin, pNmax, Bt, μt, σt, ξt, ζt, zt)

	wt, pN, Ld, output = results
	profits = output - wt*Ld

	pCt = price_index(sd, pN)
	Ld_N, Ld_T  = labor_demand(sd, wt, zt, ζt, pN)
	supply_T = TFP_T(zt, sd.pars[:Δ], ζt) * Ld_T^(sd.pars[:α_N])


	# Govt BC
	if ζt == 1
		coupons = sd.pars[:κ] * Bt
		new_debt = Bpv - (1.0 - sd.pars[:ρ]) * Bt
	else
		coupons = 0
		new_debt = 0
	end
	qg = max(minimum(sd.eq[:qᵍ]), getfrompath(p, t, :qg))

	lumpsumT = coupons + Gt - sd.pars[:τ]*wt*Ld - qg*new_debt

	
	if t % 100 == 0
		print_save("\npN = $(@sprintf("%0.3g",pN)), pN^e = $(@sprintf("%0.3g",pNg)), u = $(ceil(100*(1-Ld))) at t = $t")
		print_save("\nDiscrepancy in ∫ϕcdμ, pN = $(@sprintf("%0.3g", discrepancy)), $(@sprintf("%0.3g", log(pN) - log(pNg)))")
	end

	def_prob = 0.
	if ζt == 1
		for (jξp, ξvp) in enumerate(sd.gr[:ξ]), (jzp, zvp) in enumerate(sd.gr[:z])
			def_prob += sd.prob[:ξ][jξ, jξp] * sd.prob[:z][jz, jzp] * (1 - exp_rep[jξp, jzp])
		end
	end


	# Integrate the household's policy functions to get μ′, σ′ and C
	ϕa = zeros(N(sd,:ωf)*N(sd,:ϵ))
	ϕb = zeros(N(sd,:ωf)*N(sd,:ϵ))
	ϕc = zeros(N(sd,:ωf)*N(sd,:ϵ))
	Crate = zeros(N(sd,:ωf)*N(sd,:ϵ))

	js = 0
	ωϵv = gridmake(sd.gr[:ωf], sd.gr[:ϵ])
	λ_mat  = reshape(λt, N(sd,:ωf), N(sd,:ϵ))
	mλω = sum(λ_mat, dims=2)
	cdf_ω = cumsum(mλω[:])

	avgω_fromb = zeros(N(sd,:ωf)*N(sd,:ϵ))

	b25 = zeros(N(sd,:ωf)*N(sd,:ϵ))
	b90 = zeros(N(sd,:ωf)*N(sd,:ϵ))

	cdf_ω[1] > 0.25  ? q25 = 1 : q25 = findfirst(cdf_ω .<= 0.25)
	cdf_ω[end] < 0.9 ? q90 = length(cdf_ω) : q90 = findfirst(cdf_ω .>= 0.90)

	qhv = sd.eq[:qʰ][1]
	adjustment = sum(ergodic_ϵ(sd).*exp.(sd.gr[:ϵ]))
	for js in 1:size(ωϵv, 1)
		ωv, ϵv = ωϵv[js, :]
		yd = (wt*Ld*(1.0-sd.pars[:τ]) + profits) * exp(ϵv)/adjustment + ωv - lumpsumT
		
		sg = max(sd.pars[:ωmin], itp_ϕs(ωv, ϵv, Bt, μt, σt, ξt, ζt, zt))
		θg = min(max(0.0, itp_ϕθ(ωv, ϵv, Bt, μt, σt, ξt, ζt, zt)), 1.0)

		abc = get_abc(yd, sd.pars[:ωmin], qhv, qg, pCt, sg, θg)
		ap, bp, cc = abc[:a], abc[:b], abc[:c]

		isnan(ap) ? ap = sd.pars[:ωmin] : nothing
		isnan(bp) ? bp = 0 : nothing

		ϕa[js] = max(sd.pars[:ωmin], ap)
		ϕb[js] = max(0.0, bp)

		ϕc[js] = cc
		Crate[js] = ϕc[js] / yd

		avgω_fromb[js] = ωv * max(0.0, bp)

		jω = js[1]
		if jω < q25
			b25[js] += max(0.0, bp)
		elseif jω > q90
			b90[js] += max(0.0, bp)
		end
	end

	C = dot(λt, ϕc)
	CoY = C / output
	CoYd = dot(λt, Crate)

	p25 = dot(λt, b25) / dot(λt, ϕb)
	p90 = dot(λt, b90) / dot(λt, ϕb)

	aggC_T = C  * sd.pars[:ϖ] * (1 / pCt)^(-sd.pars[:η])
	NX = supply_T - aggC_T - Gt * (1 - sd.pars[:ϑ])

	Eω_fromb = dot(λt, avgω_fromb) / dot(λt, ϕb)

	fill_path!(p,t, Dict(:P => pN, :Pe => pNg, :Y => output, :L => Ld, :π => def_prob, :w => wt, :G => Gt, :CoY=> CoY, :CoYd => CoYd, :C => C, :T => lumpsumT, :NX => NX, :p25 => p25, :p90 => p90, :avgω => Eω_fromb))

	a  = dot(λt, ϕa)
	a2 = dot(λt, ϕa.^2)
	b  = dot(λt, ϕb)
	b2 = dot(λt, ϕb.^2)
	ab = dot(λt, ϕa.*ϕb)

	# print("\na = $a, b = $b, ab = $ab, a2 = $(a2), b2 = $(b2)")
	var_a  = a2 - a^2
	var_b  = b2 - b^2
	cov_ab = ab - a*b

	prop_domestic = b/Bpv
	Bf = Bpv - b

	μ′, σ′, q′, _ = compute_stats_logN(sd, a, b, var_a, var_b, cov_ab, itp_qᵍ, Bpv, exp_rep, jdef)

	# μ′, σ′, q′ = new_expectations(h, itp_ϕa, itp_ϕb, itp_qᵍ, Bpv, wt, thres, Bt, μt, σt, w0, ζt, zt, jdef) # This would assume that λₜ is lognormal
	# print("\n$(q′)")

	# Draw z and the reentry shock for tomorrow, deduce ζ and correct B, μ, and σ as needed, and update the distribution
	probz = cumsum(sd.prob[:z][jz,:])
	jzp = findfirst(probz .> rand())

	probξ = cumsum(sd.prob[:ξ][jξ,:])
	jξp = findfirst(probξ .> rand())

	zpv = sd.gr[:z][jzp]
	ξpv = sd.gr[:ξ][jξp]

	lμ = sd.gr[:μ][end] - sd.gr[:μ][1]
	lσ = sd.gr[:σ][end] - sd.gr[:σ][1]

	μprime = max.(min.(μ′[jξp,jzp], sd.gr[:μ][end]-0.01*lμ), sd.gr[:μ][1]+0.01*lμ)
	σprime = max.(min.(σ′[jξp,jzp], sd.gr[:σ][end]-0.01*lσ), sd.gr[:σ][1]+0.01*lσ)
	qprime = max.(min.(q′[jξp,jzp,:], 1.0), minimum(sd.eq[:qᵍ]))

	if jdef
		# Compute welfare in case of reentry and remain in default
		Wr = itp_W(Bpv, μprime[2], σprime[2], ξpv, sd.gr[:ζ][2], zpv)
		Wd = itp_W(Bpv, μprime[1], σprime[1], ξpv, sd.gr[:ζ][1], zpv)
		# Now draw reentry
		prob_reentry = sd.pars[:θ]
		reentry = (rand() <= prob_reentry)
		if reentry
			jζp = 2
			μpv = μprime[jζp]
			σpv = σprime[jζp]
			qpv = qprime[jζp]
			R = sd.pars[:κ] + (1-sd.pars[:ρ]) * qpv
		else
			jζp = 1
			μpv = μprime[jζp]
			σpv = σprime[jζp]
			qpv = qprime[jζp]
			R = (1-sd.pars[:ρ]) * qpv
		end
	else
		# Compute welfare in case repay and default
		Wr = itp_W(Bpv,					μprime[2], σprime[2], ξpv, sd.gr[:ζ][2], zpv)
		Wd = itp_W((1-sd.pars[:ℏ])*Bpv, μprime[1], σprime[1], ξpv, sd.gr[:ζ][1], zpv)
		# Now draw default
		repay_prime = (rand() <= exp_rep[jξp,jzp])
		if repay_prime
			jζp = 2
			R = sd.pars[:κ] + (1.0-sd.pars[:ρ]) * qprime
		else
			jζp = 1
			Bpv = (1.0 - sd.pars[:ℏ]) * Bpv
			μpv = μprime[jζp]
			σpv = σprime[jζp]
			qpv = qprime[jζp]
			R = (1-sd.pars[:ℏ])*(1-sd.pars[:ρ]) * qpv
		end
	end
	ζpv = sd.gr[:ζ][jζp]

	savings = ϕa + R*ϕb
	savings = max.(min.(savings, sd.pars[:ωmax]), sd.pars[:ωmin])

	basis = Basis(LinParams(sd.gr[:ωf], 0))
	Qω = BasisMatrix(basis, Expanded(), savings, 0).vals[1]
	Q = row_kron(Qϵ, Qω)

	ζt == 1.0 && ζpv == 0.0 ? new_def = 1 : new_def = 0

	λpd = Q' * λt

	M, V = unmake_logN(μpv, σpv)

	# Fill the path for next period
	if t < length(jz_series)
		jz_series[t+1] = jzp
		fill_path!(p,t+1, Dict(:B => Bpv, :μ => μpv, :σ => σpv, :ζ => ζpv, :ξ => ξpv, :z => zpv, :ψ => prop_domestic, :A => a, :Bh => b, :Bf => Bf, :Wr => Wr, :Wd => Wd, :qg => qpv, :mean => M, :var => V))
	end

	return λpd, new_def
end


function simul(sd::SOEdef; simul_length::Int64=1, burn_in::Int64=1)
	gr = sd.gr
	Random.seed!(1)

	# Setup
	T = burn_in + simul_length
	p = Path(T = T)

	jz = 1

	B0, μ0, σ0, ξ0, ζ0, z0 = mean(gr[:b]), mean(gr[:μ]), mean(gr[:σ]), gr[:ξ][1], gr[:ζ][2], gr[:z][jz]
	fill_path!(p,1, Dict(:B => B0, :μ => μ0, :σ => σ0, :w=>1.0, :ξ => ξ0, :ζ => ζ0, :z => z0))

	itp_ϕa = make_itp(sd, sd.ϕ[:a]; agg=false);
	itp_ϕb = make_itp(sd, sd.ϕ[:b]; agg=false);
	itp_ϕc = make_itp(sd, sd.ϕ[:c]; agg=false);
	itp_ϕs = make_itp(sd, sd.ϕ[:s]; agg=false);
	itp_ϕθ = make_itp(sd, sd.ϕ[:θ]; agg=false);
	itp_vf = make_itp(sd, sd.v[:v]; agg=false);

	itp_C  = make_itp(sd, sd.eq[:C]; agg=true);
	itp_B′ = make_itp(sd, sd.eq[:issuance]; agg=true);
	itp_G  = make_itp(sd, sd.eq[:spending]; agg=true);
	itp_pN = make_itp(sd, sd.eq[:pN]; agg=true);
	itp_qᵍ = make_itp(sd, sd.eq[:qᵍ]; agg=true);
	itp_W  = make_itp(sd, sd.eq[:welfare]; agg=true);


	rep_mat = reshape_long_shocks(sd, sd.gov[:repay]);
	knots = (gr[:b], gr[:μ], gr[:σ], gr[:ξ], gr[:ζ], gr[:z], gr[:ξ], gr[:z]);
	itp_repay = interpolate(knots, rep_mat, Gridded(Linear()));

	jz_series = Vector{Int64}(undef, T)
	jz_series[1] = jz

	# Initialize objects for iterating the distribution
	λ = initial_dist(sd, μ0, σ0)
	Qϵ = kron(sd.prob[:ϵ], ones(N(sd,:ωf),1))

	# Simulate
	Ndefs = 0
	for t in 1:T
		λ, new_def = iter_simul!(sd, p, t, jz_series, itp_ϕa, itp_ϕb, itp_ϕc, itp_ϕs, itp_ϕθ, itp_vf, itp_C, itp_B′, itp_G, itp_pN, itp_qᵍ, itp_repay, itp_W, λ, Qϵ)
		Ndefs += new_def
	end

	jz_series = jz_series[burn_in+1:end]

	return trim_path(p, burn_in), jz_series, Ndefs
end

function get_AR1(y::Vector)
	y_lag = y[1:end-1]
	y = y[2:end]

	data = DataFrame(yt = y, ylag = y_lag)
	OLS = glm(@formula(yt ~ ylag - 1), data, Normal(), IdentityLink())

	ρ = coef(OLS)[1]

	ϵ = y - ρ * y_lag

	σ = var(ϵ)^0.5

	return ρ, σ
end

get_MV(y::Vector) = mean(y), var(y)^0.5

function simul_stats(path::Path; nodef::Bool=false, ζ_vec::Vector=[])
	T = periods(path)
	
	if ζ_vec == []
		ζ_vec = series(path,:ζ)
	end

	conditional = (ζ_vec .> -Inf)
	if nodef
		conditional = (ζ_vec .== 1)
	end

	B_vec = series(path,:B)[conditional]
	μ_vec = series(path,:μ)[conditional]
	σ_vec = series(path,:σ)[conditional]
	w_vec = series(path,:w)[conditional]
	z_vec = exp.(series(path,:z)[conditional])
	Y_vec = series(path,:Y)[conditional]
	C_vec = series(path,:C)[conditional]
	G_vec = series(path,:G)[conditional]
	π_vec = series(path,:π)[conditional]
	ψ_vec = series(path,:ψ)[conditional]
	u_vec = 100.0 * (1.0 .- series(path, :L)[conditional])
	spr_vec = 1.0./series(path, :qg)[conditional] .- (1.04)^0.25

	print("\nT = $T")

	m_vec, sd_vec = unmake_logN(μ_vec, σ_vec)
	mean_wealth = 100*mean(m_vec./(4*Y_vec))

	ρy, σy = get_AR1(log.(Y_vec.+1e-8))
	ρc, σc = get_AR1(log.(C_vec.+1e-8))
	if var(spr_vec) > 1e-6
		ρs, σs = get_AR1(spr_vec)
	else
		ρs, σs = 1.0, 0.0
	end
	m_unemp, sd_unemp = get_MV(u_vec)
	m_debt, sd_debt = get_MV(100*B_vec./(4*Y_vec))
	m_gspend, sd_gspend = get_MV(100*G_vec./Y_vec)
	ψ_mean = 100 * median(ψ_vec)
	spr_mean = mean(spr_vec)

	v_m = [ρy; σy; ρc; σc; ρs; σs; m_debt; sd_debt; m_unemp; sd_unemp; ψ_mean; mean_wealth]

	return v_m
end
