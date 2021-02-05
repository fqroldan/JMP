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

	return C_from_λ, log(C_from_interp) - log(C_from_λ)
end

function iter_simul_t!(sd::SOEdef, p::Path, t, jz_series, itp_ϕc, itp_ϕs, itp_ϕθ, itp_C, itp_B′, itp_G, itp_pN, itp_repay, λt, discr, verbose::Bool=false; B2 = -Inf)

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

	# verbose && print("\n$(t), μ = $μt")
	if verbose && t % 10000 == 0
		print_save("\n[Bt, μt, σt, ξt, ζt, zt]\n")
		print_save("$(@sprintf("%0.3g",Bt)), ")
		print_save("$(@sprintf("%0.3g",μt)), ")
		print_save("$(@sprintf("%0.3g",σt)), ")
		print_save("$(@sprintf("%0.3g",ξt)), ")
		print_save("$(@sprintf("%0.3g",ζt)), ")
		print_save("$(@sprintf("%0.3g",zt))")
	end

	lb = sd.gr[:b][end] - sd.gr[:b][1]
	# println(B2)
	if B2 == -Inf
		Bpv = itp_B′(Bt, μt, σt, ξt, ζt, zt)
		Bpv = max.(min.(Bpv, sd.gr[:b][end]-0.01*lb), sd.gr[:b][1]+0.01*lb)
	else
		Bpv = B2
	end

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
	pCC = val_int_C * price_index(sd, pNg)

	discr[:C] = (abs(discrepancy) + (t-1) * discr[:C]) / t

	results, _ = find_prices_direct(sd, pCC, Gt, Bpv, pNg, pNmin, pNmax, Bt, μt, σt, ξt, ζt, zt)

	wt, pN, Ld, output = results
	profits = output - wt*Ld

	discr[:pN] = (abs(log(pN) - log(pNg)) + (t-1) * discr[:pN]) / t

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

	
	if verbose && t % 10000 == 0
		verbose && print_save("\npN = $(@sprintf("%0.3g",pN)), pN^e = $(@sprintf("%0.3g",pNg)), u = $(ceil(100*(1-Ld))) at t = $t")
		verbose && print_save("\nDiscrepancy in ∫ϕcdμ, pN = $(@sprintf("%0.3g", discrepancy)), $(@sprintf("%0.3g", log(pN) - log(pNg)))")
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

	quantiles_vec = [0.1, 0.25, 0.5, 0.75, 0.9]
	quantiles_ω   = [ifelse(cdf_ω[end] > quantile, findfirst(cdf_ω.>=quantile), length(cdf_ω)) for quantile in quantiles_vec]

	C_dist = zeros(length(quantiles_vec))

	sav_a_ω = zeros(length(quantiles_vec))
	sav_b_ω = zeros(length(quantiles_vec))
	prob_ϵω = zeros(length(quantiles_vec), N(sd,:ϵ))

	adjustment = sum(ergodic_ϵ(sd).*exp.(sd.gr[:ϵ]))
	qhv = sd.eq[:qʰ][1]
	for (jq, quantile) in enumerate(quantiles_ω)
		dist_ϵ = λ_mat[quantile, :]
		ωv = sd.gr[:ωf][quantile]
		sav_a = 0.0
		sav_b = 0.0
		cdist = 0.0
		for (jϵ, ϵv) in enumerate(sd.gr[:ϵ])
			yd = (wt*Ld*(1.0-sd.pars[:τ]) + profits) * exp(ϵv)/adjustment + ωv - lumpsumT
			sg = max(sd.pars[:ωmin], itp_ϕs(ωv, ϵv, Bt, μt, σt, ξt, ζt, zt))
			θg = min(max(0.0, itp_ϕθ(sg, ϵv, Bt, μt, σt, ξt, ζt, zt)), 1.0)
			abc = get_abc(yd, sd.pars[:ωmin], qhv, qg, pCt, sg, θg)
			ap, bp, cc = abc[:a], abc[:b], abc[:c]

			sav_a += ap * dist_ϵ[jϵ] / sum(dist_ϵ)
			sav_b += bp * dist_ϵ[jϵ] / sum(dist_ϵ)

			cdist += cc * dist_ϵ[jϵ] / sum(dist_ϵ)
		end

		C_dist[jq] = cdist
		sav_a_ω[jq] = sav_a
		sav_b_ω[jq] = sav_b
		prob_ϵω[jq,:] = dist_ϵ / sum(dist_ϵ)
	end

	avgω_fromb = zeros(N(sd,:ωf)*N(sd,:ϵ))

	b25 = zeros(N(sd,:ωf)*N(sd,:ϵ))
	b90 = zeros(N(sd,:ωf)*N(sd,:ϵ))

	cdf_ω[1] > 0.25  ? q25 = 1 : q25 = findfirst(cdf_ω .>= 0.25)
	cdf_ω[end] < 0.9 ? q90 = length(cdf_ω) : q90 = findfirst(cdf_ω .>= 0.90)

	for js in 1:size(ωϵv, 1)
		ωv, ϵv = ωϵv[js, :]
		yd = (wt*Ld*(1.0-sd.pars[:τ]) + profits) * exp(ϵv)/adjustment + ωv - lumpsumT
		
		sg = max(sd.pars[:ωmin], itp_ϕs(ωv, ϵv, Bt, μt, σt, ξt, ζt, zt))
		θg = min(max(0.0, itp_ϕθ(sg, ϵv, Bt, μt, σt, ξt, ζt, zt)), 1.0)

		abc = get_abc(yd, sd.pars[:ωmin], qhv, qg, pCt, sg, θg)
		ap, bp, cc = abc[:a], abc[:b], abc[:c]

		isnan(ap) ? ap = sd.pars[:ωmin] : nothing
		isnan(bp) ? bp = 0 : nothing

		ϕa[js] = max(sd.pars[:ωmin], ap)
		ϕb[js] = max(0.0, bp)

		ϕc[js] = cc
		Crate[js] = ϕc[js] / yd

		avgω_fromb[js] = ωv * max(0.0, bp)

		if ωv < cdf_ω[q25]
			b25[js] += max(0.0, bp)
		elseif ωv > cdf_ω[q90]
			b90[js] += max(0.0, bp)
		end
	end

	# Compute Gini index
	# S_i = pdf times wealth from 1 to i
	S = [sum( mλω[jj]*ωϵv[jj, 1] for jj in 1:i ) for i in 1:length(mλω)]

	# S_0 = 0
	S = [0;S]

	Gini = 1 - (sum( mλω[jj] * (S[jj] + S[jj+1]) for jj in 1:length(mλω) ))/S[end]

	C = dot(λt, ϕc)
	CoY = C / output
	CoYd = dot(λt, Crate)

	p25 = dot(λt, b25) / dot(λt, ϕb)
	p90 = dot(λt, b90) / dot(λt, ϕb)

	aggC_T = C  * (1-sd.pars[:ϖ]) * pCt^sd.pars[:η] # = (1 / pCt)^(-sd.pars[:η])
	NX = supply_T - aggC_T - Gt * (1 - sd.pars[:ϑ])

	Eω_fromb = dot(λt, avgω_fromb) / dot(λt, ϕb)

	meanω = dot(λt, ωϵv[:,1])
	mω2 = dot(λt, (ωϵv[:,1]).^2)
	varω = mω2 - meanω^2

	fill_path!(p,t, Dict(:P => pN, :Pe => pNg, :Y => output, :L => Ld, :π => def_prob, :w => wt, :G => Gt, :CoY=> CoY, :CoYd => CoYd, :C => C, :T => lumpsumT, :NX => NX, :p25 => p25, :p90 => p90, :mean => meanω, :var => varω, :avgω => Eω_fromb, :Gini => Gini, :C10 => C_dist[1], :C25 => C_dist[2], :C50 => C_dist[3], :C75 => C_dist[4], :C90 => C_dist[5]))
	
	return ϕa, ϕb, Bpv, exp_rep, quantiles_ω, sav_a_ω, sav_b_ω, prob_ϵω, jdef
end

function iter_simul_tp!(sd::SOEdef, p::Path, t, jz_series, λt, ϕa, ϕb, Bpv, quantiles_ω, sav_a_ω, sav_b_ω, prob_ϵω, jdef, itp_repay, itp_qᵍ, itp_vf, Qϵ)

	Bt = getfrompath(p, t, :B)
	μt = getfrompath(p, t, :μ)
	σt = getfrompath(p, t, :σ)
	ξt = getfrompath(p, t, :ξ)
	jξ = findfirst(sd.gr[:ξ] .== ξt)
	ζt = getfrompath(p, t, :ζ)
	zt = getfrompath(p, t, :z)

	jz = findfirst(sd.gr[:z] .== zt)

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

	qguess = itp_qᵍ(Bt, μt, σt, ξt, ζt, zt)

	exp_rep = zeros(N(sd,:ξ), N(sd,:z))
	# itp_repay = extrapolate(itp_repay, Interpolations.Flat())
	for (jξp, ξpv) in enumerate(sd.gr[:ξ]), (jzp, zpv) in enumerate(sd.gr[:z])
		exp_rep[jξp, jzp] = max(0, min(1, itp_repay(Bt, μt, σt, ξt, ζt, zt, ξpv, zpv)))
	end

	μ′, σ′, q′, _ = compute_stats_logN(sd, a, b, var_a, var_b, cov_ab, itp_qᵍ, Bpv, exp_rep, jdef, qguess)

	# μ′, σ′, q′ = new_expectations(h, itp_qᵍ, Bpv, wt, thres, Bt, μt, σt, w0, ζt, zt, jdef) # This would assume that λₜ is lognormal
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
		Bpvec = [Bpv, Bpv]
		repay_prime = (rand() <= sd.pars[:θ])
		Rvec = [qprime[1], sd.pars[:κ] + (1-sd.pars[:ρ]) * qprime[2]]
	else
		Bpvec = [(1-sd.pars[:ℏ])*Bpv, Bpv]
		repay_prime = (rand() <= exp_rep[jξp,jzp])
		Rvec = [(1-sd.pars[:ℏ]) * qprime[1], sd.pars[:κ] + (1.0-sd.pars[:ρ]) * qprime[2]]
	end

	Wr = 0.0 # itp_W(Bpvec[2], μprime[2], σprime[2], ξpv, sd.gr[:ζ][2], zpv)
	Wd = 0.0 # itp_W(Bpvec[1], μprime[1], σprime[1], ξpv, sd.gr[:ζ][1], zpv)

	Wr_vec = [sum(itp_vf(sav_a_ω[jq] + Rvec[2]*sav_b_ω[jq], ϵpv, Bpvec[2], μprime[2], σprime[2], ξpv, sd.gr[:ζ][2], zpv) * prob_ϵω[jq,jϵp] for (jϵp, ϵpv) in enumerate(sd.gr[:ϵ])) for jq in 1:length(quantiles_ω)]
	Wd_vec = [sum(itp_vf(sav_a_ω[jq] + Rvec[1]*sav_b_ω[jq], ϵpv, Bpvec[1], μprime[1], σprime[1], ξpv, sd.gr[:ζ][1], zpv) * prob_ϵω[jq,jϵp] for (jϵp, ϵpv) in enumerate(sd.gr[:ϵ])) for jq in 1:length(quantiles_ω)]

	λpd = zeros(size(λt))
	for jζp in 1:2
		savings = ϕa + Rvec[jζp]*ϕb
		savings = max.(min.(savings, sd.pars[:ωmax]), sd.pars[:ωmin])

		basis = Basis(LinParams(sd.gr[:ωf], 0))
		Qω = BasisMatrix(basis, Expanded(), savings, 0).vals[1]
		Q = row_kron(Qϵ, Qω)

		λp = Q' * λt
		λ_matp = reshape(λp, N(sd,:ωf), N(sd,:ϵ))
		if jζp == 2 # repayment
			Wr = sum(itp_vf(ωpv, ϵpv, Bpvec[jζp], μprime[jζp], σprime[jζp], ξpv, sd.gr[:ζ][jζp], zpv) * λ_matp[jωp, jϵp] for (jωp, ωpv) in enumerate(sd.gr[:ωf]), (jϵp, ϵpv) in enumerate(sd.gr[:ϵ]))
			# println(Wr)
		else
			Wd = sum(itp_vf(ωpv, ϵpv, Bpvec[jζp], μprime[jζp], σprime[jζp], ξpv, sd.gr[:ζ][jζp], zpv) * λ_matp[jωp, jϵp] for (jωp, ωpv) in enumerate(sd.gr[:ωf]), (jϵp, ϵpv) in enumerate(sd.gr[:ϵ]))
			# println(Wd)
		end

		if repay_prime && jζp == 2
			λpd = copy(λp)
		elseif !repay_prime && jζp == 1
			λpd = copy(λp)
		end
	end

	if repay_prime
		jζp = 2
	else
		jζp = 1
	end

	R = Rvec[jζp]
	μpv = μprime[jζp]
	σpv = σprime[jζp]
	qpv = qprime[jζp]

	spr = get_spr(qpv, sd.pars[:κ])

	if !jdef && !repay_prime
		Bpv = (1.0 - sd.pars[:ℏ]) * Bpv
	end


	# if jdef
	# 	# Compute welfare in case of reentry and remain in default
	# 	Wr = itp_W(Bpv, μprime[2], σprime[2], ξpv, sd.gr[:ζ][2], zpv)
	# 	Wd = itp_W(Bpv, μprime[1], σprime[1], ξpv, sd.gr[:ζ][1], zpv)
		
	# 	# Now draw reentry
	# 	reentry = (rand() <= sd.pars[:θ])
	# 	if reentry
	# 		jζp = 2
	# 		μpv = μprime[jζp]
	# 		σpv = σprime[jζp]
	# 		qpv = qprime[jζp]
	# 		R = sd.pars[:κ] + (1-sd.pars[:ρ]) * qpv
	# 	else
	# 		jζp = 1
	# 		μpv = μprime[jζp]
	# 		σpv = σprime[jζp]
	# 		qpv = qprime[jζp]
	# 		R = (1-sd.pars[:ρ]) * qpv
	# 	end
	# else
	# 	# Compute welfare in case repay and default
	# 	Wr = itp_W(Bpv,					μprime[2], σprime[2], ξpv, sd.gr[:ζ][2], zpv)
	# 	Wd = itp_W((1-sd.pars[:ℏ])*Bpv, μprime[1], σprime[1], ξpv, sd.gr[:ζ][1], zpv)
	# 	# Now draw default
	# 	repay_prime = (rand() <= exp_rep[jξp,jzp])
	# 	if repay_prime
	# 		jζp = 2
	# 		μpv = μprime[jζp]
	# 		σpv = σprime[jζp]
	# 		qpv = qprime[jζp]
	# 		R = sd.pars[:κ] + (1.0-sd.pars[:ρ]) * qpv
	# 	else
	# 		jζp = 1
	# 		Bpv = (1.0 - sd.pars[:ℏ]) * Bpv
	# 		μpv = μprime[jζp]
	# 		σpv = σprime[jζp]
	# 		qpv = qprime[jζp]
	# 		R = (1-sd.pars[:ℏ])*(1-sd.pars[:ρ]) * qpv
	# 	end
	# end
	ζpv = sd.gr[:ζ][jζp]

	ζt == 1.0 && ζpv == 0.0 ? new_def = 1 : new_def = 0

	# λpd = Q' * λt

	# M, V = unmake_logN(μpv, σpv)
	# println(t)
	# println(spr)

	# Fill the path for next period
	if t < length(jz_series)
		jz_series[t+1] = jzp
		fill_path!(p,t+1, Dict(:B => Bpv, :μ => μpv, :σ => σpv, :ζ => ζpv, :ξ => ξpv, :z => zpv, :ψ => prop_domestic, :A => a, :Bh => b, :Bf => Bf, :Wr => Wr, :Wd => Wd, :qg => qpv, :spread=>spr))
		fill_path!(p, t+1, Dict(:Wr10 => Wr_vec[1], :Wr25 => Wr_vec[2], :Wr50 => Wr_vec[3], :Wr75 => Wr_vec[4], :Wr90 => Wr_vec[5], :Wd10 => Wd_vec[1], :Wd25 => Wd_vec[2], :Wd50 => Wd_vec[3], :Wd75 => Wd_vec[4], :Wd90 => Wd_vec[5]))
	end

	return λpd, new_def
end

function iter_simul!(sd::SOEdef, p::Path, t, jz_series, itp_ϕc, itp_ϕs, itp_ϕθ, itp_vf, itp_C, itp_B′, itp_G, itp_pN, itp_qᵍ, itp_repay, itp_W, λt, Qϵ, discr, verbose::Bool=false; B2 = -Inf)

	ϕa, ϕb, Bpv, _, quantiles_ω, sav_a_ω, sav_b_ω, prob_ϵω, jdef = iter_simul_t!(sd, p, t, jz_series, itp_ϕc, itp_ϕs, itp_ϕθ, itp_C, itp_B′, itp_G, itp_pN, itp_repay, λt, discr, verbose, B2 = B2)

	λpd, new_def = iter_simul_tp!(sd, p, t, jz_series, λt, ϕa, ϕb, Bpv, quantiles_ω, sav_a_ω, sav_b_ω, prob_ϵω, jdef, itp_repay, itp_qᵍ, itp_vf, Qϵ)

	return λpd, new_def
end

function simul(sd::SOEdef, jk=1, simul_length::Int64=1, burn_in::Int64=1; ϕ=sd.ϕ, verbose::Bool=false)
	gr = sd.gr
	Random.seed!(jk)

	# Setup
	T = burn_in + simul_length
	p = Path(T = T)

	jz = Int(ceil(length(sd.gr[:z])/2))

	B0, μ0, σ0, ξ0, ζ0, z0 = mean(sd.gr[:b]), mean(sd.gr[:μ]), mean(sd.gr[:σ]), sd.gr[:ξ][1], sd.gr[:ζ][2], sd.gr[:z][jz]
	fill_path!(p,1, Dict(:B => B0, :μ => μ0, :σ => σ0, :w=>1.0, :ξ => ξ0, :ζ => ζ0, :z => z0))

	itp_ϕc = make_itp(sd, ϕ[:c]; agg=false);
	itp_ϕs = make_itp(sd, ϕ[:s]; agg=false);
	itp_ϕθ = make_itp(sd, ϕ[:θ]; agg=false);
	itp_vf = make_itp(sd, sd.v[:v]; agg=false);

	itp_C  = make_itp(sd, sd.eq[:C]; agg=true);
	itp_B′ = make_itp(sd, sd.eq[:issuance]; agg=true);
	itp_G  = make_itp(sd, sd.eq[:spending]; agg=true);
	itp_pN = make_itp(sd, sd.eq[:pN]; agg=true);
	itp_qᵍ = make_itp(sd, sd.eq[:qᵍ]; agg=true);
	itp_W  = make_itp(sd, sd.eq[:welfare]; agg=true);


	rep_mat = reshape_long_shocks(sd, sd.gov[:repay]);
	knots = (gr[:b], gr[:μ], gr[:σ], gr[:ξ], gr[:ζ], gr[:z], gr[:ξ], gr[:z]);
	itp_repay = extrapolate(interpolate(knots, rep_mat, Gridded(Linear())), Interpolations.Flat());

	jz_series = Vector{Int64}(undef, T)
	jz_series[1] = jz

	# Initialize objects for iterating the distribution
	λ = initial_dist(sd, μ0, σ0);
	Qϵ = kron(sd.prob[:ϵ], ones(N(sd,:ωf),1))

	# Simulate
	t0 = time()
	Ndefs = 0
	discr = Dict(:C => 0.0, :pN => 0.0, :μ => 0.0, :σ => 0.0)
	for t in 1:T
		λ, new_def = iter_simul!(sd, p, t, jz_series, itp_ϕc, itp_ϕs, itp_ϕθ, itp_vf, itp_C, itp_B′, itp_G, itp_pN, itp_qᵍ, itp_repay, itp_W, λ, Qϵ, discr, verbose)
		Ndefs += new_def
	end

	# verbose && print_save("\nTime in simulation: $(time_print(time()-t0))")

	jz_series = jz_series[burn_in+1:end]

	return trim_path(p, burn_in), jz_series, Ndefs, discr, λ
end

function iter_simul_switch!(sd1::SOEdef, sd2::SOEdef, p::Path, t, jz_series, itp_ϕc, itp_ϕs, itp_ϕθ, itp_vf, itp_C, itp_B′, itp_G, itp_pN, itp_qᵍ, itp_repay1, itp_repay2, itp_W, λt, Qϵ, discr, verbose::Bool=false; B2=-Inf)

	ϕa, ϕb, Bpv, _, quantiles_ω, sav_a_ω, sav_b_ω, prob_ϵω, jdef = iter_simul_t!(sd1, p, t, jz_series, itp_ϕc, itp_ϕs, itp_ϕθ, itp_C, itp_B′, itp_G, itp_pN, itp_repay1, λt, discr, verbose, B2 = B2)

	λpd, new_def = iter_simul_tp!(sd2, p, t, jz_series, λt, ϕa, ϕb, Bpv, quantiles_ω, sav_a_ω, sav_b_ω, prob_ϵω, jdef, itp_repay2, itp_qᵍ, itp_vf, Qϵ)

	return λpd, new_def
end

function simul_switch!(p::Path, sd1::SOEdef, sd2::SOEdef, jk, length1, length2, length3, B0, μ0, σ0, ξ0, ζ0, z0, λ0, itp_ϕc1, itp_ϕs1, itp_ϕθ1, itp_vf1, itp_C1, itp_B′1, itp_G1, itp_pN1, itp_qᵍ1, itp_W1, itp_repay1, itp_ϕc2, itp_ϕs2, itp_ϕθ2, itp_vf2, itp_C2, itp_B′2, itp_G2, itp_pN2, itp_qᵍ2, itp_W2, itp_repay2, verbose = false; Bvec = [-Inf for jt in 1:length1+length2+length3])
	
	Random.seed!(100+jk)

	# Setup
	T = length1 + length2 + length3

	jz = findfirst(sd1.gr[:z] .== z0)

	fill_path!(p,1, Dict(:B => B0, :μ => μ0, :σ => σ0, :w=>1.0, :ξ => ξ0, :ζ => ζ0, :z => z0))

	jz_series = Vector{Int64}(undef, T)
	jz_series[1] = jz

	# Initialize objects for iterating the distribution
	λ = λ0
	Qϵ = kron(sd1.prob[:ϵ], ones(N(sd1,:ωf),1));

	# Simulate
	t0 = time()
	Ndefs = 0
	discr = Dict(:C => 0.0, :pN => 0.0, :μ => 0.0, :σ => 0.0)
	for t in 1:T
		if t < length1 # simulate in 1
			λ, new_def = iter_simul!(sd1, p, t, jz_series, itp_ϕc1, itp_ϕs1, itp_ϕθ1, itp_vf1, itp_C1, itp_B′1, itp_G1, itp_pN1, itp_qᵍ1, itp_repay1, itp_W1, λ, Qϵ, discr, verbose, B2=Bvec[t])
		elseif t == length1 # switch from 1 to 2
			λ, new_def = iter_simul_switch!(sd1, sd2, p, t, jz_series, itp_ϕc1, itp_ϕs1, itp_ϕθ1, itp_vf2, itp_C1, itp_B′1, itp_G1, itp_pN1, itp_qᵍ2, itp_repay1, itp_repay2, itp_W1, λ, Qϵ, discr, verbose, B2=Bvec[t])
		elseif t < length1 + length2 # simulate in 2
			λ, new_def = iter_simul!(sd2, p, t, jz_series, itp_ϕc2, itp_ϕs2, itp_ϕθ2, itp_vf2, itp_C2, itp_B′2, itp_G2, itp_pN2, itp_qᵍ2, itp_repay2, itp_W2, λ, Qϵ, discr, verbose, B2=Bvec[t])
		elseif t == length1 + length2 # switch back to 1
			λ, new_def = iter_simul_switch!(sd2, sd1, p, t, jz_series, itp_ϕc2, itp_ϕs2, itp_ϕθ2, itp_vf1, itp_C2, itp_B′2, itp_G2, itp_pN2, itp_qᵍ1, itp_repay1, itp_repay2, itp_W2, λ, Qϵ, discr, verbose)
		else # simulate in 1
			λ, new_def = iter_simul!(sd1, p, t, jz_series, itp_ϕc1, itp_ϕs1, itp_ϕθ1, itp_vf1, itp_C1, itp_B′1, itp_G1, itp_pN1, itp_qᵍ1, itp_repay1, itp_W1, λ, Qϵ, discr, verbose)
		end
		Ndefs += new_def
	end

	# verbose && print_save("\nTime in simulation: $(time_print(time()-t0))")

	# jz_series = jz_series[burn_in+1:end]

	return trim_path(p, 1), jz_series, Ndefs, discr
end

function IRF_default(sd::SOEdef, sd_nodef::SOEdef, length1, length2, length3; same_pol=false, burn_in = 4*100, B0 = mean(sd.gr[:b]), K, cond_Y = -Inf)
	pv  = Vector{Path}(undef, K)
	pv2 = Vector{Path}(undef, K)
	pv3 = Vector{Path}(undef, K)

	print_save("Starting $K simulations at $(Dates.format(now(), "HH:MM:SS")).\n")

	t0 = time()
	Threads.@threads for jk in 1:K
		p, _, _, _ = simul_switch(sd_nodef, sd, jk, length1, length2, length3, B0=B0, T = 1+4*100)
		Bpv = series(p, :B)
		p2, _, _, _ = simul_switch(sd_nodef, sd_nodef, jk, length1, length2, length3, B0=B0, T=1+4*100)
		p3, _, _, _ = simul_switch(sd_nodef, sd_nodef, jk, length1, length2, length3, B0=B0, T=1+4*100, Bvec = Bpv)
		
		pv[jk]  = p
		pv2[jk] = p2
		pv3[jk] = p3
	end

	pv2 = [pv2[jk] for jk in eachindex(pv2) if minimum(series(pv[jk], :ζ)) == 1]
	pv3 = [pv3[jk] for jk in eachindex(pv3) if minimum(series(pv[jk], :ζ)) == 1]
	pv = [pv[jk] for jk in eachindex(pv) if minimum(series(pv[jk], :ζ)) == 1]

	print_save("Time in simulation: $(time_print(time()-t0)). Left with $(length(pv)) simulations.\n")

	return pv, length1, length1+length2, pv2, pv3
end


function simul_switch(sd1::SOEdef, sd2::SOEdef, jk, length1, length2, length3; B0=mean(sd1.gr[:b]), μ0=mean(sd1.gr[:μ]), σ0=mean(sd1.gr[:σ]), ξ0=sd1.gr[:ξ][1], ζ0=sd1.gr[:ζ][2], z0=mean(sd1.gr[:z]), T, Bvec=[-Inf for jt in 1:length1+length2+length3])

	_, _, _, _, λ = simul(sd1, jk, T, 1)

	simul_switch(sd1, sd2, jk, length1, length2, length3, λ, B0, μ0, σ0, ξ0, ζ0, z0, Bvec=Bvec)
end


function simul_switch(sd1::SOEdef, sd2::SOEdef, jk, length1, length2, length3, λ0, B0=mean(sd1.gr[:b]), μ0=mean(sd1.gr[:μ]), σ0=mean(sd1.gr[:σ]), ξ0=sd1.gr[:ξ][1], ζ0=sd1.gr[:ζ][2], z0=mean(sd1.gr[:z]); Bvec=[-Inf for jt in 1:length1+length2+length3])

	length1 += 1

	T = length1 + length2 + length3
	p = Path(T = T)

	itp_ϕc1 = make_itp(sd1, sd1.ϕ[:c]; agg=false);
	itp_ϕs1 = make_itp(sd1, sd1.ϕ[:s]; agg=false);
	itp_ϕθ1 = make_itp(sd1, sd1.ϕ[:θ]; agg=false);
	itp_vf1 = make_itp(sd1, sd1.v[:v]; agg=false);

	itp_C1  = make_itp(sd1, sd1.eq[:C]; agg=true);
	itp_B′1 = make_itp(sd1, sd1.eq[:issuance]; agg=true);
	itp_G1  = make_itp(sd1, sd1.eq[:spending]; agg=true);
	itp_pN1 = make_itp(sd1, sd1.eq[:pN]; agg=true);
	itp_qᵍ1 = make_itp(sd1, sd1.eq[:qᵍ]; agg=true);
	itp_W1  = make_itp(sd1, sd1.eq[:welfare]; agg=true);

	rep_mat1 = reshape_long_shocks(sd1, sd1.gov[:repay]);
	knots = (sd1.gr[:b], sd1.gr[:μ], sd1.gr[:σ], sd1.gr[:ξ], sd1.gr[:ζ], sd1.gr[:z], sd1.gr[:ξ], sd1.gr[:z]);
	itp_repay1 = extrapolate(interpolate(knots, rep_mat1, Gridded(Linear())), Interpolations.Flat());

	itp_ϕc2 = make_itp(sd2, sd2.ϕ[:c]; agg=false);
	itp_ϕs2 = make_itp(sd2, sd2.ϕ[:s]; agg=false);
	itp_ϕθ2 = make_itp(sd2, sd2.ϕ[:θ]; agg=false);
	itp_vf2 = make_itp(sd2, sd2.v[:v]; agg=false);

	itp_C2  = make_itp(sd2, sd2.eq[:C]; agg=true);
	itp_B′2 = make_itp(sd2, sd2.eq[:issuance]; agg=true);
	itp_G2  = make_itp(sd2, sd2.eq[:spending]; agg=true);
	itp_pN2 = make_itp(sd2, sd2.eq[:pN]; agg=true);
	itp_qᵍ2 = make_itp(sd2, sd2.eq[:qᵍ]; agg=true);
	itp_W2  = make_itp(sd2, sd2.eq[:welfare]; agg=true);

	rep_mat2 = reshape_long_shocks(sd2, sd2.gov[:repay]);
	knots = (sd2.gr[:b], sd2.gr[:μ], sd2.gr[:σ], sd2.gr[:ξ], sd2.gr[:ζ], sd2.gr[:z], sd2.gr[:ξ], sd2.gr[:z]);
	itp_repay2 = extrapolate(interpolate(knots, rep_mat2, Gridded(Linear())), Interpolations.Flat());

	simul_switch!(p, sd1, sd2, jk, length1, length2, length3, B0, μ0, σ0, ξ0, ζ0, z0, λ0, itp_ϕc1, itp_ϕs1, itp_ϕθ1, itp_vf1, itp_C1, itp_B′1, itp_G1, itp_pN1, itp_qᵍ1, itp_W1, itp_repay1, itp_ϕc2, itp_ϕs2, itp_ϕθ2, itp_vf2, itp_C2, itp_B′2, itp_G2, itp_pN2, itp_qᵍ2, itp_W2, itp_repay2, Bvec = Bvec)
end

function parsimul(sd::SOEdef; ϕ=sd.ϕ, simul_length::Int64=1, burn_in::Int64=1)
	K = Threads.nthreads()
	simul_length = ceil(Int, simul_length/K)

	pv = Vector{Path}(undef, K)
	Ndefs = Vector{Int64}(undef, K)

	print_save("Starting $K simulations of $(floor(Int, simul_length/4)) years at $(Dates.format(now(), "HH:MM")).\n")

	discr_tot = Dict(:C => 0.0, :pN => 0.0)

	t0 = time()
	Threads.@threads for jk in 1:K
		pp, jzs, N, discr, _ = simul(sd, jk, simul_length, burn_in; ϕ=ϕ, verbose=(jk==1))
		pv[jk] = pp
		Ndefs[jk] = N
		for (key, val) in discr_tot
			discr_tot[key] += discr[key] / K
		end
	end
	print_save("Time in simulation: $(time_print(time()-t0))\n")

	for (key, val) in discr_tot
		print_save("Average discrepancy in $(key): $(@sprintf("%0.3g", 100*discr_tot[key]))%\n")
	end

	N = sum(Ndefs)
	return pv, N
end

function get_AR1(y::Vector)
	y_lag = y[1:end-1]
	y = y[2:end]

	data = DataFrame(yt = y, ylag = y_lag)
	OLS = lm(@formula(yt ~ ylag), data)

	ρ = coef(OLS)[2]

	ϵ = y - predict(OLS)
	σ2 = sum(ϵ.^2)/(length(ϵ))
	σϵ = sqrt(σ2)

	σy = sqrt(σ2 / (1-ρ^2))
	return ρ, σϵ
end

function get_AR1(y::Vector, ζ::Vector)

	y_lag   = y[1:end-1]
	new_def = [ζ[tt]-ζ[tt-1] == -1 for tt in 2:length(ζ)]
	y       = y[2:end]

	data = DataFrame(yt = y, ylag = y_lag, X = new_def)
		
	if maximum(new_def) == 0
		OLS = lm(@formula(yt ~ ylag), data)
	else
		OLS = lm(@formula(yt ~ ylag + X), data)
	end

	ρ = coef(OLS)[2]

	ϵ = y - predict(OLS)
	σ2 = sum(ϵ.^2)/(length(ϵ))
	σϵ = sqrt(σ2)

	σy = sqrt(σ2 / (1-ρ^2))
	return ρ, σϵ
end

get_MV(y::Vector) = mean(y), var(y)^0.5

function simul_stats(pv::Vector{T}) where T <: AbstractPath
	K = length(pv)
	v_m = simul_stats(first(pv), verbose=true)
	
	for jp in 2:K
		v_new = simul_stats(pv[jp])

		v_m = (v_m * (jp-1) + v_new) / jp # Recursive average
	end
	return v_m
end

get_spr(q,κ) = 10000 * ( (1+(κ * (1/q - 1)))^4 - 1 )

function simul_stats(path::Path; nodef::Bool=false, ζ_vec::Vector=[], verbose::Bool=false, κ = 0.05985340654896883)
	T = periods(path)
	
	if ζ_vec == []
		ζ_vec = series(path,:ζ)
	end

	conditional = (ζ_vec .> -Inf)
	if nodef
		conditional = (ζ_vec .== 1)
	end

	B_vec = series(path,:B)
	μ_vec = series(path,:μ)
	σ_vec = series(path,:σ)
	w_vec = series(path,:w)
	z_vec = exp.(series(path,:z))
	Y_vec = series(path,:Y)
	C_vec = series(path,:C)
	G_vec = series(path,:G)
	π_vec = series(path,:π)
	ψ_vec = series(path,:ψ)
	mean_vec = series(path, :mean)
	u_vec = 100.0 * (1.0 .- series(path, :L))
	spr_vec = get_spr.(series(path, :qg), κ)
	if haskey(path.data, :spread)
		sprvec = series(path, :spread)
		if norm(spr_vec .- sprvec).>1e-8
			print("WARNING: Difference in calculation of spreads\n")
		end
	end
	Gini_vec = 100*series(path,:Gini)

	# verbose && print("T = $T\n")

	# m_vec, sd_vec = unmake_logN(μ_vec, σ_vec)
	# mean_wealth = 100*mean(m_vec./(4*Y_vec))
	mean_wealth = 25 * mean( mean_vec ./ Y_vec )

	ρy, σy = get_AR1(log.(Y_vec.+1e-8), ζ_vec)
	ρc, σc = get_AR1(log.(C_vec.+1e-8), ζ_vec)
	if var(spr_vec) > 1e-6
		ρs, σs = get_AR1(spr_vec, ζ_vec)
	else
		ρs, σs = 1.0, 0.0
	end
	m_unemp, sd_unemp = get_MV(u_vec)
	m_debt, sd_debt = get_MV(25*B_vec./Y_vec)
	m_gspend, sd_gspend = get_MV(100*G_vec./Y_vec)
	ψ_mean = 100 * median(ψ_vec)
	spr_mean = mean(spr_vec)
	Gini_mean = mean(Gini_vec)

	v_m = [ρy, 100*σy, ρc, 100*σc, ρs, σs, m_debt, sd_debt, m_unemp, sd_unemp, ψ_mean, mean_wealth, Gini_mean]

	return v_m
end

function find_crises(pp::Path, thres::Number, sym::Symbol=:π, k::Int64=7, k_back=2k, thres_back::Number=Inf)
	sym_vec = series(pp,sym)
	ζ_vec = series(pp,:ζ)

	T = periods(pp)

	# println(findall(π_vec .>= πthres))

	Nc = 0
	jt = k_back
	tvec = Int64[]
	while jt < T - k
		# println(jt)
		jt += 1
		if minimum(sym_vec[jt-k:jt+k]) >= thres && sym_vec[jt-k_back] <= thres_back && minimum(ζ_vec[jt-k_back:jt+k]) == 1 # Not default
			Nc += 1
			push!(tvec, jt)
			jt += k
		end
	end

	return Nc, tvec
end

function find_defaults(pp::Path, k::Int64=7)
	ζ_vec = series(pp,:ζ)

	T = periods(pp)

	Nc = 0
	jt = k
	tvec = Int64[]
	while jt < T
		jt += 1
		if ζ_vec[jt-1] == 1 && ζ_vec[jt] == 0 # New default at t
			Nc += 1
			push!(tvec, jt)
			jt += 1
		end
	end
	return Nc, tvec
end

function get_crises(pv::Vector{T}, thres::Number, sym::Symbol=:π, k::Int64=7, k_back=15, thres_back::Number=Inf; type="highspreads") where T <: AbstractPath
	Nc = 0
	tmat = Vector{Vector{Int64}}(undef,0)
	for (jp, pp) in enumerate(pv)
		if type == "highspreads"
			N_new, tvec = find_crises(pp, thres, sym, k, k_back, thres_back)
		elseif type == "default"
			N_new, tvec = find_defaults(pp)
		else
			throw(error("Must choose 'highspreads' or 'default'"))
		end
		Nc += N_new

		tmat = push!(tmat, tvec)
	end

	return Nc, tmat
end

get_crises(pp::Path, thres::Number, k::Int64=7) = get_crises([pp], πthres, k)

series_crises(pp::Path, tvv::Vector{Vector{Int64}}, key::Symbol, k::Int64=9, k_back=k, relative=false) = series_crises([pp], tvv, key, k, k_back, relative=relative)
function series_crises(pv::Vector{T}, tvv::Vector{Vector{Int64}}, key::Symbol, k::Int64=9, k_back = k; relative=false) where T <: AbstractPath
	
	Nc = sum([length(tv) for tv in tvv])
	# symmetric ? k_back = k : k_back = 2k

	ymat = Matrix{Float64}(undef, (k+k_back)+1, Nc)

	jc = 0
	for (jv,tv) in enumerate(tvv)
		if length(tv) > 0
			Y = series(pv[jv], key)
			for (jt, tt) in enumerate(tv)
				if tt-k_back > 0
					jc += 1
					ymat[:, jc] = Y[tt-k_back:tt+k]
					if relative
						ymat[:,jc] = (ymat[:,jc] .- Y[tt-k_back]) / Y[tt-k_back]
					end
				end
			end
		end
	end

	y_up = [quantile(ymat[jt,:], 0.75) for jt in 1:size(ymat,1)]
	y_me = [quantile(ymat[jt,:], 0.50) for jt in 1:size(ymat,1)]
	y_lo = [quantile(ymat[jt,:], 0.25) for jt in 1:size(ymat,1)]
	y_av = mean(ymat, dims=2)[:]

	return ymat, y_up, y_me, y_lo, y_av
end

function MIT_shock(sd::SOEdef, B0 = mean(sd.gr[:b]), ϵb = 0.05; K=100, T=4*10, burn_in = 4*100, verbose=false)
	Random.seed!(1)
	gr = sd.gr

	μ0, σ0 = [mean(sd.gr[key]) for key in [:μ, :σ]]
	ζ0 = gr[:ζ][2] # Start from no default

	jξ = ceil(Int, N(sd,:ξ)/2)
	ξ0 = sd.gr[:ξ][jξ]
	jz = ceil(Int, N(sd,:z)/2)
	z0 = sd.gr[:z][jz]

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

	jz_series = Vector{Int64}(undef, burn_in)
	jz_series[1] = jz

	# Initialize objects for iterating the distribution
	λ = initial_dist(sd, μ0, σ0)
	Qϵ = kron(sd.prob[:ϵ], ones(N(sd,:ωf),1))

	p = Path(T=burn_in)

	fill_path!(p,1, Dict(:B => B0, :μ => μ0, :σ => σ0, :w=>1.0, :ξ => ξ0, :ζ => ζ0, :z => z0))

	# Simulate for 'burn_in' periods to initialize the distirbution
	discr = Dict(:C => 0.0, :pN => 0.0)
	for t in 1:burn_in
		λ, new_def = iter_simul!(sd, p, t, jz_series, itp_ϕc, itp_ϕs, itp_ϕθ, itp_vf, itp_C, itp_B′, itp_G, itp_pN, itp_qᵍ, itp_repay, itp_W, λ, Qϵ, discr, verbose)
	end

	μ0, σ0 = [p.data[key][end] for key in [:μ, :σ]]
	λhigh = copy(λ)

	pv = Vector{Path}(undef, K)
	pv_high = Vector{Path}(undef, K)
	Threads.@threads for jk in 1:K
		Random.seed!(jk)
		pv[jk] = Path(T=T)
		fill_path!(pv[jk],1, Dict(:B => B0, :μ => μ0, :σ => σ0, :w=>1.0, :ξ => ξ0, :ζ => ζ0, :z => z0))
		jz_series = Vector{Int64}(undef, T)
		jz_series[1] = jz
		discr = Dict(:C => 0.0, :pN => 0.0)
		for t in 1:T
			λ, new_def = iter_simul!(sd, pv[jk], t, jz_series, itp_ϕc, itp_ϕs, itp_ϕθ, itp_vf, itp_C, itp_B′, itp_G, itp_pN, itp_qᵍ, itp_repay, itp_W, λ, Qϵ, discr, verbose)
		end
		pv_high[jk] = Path(T=T)
		Bhigh = (1+ϵb)*B0
		Random.seed!(jk)
		fill_path!(pv_high[jk],1, Dict(:B => Bhigh, :μ => μ0, :σ => σ0, :w=>1.0, :ξ => ξ0, :ζ => ζ0, :z => z0))
		discr = Dict(:C => 0.0, :pN => 0.0)
		for t in 1:T
			λhigh, new_def = iter_simul!(sd, pv_high[jk], t, jz_series, itp_ϕc, itp_ϕs, itp_ϕθ, itp_vf, itp_C, itp_B′, itp_G, itp_pN, itp_qᵍ, itp_repay, itp_W, λhigh, Qϵ, discr, verbose)
		end
	end

	return pv, pv_high
end


function simul_stats_nodef(pv::Vector{T}) where T <: AbstractPath

	pvec = Vector{AbstractPath}(undef, 0 )

	for pp in pv
		def_state = (series(pp,:ζ) .== 0)

		defs = findall([def_state[tt]-def_state[tt-1].==1  for tt in 2:length(def_state)])
		reen = findall([def_state[tt]-def_state[tt-1].==-1 for tt in 2:length(def_state)])

		if defs[1] > reen[1]
			print("WARNING: started in default")
		else
			reen = [1; reen]
		end

		if defs[end] == periods(pp)
		elseif length(reen) > length(defs)
			defs = [defs; periods(pp)]
		end

		for (jd, tt) in enumerate(defs)
			len = defs[jd] - reen[jd] + 1
			
			if len > 2
				path_between = Path{len}(Dict(key => val[reen[jd]:defs[jd]] for (key, val) in pp.data))

				push!(pvec, path_between)
			end
		end
	end

	println("Left with $(length(pvec)) paths")

	T_vec = [periods(pp) for pp in pvec]
	V_vec = [simul_stats(pp) for pp in pvec]

	mean_V = sum(V_vec .* T_vec / sum(T_vec))
end

rel_var(p_bench::Path, p_nodef::Path, key::Symbol) = rel_var([p_bench], [p_nodef], key)
function rel_var(p_bench::Vector{T}, p_nodef::Vector{T}, key::Symbol) where T <: AbstractPath
	rel(x,y) = (y-x) / x
	return rel(mean([mean(series(p, key)) for p in p_bench]), mean([mean(series(p, key)) for p in p_nodef]))
end
