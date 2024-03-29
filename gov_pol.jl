function integrate_itp(sd::SOEdef, bv, μv, σv, ξv, ζv, zv, itp_obj)
	pars, gr = sd.pars, sd.gr

	λϵ = ergodic_ϵ(sd)

	ωmin_int, ωmax_int = quantile.(LogNormal(μv, σv), [1e-6; 1-1e-6]) .+ pars[:ωmin]
	ωmax_int = min(ωmax_int, maximum(gr[:ω]))
	# ωmin_int = min(ωmin_int, maximum(gr[:ω]) - 1)
	# ωmin_int = pars[:ωmin]
	# itp_obj = extrapolate(itp_obj, Interpolations.Line())
	W, sum_prob = 0.0, 0.0
	for (jϵ, ϵv) in enumerate(gr[:ϵ])
		f_pdf(ω) = pdf(LogNormal(μv, σv), ω.-pars[:ωmin])
		(val_pdf, err) = hquadrature(f_pdf, ωmin_int, ωmax_int, rtol=1e-10, atol=1e-12, maxevals=0)
		sum_prob += val_pdf * λϵ[jϵ]

	    f(ω) = f_pdf(ω) * itp_obj(ω, ϵv, bv, μv, σv, ξv, ζv, zv)
	    (val, err) = hquadrature(f, ωmin_int, ωmax_int, rtol=1e-10, atol=1e-12, maxevals=0)
	    W += val * λϵ[jϵ]
	end
	W = W / sum_prob

	return W
end

function update_govpol(sd::SOEdef)
	B′ = sd.eq[:issuance]	

	itp_W = make_itp(sd, sd.eq[:welfare]; agg=true)

	# More μ means default less often
	# μ_gov = 0.016 #* 0.0
	μ_gov = sd.pars[:μ_gov]
	σ_gov = 0.0015

	Jgrid = agg_grid(sd);
	rep_prob = zeros(N(sd,:b), N(sd,:μ), N(sd,:σ), N(sd,:ξ), N(sd,:ζ), N(sd,:z), N(sd,:ξ), N(sd,:z))
	for js in axes(Jgrid, 1)
		μ′_arr = sd.LoM[:μ][js,:,:]
		σ′_arr = sd.LoM[:σ][js,:,:]

		jb, jμ, jσ, jξ, jζ, jz = Jgrid[js, :]

		ζv = sd.gr[:ζ][Jgrid[js, 5]]
		
		bpv = B′[js]
		for (jξp, ξpv) in enumerate(sd.gr[:ξ]), (jzp, zpv) in enumerate(sd.gr[:z])
			if ζv == 1 # NO default at t
				jζp = 1 # Default at t+1
				ζpv = sd.gr[:ζ][jζp]
				μpv = μ′_arr[jξp, jzp][jζp]
				σpv = σ′_arr[jξp, jzp][jζp]
				Wd = itp_W((1.0-sd.pars[:ℏ])*bpv, μpv, σpv, ξpv, ζpv, zpv)

				jζp = 2 # No default at t+1
				ζpv = sd.gr[:ζ][jζp]
				μpv = μ′_arr[jξp, jzp][jζp]
				σpv = σ′_arr[jξp, jzp][jζp]
				Wr = itp_W(bpv, μpv, σpv, ξpv, ζpv, zpv)

				rep_prob[jb, jμ, jσ, jξ, jζ, jz, jξp, jzp] = 1.0 - cdf(Normal(μ_gov, σ_gov), Wd-Wr)
			else # If default at t, no decision to make
				rep_prob[jb, jμ, jσ, jξ, jζ, jz, jξp, jzp] = 0.
			end
		end
	end

	rep_new = reshape(rep_prob, length(rep_prob))

	return rep_new
end

function update_W(sd::SOEdef)
	itp_vf = make_itp(sd, sd.v[:v]; agg=false)

	Jgrid = agg_grid(sd)
	W = zeros(size(Jgrid, 1))
	for js in axes(Jgrid, 1)
		bv = sd.gr[:b][Jgrid[js, 1]]
		μv = sd.gr[:μ][Jgrid[js, 2]]
		σv = sd.gr[:σ][Jgrid[js, 3]]
		ξv = sd.gr[:ξ][Jgrid[js, 4]]
		ζv = sd.gr[:ζ][Jgrid[js, 5]]
		zv = sd.gr[:z][Jgrid[js, 6]]
		
		W[js] = integrate_itp(sd, bv, μv, σv, ξv, ζv, zv, itp_vf)
	end

	return W
end


function mpe_iter!(sd::SOEdef; maxiter::Int64=500, tol::Float64=5e-4, nodef::Bool=sd.opt[:nodef], noΔ::Bool=sd.opt[:noΔ], rep_agent::Bool=sd.opt[:rep_agent], run_number::Int64=1, save_copies::Bool=false, nob::Bool=false, verbose::Bool=false)

	init_msg = "\nIterating on the government's policy"
	if nodef
		init_msg *= " in nodef"
	elseif noΔ
		init_msg *= " in noΔ"
		sd.pars[:Δ] = 0
	elseif nob
		init_msg *= " in nob"
	end
	init_msg*= ":"
	
	print_save(init_msg)
	time_init = time()
	t_old = time_init
	iter = 0
	dist = 1+tol

	upd_η = 1.
	upd_ηR = 0.14

	tol_eqm = 0.1
	maxiter_CE = 100

	if nodef
		# Make sure the government never defaults
		sd.gov[:repay] = ones(size(sd.gov[:repay]))
	end
	if nob
		sd.opt[:nob] = true
	end
	if noΔ
		sd.opt[:noΔ] = true
	end
	Jgrid = agg_grid(sd);

	while iter < 2 || (dist > tol && iter < maxiter)
		iter += 1
		print_save("\n\nOuter Iteration $iter (run $(run_number)) with upd_ηR = $(@sprintf("%0.3g",upd_ηR)) at $(Dates.format(now(), "HH:MM"))")
		print_save("\nqᵍ between $(round(minimum(sd.eq[:qᵍ][Jgrid[:,5].==1]),digits=4)) and $(round(maximum(sd.eq[:qᵍ]),digits=4)). risk-free is $(round(mean(sd.eq[:qʰ]),digits=4))")
		print_save(" (spread between $(floor(Int,10000*minimum(sd.eq[:spread]))) bps and $(floor(Int,10000*maximum(sd.eq[:spread][Jgrid[:,5].==1]))) bps)")

		""" RUN COMP_EQM LOOP """
		dist_CE1 = comp_eqm!(sd, verbose = verbose, tol = tol_eqm, maxiter = maxiter_CE)
		dist_CE = min(dist_CE1, tol_eqm)

		verbose && print_save("\nNew pN = [$(@sprintf("%0.3g",minimum(sd.eq[:pN]))), $(@sprintf("%0.3g",maximum(sd.eq[:pN])))]")
		verbose && print_save("\nNew w  = [$(@sprintf("%0.3g",minimum(sd.eq[:wage]))), $(@sprintf("%0.3g",maximum(sd.eq[:wage])))]")
		verbose && print_save("\nNew Ld = [$(@sprintf("%0.3g",minimum(sd.eq[:Ld]))), $(@sprintf("%0.3g",maximum(sd.eq[:Ld])))]")


		""" UPDATES """
		W_new = update_W(sd)
		sd.eq[:welfare] = upd_η * W_new + (1.0-upd_η) * sd.eq[:welfare]

		if isnan.(sum(sd.eq[:welfare]))
			print_save("\nWARNING: ‖welf‖ = NaN")
		end

		old_rep = copy(sd.gov[:repay])

		if nodef || noΔ || nob
			# Keep the same default policy
			dist = 0.0
		else
			# Update the default policy
			new_rep = update_govpol(sd)
			old_norm = sqrt.(sum(old_rep.^2))
			# print_save("\n||rep₀, repₜ|| = $(@sprintf("%0.3g",old_norm))")
			if isapprox(old_norm, 0.0)
				old_norm = 1.0
			end
			new_norm = sqrt.(sum(new_rep.^2))
			# print_save(", $((@sprintf("%0.3g",new_norm)))")
			dist = sqrt.(sum( (new_rep - old_rep).^2 )) / old_norm
			sd.gov[:repay] = upd_ηR * new_rep + (1.0-upd_ηR) * old_rep
		end

		if save_copies && iter % 50 == 0
			save(pwd() * "/../Output/run$(run_number)/SOEdef.jld", "sd", sd)
		end

		tol_eqm = max(exp(0.9*log(1+min(dist_CE1,tol_eqm)))-1, dist, tol)
		upd_ηR = max(upd_ηR * 0.995, 5e-2)
		t_new = time()
		print_save("\nDistance = $(@sprintf("%0.3g",dist)) after $(time_print(t_new-t_old)) and $iter iterations. New tol = $(@sprintf("%0.3g",tol_eqm))")

		dist = max(dist, dist_CE)

		time_old = time()

		if nodef == false && noΔ == false && nob == false
			maxiter_CE = 25
		end
	end
	if save_copies
		save(pwd() * "/../Output/run$(run_number)/SOEdef.jld", "sd", sd)
	end
	if dist <= tol
		print_save("\n\nConverged in $iter iterations. ")
	else
		print_save("\n\nStopping at distance $(@sprintf("%0.3g",dist)). ")
	end


	print_save("\nTotal time: $(time_print(time()-time_init))\n")
end
