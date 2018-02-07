utility(h::Hank, c::Float64) = ifelse(c > 1e-10, c^(1.0 - h.γ) / (1.0 - h.γ), -1e10)

function uprime(h::Hank, c_vec)
	u = zeros(size(c_vec))
	for (jc, cv) in enumerate(c_vec)
		if cv > 0
			u[jc] = cv.^(-h.γ)
		else
			u[jc] = 1e10
		end
	end
	if length(c_vec) == 1
		return u[1]
	else
		return u
	end
end

function uprime_inv(h::Hank, c_vec::Vector)
	u = zeros(size(c_vec))
	for (jc, cv) in enumerate(c_vec)
		if cv > 0
			u[jc] = cv.^(-1./h.γ)
		else
			u[jc] = 1e-10
		end
	end
	if length(c_vec) == 1
		return u[1]
	else
		return u
	end
end

function get_abc(RHS::Float64, ωmin::Float64, qʰ::Float64, qᵍ::Float64, pC::Float64, sp::Float64, θp::Float64)
	""" Recovers private and public debt purchases and consumption from savings decisions """
	ap = ωmin + θp * (sp - qʰ*ωmin) / qʰ
	bp = (1.0-θp)  * (sp - qʰ*ωmin) / qᵍ
	C  = (RHS - sp) / pC

	return ap, bp, C
end

function value(h::Hank, sp::Float64, θp::Float64, itp_vf, jϵ, jz, thres, RHS, qʰ, qᵍ, pC, bpv, μ′, σ′, wpv, itp_R, jdefault)

	ap, bp, C = get_abc(RHS, h.ωmin, qʰ, qᵍ, pC, sp, θp)

	""" CHANGE THIS FOR GHH """
	ℓ = 0

	u = utility(h, C - ℓ)

	# Basis matrix for continuation values
	check, Ev, test, ut = 0., 0., 0, 0.

	# itp_vf = itp_obj_vf
	for (jzp, zpv) in enumerate(h.zgrid)

		μpv = μ′[jzp, 1]
		σpv = σ′[jzp, 1]

		ζpv = 1.0
		if jdefault == false && zpv <= thres
			ζpv = 2.0
		end

		ωpv = ap + bp * itp_R[bpv, μpv, σpv, wpv, ζpv, zpv]
		ωpv = min(ωpv, maximum(h.ωgrid))		
		
		for (jϵp, ϵpv) in enumerate(h.ϵgrid)
			prob =  h.Pz[jz, jzp] * h.Pϵ[jϵ, jϵp]
			if jdefault
				ζ_reent, ζ_cont = 1.0, 3.0
				Ev += (itp_vf[ωpv, ϵpv, bpv, μ′[jzp, 1], σ′[jzp, 1], wpv, ζ_reent, zpv])^(1.0-h.γ) * prob * h.θ

				Ev += (itp_vf[ωpv, ϵpv, bpv, μ′[jzp, 2], σ′[jzp, 2], wpv, ζ_cont,  zpv])^(1.0-h.γ) * prob * (1.0 - h.θ)
			else
				Ev += (itp_vf[ωpv, ϵpv, bpv, μpv, σpv, wpv, ζpv, zpv])^(1.0-h.γ) * prob
			end
			check += prob
		end
		# isnan(Ev)? print_save("\n $test NaNs out of $(length(h.zgrid))"): Void
	end
	# test > 0? print_save("\n $test NaNs out of $(length(h.zgrid))"): Void
	Tv = Ev^(1.0/(1.0-h.γ))


	isapprox(check, 1) || warn("wrong expectation operator")

	# Compute value
	if h.EpsteinZin
		C > 1e-10? ut = C^((h.ψ-1.0)/h.ψ): ut = 1e-10
		vf = (1.0 - h.β) * ut + h.β * Tv^((h.ψ-1.0)/h.ψ)
		vf = vf^((h.ψ)/(h.ψ-1.0))
		return vf
	else
		vf = (1.0 - h.β) * u + h.β * Ev
		return vf
	end
	Void
end

function solve_optvalue(h::Hank, guess::Vector, itp_vf, jϵ, jz, thres, RHS, qʰv, qᵍv, pCv, bpv, μpv, σpv, wpv, itp_R, jdef, ωmax)

	both = false

	ap, bp, cmax, fmax = 0., 0., 0., 0.
	if both
		wrap_value(x::Vector, grad::Vector=similar(x)) = value(h, x[1], x[2], itp_vf, jϵ, jz, thres, RHS, qʰv, qᵍv, pCv, bpv, μpv, σpv, wpv, itp_R, jdef)

		opt = Opt(:LN_COBYLA, length(guess))
		# upper_bounds!(opt, [ωmax, 1e-2])
		upper_bounds!(opt, [ωmax, 1.0])
		lower_bounds!(opt, [qʰv*ωmin, 0.0])
		# xtol_abs!(opt, 1e-16)
		# ftol_abs!(opt, 1e-16)

		max_objective!(opt, wrap_value)
		(fmax, xmax, ret) = NLopt.optimize(opt, guess)
		# ret == :SUCCESS || println(ret)

		sp = xmax[1]
		θp = xmax[2]
		ap, bp, cmax = get_abc(RHS, h.ωmin, qʰv, qᵍv, pCv, sp, θp)
	else
		curr_min = 1e10
		for θp in 0:0.25:1
			res = Optim.optimize(
				sp -> -value(h, sp, θp, itp_vf, jϵ, jz, thres, RHS, qʰv, qᵍv, pCv, bpv, μpv, σpv, wpv, itp_R, jdef),
					qʰv*h.ωmin, ωmax, GoldenSection()
				)
			if res.minimum < curr_min
				sp = res.minimizer
				ap, bp, cmax = get_abc(RHS, h.ωmin, qʰv, qᵍv, pCv, sp, θp)
				fmax = -res.minimum
				curr_min = res.minimum
			end
		end
	end
		
	return ap, bp, cmax, fmax
end


function opt_value(h::Hank, qʰ_mat, qᵍ_mat, wL_mat, T_mat, pC_mat, itp_R, itp_vf; resolve::Bool = true)

	vf = Array{Float64}(size(h.vf))
	ϕa = Array{Float64}(size(h.ϕa))
	ϕb = Array{Float64}(size(h.ϕb))
	ϕc = Array{Float64}(size(h.ϕc))
	for js in 1:size(h.Jgrid,1)
		jb = h.Jgrid[js, 1]
		jμ = h.Jgrid[js, 2]
		jσ = h.Jgrid[js, 3]
		jw = h.Jgrid[js, 4]
		jζ = h.Jgrid[js, 5]
		jz = h.Jgrid[js, 6]

		qʰv = qʰ_mat[jb, jμ, jσ, jw, jζ, jz]
		qᵍv = qᵍ_mat[jb, jμ, jσ, jw, jζ, jz]
		wv  = wL_mat[jb, jμ, jσ, jw, jζ, jz]
		Tv  = T_mat[jb, jμ, jσ, jw, jζ, jz]
		pCv = pC_mat[jb, jμ, jσ, jw, jζ, jz]

		bpv = h.issuance[js]
		μpv = h.μ′[js,:,:]
		σpv = h.σ′[js,:,:]
		wpv = h.w′[js]
		thres = h.def_thres[js]

		minimum(μpv) < minimum(h.μgrid) || maximum(μpv) > maximum(h.μgrid)? print_save("μ out of bounds"): Void
		minimum(σpv) < minimum(h.σgrid) || maximum(σpv) > maximum(h.σgrid)? print_save("σ out of bounds"): Void
		bpv < minimum(h.bgrid) || bpv > maximum(h.bgrid)? print_save("b out of bounds"): Void
		wpv < minimum(h.wgrid) || wpv > maximum(h.wgrid)? print_save("w out of bounds"): Void


		jdef = (h.ζgrid[jζ] != 1.0)

		for (jϵ, ϵv) in enumerate(h.ϵgrid), (jω, ωv) in enumerate(h.ωgrid)

			RHS = ωv + wv * ϵv - Tv

			ap, bp, fmax = 0., 0., 0.
			ag, bg = h.ϕa[jω, jϵ, jb, jμ, jσ, jw, jζ, jz], h.ϕb[jω, jϵ, jb, jμ, jσ, jw, jζ, jz]

			ωg = qʰv * ag + qᵍv * bg
			θg = qʰv * (ag - h.ωmin) / (ωg - qʰv*h.ωmin)
			if resolve
				# θg = 1.0

				ωmax = RHS - 1e-10

				if ωmax < qʰv * h.ωmin
					print_save("\nCan't afford positive consumption at $([jb, jμ, jσ, jw, jζ, jz]) with w*Lᵈ=$(round(wv,2)), T=$(round(Tv,2))")
					ap, bp, cmax = h.ωmin, 0., 1e-10
					fmax = value(h, qʰv*ap, 0., itp_vf, jϵ, jz, thres, RHS, qʰv, qᵍv, pCv, bpv, μpv, σpv, wpv, itp_R, jdef)
				else
					if ωg > ωmax
						ωg = max(ωmax - 1e-2, 0)
					end
					guess = [ωg, θg]

					ap, bp, cmax, fmax = solve_optvalue(h, guess, itp_vf, jϵ, jz, thres, RHS, qʰv, qᵍv, pCv, bpv, μpv, σpv, wpv, itp_R, jdef, ωmax)
				end
			else
				ap = h.ϕa[jω, jϵ, jb, jμ, jσ, jw, jζ, jz]
				bp = h.ϕb[jω, jϵ, jb, jμ, jσ, jw, jζ, jz]
				cmax = h.ϕc[jω, jϵ, jb, jμ, jσ, jw, jζ, jz]
				fmax = value(h, ωg, θp, itp_vf, jϵ, jz, thres, RHS, qʰv, qᵍv, pCv, bpv, μpv, σpv, wpv, itp_R, jdef)
			end
			cmax < 0? warn("c = $cmax"): Void
			
			ϕa[jω, jϵ, jb, jμ, jσ, jw, jζ, jz] = ap
			ϕb[jω, jϵ, jb, jμ, jσ, jw, jζ, jz] = bp
			ϕc[jω, jϵ, jb, jμ, jσ, jw, jζ, jz] = cmax
			vf[jω, jϵ, jb, jμ, jσ, jw, jζ, jz] = fmax
		end
	end

	return vf, ϕa, ϕb, ϕc
end

function bellman_iteration!(h::Hank, qʰ_mat, qᵍ_mat, wL_mat, T_mat, R_mat, pC_mat; resolve::Bool=true)
	# Interpolate the value function
	all_knots = (h.ωgrid, h.ϵgrid, h.bgrid, h.μgrid, h.σgrid, h.wgrid, 1:3, h.zgrid)

	sum(isnan(h.vf)) > 0? print_save("\n $(sum(isnan(h.vf))) NaNs detected in h.vf"): Void
	sum(h.vf .<= 0) > 0? print_save("\n $(sum(h.vf .<= 0)) negative entries detected in h.vf"): Void

	itp_vf = interpolate(all_knots, h.vf, (Gridded(Linear()), Gridded(Linear()), Gridded(Linear()), Gridded(Linear()), Gridded(Linear()), Gridded(Linear()), NoInterp(), Gridded(Linear())))
	itp_R  = interpolate((h.bgrid, h.μgrid, h.σgrid, h.wgrid, 1:3, h.zgrid), R_mat, (Gridded(Linear()), Gridded(Linear()), Gridded(Linear()), Gridded(Linear()), NoInterp(), Gridded(Linear())))

	# Compute values
	vf, ϕa, ϕb, ϕc = opt_value(h, qʰ_mat, qᵍ_mat, wL_mat, T_mat, pC_mat, itp_R, itp_vf, resolve = resolve)

	sum(isnan.(vf)) > 0? print_save("\n$(sum(isnan.(vf))) found in vf"): Void
	sum(isnan.(ϕa)) > 0? print_save("$(sum(isnan.(ϕa))) found in ϕa"): Void
	sum(isnan.(ϕb)) > 0? print_save("$(sum(isnan.(ϕb))) found in ϕb"): Void
	sum(isnan.(ϕc)) > 0? print_save("$(sum(isnan.(ϕc))) found in ϕc"): Void

	# Store results in the type
	h.ϕa = ϕa
	h.ϕb = ϕb
	h.ϕc = ϕc
	h.vf = vf

	Void
end