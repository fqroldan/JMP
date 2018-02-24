utility(h::Hank, c::Float64) = ifelse(c > 1e-10, c^(1.0 - h.ψ) / (1.0 - h.ψ), -1e10)

function G(h::Hank, v::Float64)
	if h.EpsteinZin
		if h.γ != 1
			return v^(1.0-h.γ)
		else
			return log(v)
		end
	else
		return v
	end
end

function T(h::Hank, Ev::Float64)
	if h.γ != 1
		return Ev^(1.0/(1.0-h.γ))
	else
		return exp(Ev)
	end
end

function uprime(h::Hank, c_vec)
	u = zeros(size(c_vec))
	for (jc, cv) in enumerate(c_vec)
		if cv > 0
			u[jc] = cv.^(-h.ψ)
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
			u[jc] = cv.^(-1./h.ψ)
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

function value(h::Hank, sp::Float64, θp::Float64, itp_obj_vf, jϵ, jz, thres, RHS, qʰ, qᵍ, pC, bpv, μ′, σ′, wpv, itp_qᵍ, jdefault)

	ap, bp, C = get_abc(RHS, h.ωmin, qʰ, qᵍ, pC, sp, θp)

	""" CHANGE THIS FOR GHH """
	ℓ = 0

	u = utility(h, C - ℓ)

	# Basis matrix for continuation values
	check, Ev, test, ut = 0., 0., 0, 0.

	if jdefault
		for (jzp, zpv) in enumerate(h.zgrid)
			for (jϵp, ϵpv) in enumerate(h.ϵgrid)
				prob =  h.Pz[jz, jzp] * h.Pϵ[jϵ, jϵp]

				# Reentry
				ζpv = 1
				μpv, σpv = μ′[jzp, 1], σ′[jzp, 1]
				R = h.κ + (1.0 - h.ρ) * itp_qᵍ[bpv, μpv, σpv, wpv, ζpv, jzp]
				ωpv = ap + bp * R
				ωpv = min(ap + bp * R, h.ωmax)
				if ωpv > h.ωmax
					itp_vf = extrapolate(itp_obj_vf, Linear())
				else
					itp_vf = itp_obj_vf
				end
				Ev += G(h, itp_vf[ωpv, jϵp, bpv, μpv, σpv, wpv, ζpv, jzp]) * prob * h.θ
				
				# Continue in default
				ζpv = 2
				μpv, σpv = μ′[jzp, 2], σ′[jzp, 2]
				R = (1.0 - h.ρ) * itp_qᵍ[bpv, μpv, σpv, wpv, ζpv, jzp]
				ωpv = ap + bp * R
				ωpv = min(ap + bp * R, h.ωmax)
				if ωpv > h.ωmax
					itp_vf = extrapolate(itp_obj_vf, Linear())
				else
					itp_vf = itp_obj_vf
				end
				Ev += G(h, itp_vf[ωpv, jϵp, bpv, μpv, σpv, wpv, ζpv, jzp]) * prob * (1.0 - h.θ)
				check += prob
			end
		end
	else
		for (jzp, zpv) in enumerate(h.zgrid)
			for (jϵp, ϵpv) in enumerate(h.ϵgrid)
				prob =  h.Pz[jz, jzp] * h.Pϵ[jϵ, jϵp]
				check += prob

				μpv = μ′[jzp, 1]
				σpv = σ′[jzp, 1]				
				if zpv > thres
					ζpv = 1
					R = h.κ + (1.0 - h.ρ) * itp_qᵍ[bpv, μpv, σpv, wpv, ζpv, jzp]
					ωpv = ap + bp * R
					ωpv = min(ap + bp * R, h.ωmax)
					if ωpv > h.ωmax
						itp_vf = extrapolate(itp_obj_vf, Linear())
					else
						itp_vf = itp_obj_vf
					end
					Ev += G(h, itp_vf[ωpv, jϵp, bpv, μpv, σpv, wpv, ζpv, jzp]) * prob
				else
					ζpv = 2
					R = (1.0 - h.ℏ) * (1.0 - h.ρ) * itp_qᵍ[(1.0 - h.ℏ)*bpv, μpv, σpv, wpv, ζpv, jzp]
					ωpv = ap + bp * R
					ωpv = min(ap + bp * R, h.ωmax)
					if ωpv > h.ωmax
						itp_vf = extrapolate(itp_obj_vf, Linear())
					else
						itp_vf = itp_obj_vf
					end
					Ev += G(h, itp_vf[ωpv, jϵp, (1.0 - h.ℏ)*bpv, μpv, σpv, wpv, ζpv, jzp]) * prob				
				end
			end
		end
	end
	Tv = T(h, Ev)


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

function solve_optvalue(h::Hank, guess::Vector, itp_vf, jϵ, jz, thres, RHS, qʰv, qᵍv, pCv, bpv, μpv, σpv, wpv, itp_qᵍ, jdef, ωmax)

	both = false

	ap, bp, cmax, fmax = 0., 0., 0., 0.
	if both
		wrap_value(x::Vector, grad::Vector=similar(x)) = value(h, x[1], x[2], itp_vf, jϵ, jz, thres, RHS, qʰv, qᵍv, pCv, bpv, μpv, σpv, wpv, itp_qᵍ, jdef)

		opt = Opt(:LN_BOBYQA, length(guess))
		# upper_bounds!(opt, [ωmax, 1e-2])
		upper_bounds!(opt, [ωmax, 1.0])
		lower_bounds!(opt, [qʰv*h.ωmin, 0.0])
		# xtol_abs!(opt, 1e-16)
		ftol_rel!(opt, 1e-4)

		guess[1] = max(min(guess[1], ωmax), qʰv*h.ωmin)
		guess[2] = max(min(guess[2], 1.0), 0.0)

		max_objective!(opt, wrap_value)
		(fmax, xmax, ret) = NLopt.optimize(opt, guess)
		# ret == :SUCCESS || println(ret)

		sp = xmax[1]
		θp = xmax[2]
		ap, bp, cmax = get_abc(RHS, h.ωmin, qʰv, qᵍv, pCv, sp, θp)
	else
		curr_min = 1e10
		θp_grid = 1-[0 0.0125 0.025 0.05 0.1 0.2 0.5 1]
		for θp in θp_grid
			res = Optim.optimize(
				sp -> -value(h, sp, θp, itp_vf, jϵ, jz, thres, RHS, qʰv, qᵍv, pCv, bpv, μpv, σpv, wpv, itp_qᵍ, jdef),
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


function opt_value(h::Hank, qʰ_mat, qᵍ_mat, wL_mat, T_mat, pC_mat, itp_qᵍ, itp_vf; resolve::Bool = true)

	vf = SharedArray{Float64}(size(h.vf))
	ϕa = SharedArray{Float64}(size(h.ϕa))
	ϕb = SharedArray{Float64}(size(h.ϕb))
	ϕc = SharedArray{Float64}(size(h.ϕc))
	@sync @parallel for js in 1:size(h.Jgrid,1)
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

			RHS = ωv + wv * exp(ϵv) - Tv

			ap, bp, fmax = 0., 0., 0.
			ag, bg = h.ϕa[jω, jϵ, jb, jμ, jσ, jw, jζ, jz], h.ϕb[jω, jϵ, jb, jμ, jσ, jw, jζ, jz]

			ωg = qʰv * ag + qᵍv * bg
			θg = qʰv * (ag - h.ωmin) / (ωg - qʰv*h.ωmin)
			# print_save("a,b,s,θ = $([ag, bg, ωg, θg])")
			if resolve
				# θg = 1.0

				ωmax = RHS - 1e-10

				if ωmax < qʰv * h.ωmin
					print_save("\nCan't afford positive consumption at $([jb, jμ, jσ, jw, jζ, jz]) with w*Lᵈ=$(round(wv,2)), T=$(round(Tv,2))")
					ap, bp, cmax = h.ωmin, 0., 1e-10
					fmax = value(h, qʰv*ap, 0., itp_vf, jϵ, jz, thres, RHS, qʰv, qᵍv, pCv, bpv, μpv, σpv, wpv, itp_qᵍ, jdef)
				else
					if ωg > ωmax
						ωg = max(ωmax - 1e-2, 0)
					end
					isapprox(θg, 1) && θg > 1? θg = 1.0: Void
					guess = [ωg, θg]

					ap, bp, cmax, fmax = solve_optvalue(h, guess, itp_vf, jϵ, jz, thres, RHS, qʰv, qᵍv, pCv, bpv, μpv, σpv, wpv, itp_qᵍ, jdef, ωmax)
				end
			else
				ap = h.ϕa[jω, jϵ, jb, jμ, jσ, jw, jζ, jz]
				bp = h.ϕb[jω, jϵ, jb, jμ, jσ, jw, jζ, jz]
				cmax = h.ϕc[jω, jϵ, jb, jμ, jσ, jw, jζ, jz]
				fmax = value(h, ωg, θp, itp_vf, jϵ, jz, thres, RHS, qʰv, qᵍv, pCv, bpv, μpv, σpv, wpv, itp_qᵍ, jdef)
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

function bellman_iteration!(h::Hank, qʰ_mat, qᵍ_mat, wL_mat, T_mat, pC_mat; resolve::Bool=true)
	# Interpolate the value function
	# all_knots = (h.ωgrid, 1:h.Nϵ, h.bgrid, h.μgrid, h.σgrid, h.wgrid, 1:h.Nζ, 1:h.Nz)
	# agg_knots = (h.bgrid, h.μgrid, h.σgrid, h.wgrid, 1:h.Nζ, 1:h.Nz)

	# sum(isnan.(h.vf)) > 0? print_save("\n $(sum(isnan.(h.vf))) NaNs detected in h.vf"): Void
	# sum(h.vf .<= 0) > 0? print_save("\n $(sum(h.vf .<= 0)) negative entries detected in h.vf"): Void

	# itp_vf = interpolate(all_knots, h.vf, (Gridded(Linear()), NoInterp(), Gridded(Linear()), Gridded(Linear()), Gridded(Linear()), Gridded(Linear()), NoInterp(), NoInterp()))
	# itp_qᵍ  = interpolate(agg_knots, qᵍ_mat, (Gridded(Linear()), Gridded(Linear()), Gridded(Linear()), Gridded(Linear()), NoInterp(), NoInterp()))

	ωrange = linspace(h.ωgrid[1], h.ωgrid[end], h.Nω)
	brange = linspace(h.bgrid[1], h.bgrid[end], h.Nb)
	μrange = linspace(h.μgrid[1], h.μgrid[end], h.Nμ)
	σrange = linspace(h.σgrid[1], h.σgrid[end], h.Nσ)
	wrange = linspace(h.wgrid[1], h.wgrid[end], h.Nw)

	unscaled_itp_vf = interpolate(h.vf, (BSpline(Quadratic(Line())), NoInterp(), BSpline(Linear()), BSpline(Linear()), BSpline(Linear()), BSpline(Linear()), NoInterp(), NoInterp()), OnGrid())
	unscaled_itp_qᵍ  = interpolate(qᵍ_mat, (BSpline(Linear()), BSpline(Linear()), BSpline(Linear()), BSpline(Linear()), NoInterp(), NoInterp()), OnGrid())
	itp_vf = Interpolations.scale(unscaled_itp_vf, ωrange, 1:h.Nϵ, brange, μrange, σrange, wrange, 1:h.Nζ, 1:h.Nz)
	itp_qᵍ = Interpolations.scale(unscaled_itp_qᵍ, brange, μrange, σrange, wrange, 1:h.Nζ, 1:h.Nz)

	# Compute values
	vf, ϕa, ϕb, ϕc = opt_value(h, qʰ_mat, qᵍ_mat, wL_mat, T_mat, pC_mat, itp_qᵍ, itp_vf, resolve = resolve)

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