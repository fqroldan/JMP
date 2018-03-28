utility(h::Hank, c::Float64) = ifelse(c > 1e-10, c^(1.0 - h.ψ) / (1.0 - h.ψ), -1e10)

function EZ_G(h::Hank, v::Float64)
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

function EZ_T(h::Hank, Ev::Float64)
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

function value(h::Hank, sp::Float64, θp::Float64, itp_vf_s::Array{Interpolations.ScaledInterpolation{Float64,2,Interpolations.BSplineInterpolation{Float64,2,Array{Float64,2},Tuple{Interpolations.BSpline{Interpolations.Quadratic{Interpolations.Line}},Interpolations.NoInterp},Interpolations.OnGrid,(1, 0)},Tuple{Interpolations.BSpline{Interpolations.Quadratic{Interpolations.Line}},Interpolations.NoInterp},Interpolations.OnGrid,Tuple{StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}},UnitRange{Int64}}},2}, jϵ, jz, thres, RHS, qʰ, qᵍ, qᵍp, pC, jdefault)

	ap, bp, C = get_abc(RHS, h.ωmin, qʰ, qᵍ, pC, sp, θp)

	# Basis matrix for continuation values
	check, Ev, test, ut = 0., 0., 0, 0.

	if jdefault
		for (jzp, zpv) in enumerate(h.zgrid)
			for (jϵp, ϵpv) in enumerate(h.ϵgrid)
				prob =  h.Pz[jz, jzp] * h.Pϵ[jϵ, jϵp]

				# Reentry
				ζpv = 1
				R = h.κ + (1.0 - h.ρ) * qᵍp[jzp, 1]
				ωpv = ap + bp * R
				ωpv = min(h.ωmax, ωpv)
				v = itp_vf_s[jzp, 1][ωpv, jϵp]::Float64
				Ev += EZ_G(h, v) * prob * h.θ
				
				# Continue in default
				ζpv = 2
				R = (1.0 - h.ρ) * qᵍp[jzp, 2]
				ωpv = ap + bp * R
				ωpv = min(h.ωmax, ωpv)
				v = itp_vf_s[jzp, 2][ωpv, jϵp]::Float64
				Ev += EZ_G(h, v) * prob * (1.0 - h.θ)
				check += prob
			end
		end
	else
		for (jzp, zpv) in enumerate(h.zgrid)
			for (jϵp, ϵpv) in enumerate(h.ϵgrid)
				prob =  h.Pz[jz, jzp] * h.Pϵ[jϵ, jϵp]
				check += prob

				if zpv > thres
					ζpv = 1
					R = h.κ + (1.0 - h.ρ) * qᵍp[jzp, 1]
					ωpv = ap + bp * R
					ωpv = min(h.ωmax, ωpv)
					v = itp_vf_s[jzp, 1][ωpv, jϵp]::Float64
					Ev += EZ_G(h, v) * prob
				else
					ζpv = 2
					R = (1.0 - h.ρ) * qᵍp[jzp, 3]
					ωpv = ap + bp * R
					ωpv = min(h.ωmax, ωpv)
					v = itp_vf_s[jzp, 3][ωpv, jϵp]::Float64
					Ev += EZ_G(h, v) * prob
				end
			end
		end
	end

	isapprox(check, 1) || print_save("\nwrong expectation operator")

	# Compute value
	if h.EpsteinZin
		Tv = EZ_T(h, Ev)

		EZ_exp = (h.ψ-1.0)/h.ψ
		C > 1e-10? ut = C^(EZ_exp): ut = 1e-10

		vf = (1.0 - h.β) * ut + h.β * Tv^(EZ_exp)
		vf = vf^(1.0/EZ_exp)
		return vf
	else
		""" CHANGE THIS FOR GHH """
		ℓ = 0
		u = utility(h, C - ℓ)
		
		vf = (1.0 - h.β) * u + h.β * Ev
		return vf
	end
	Void
end

function solve_optvalue(h::Hank, guess::Vector, itp_vf_s, jϵ, jz, thres, RHS, qʰv, qᵍv, qᵍp, pCv, jdef, ωmax)


	optim_type = "multivariate"
	
	minθ = min(max(0.0, guess[2]-0.2), 0.8)
	maxθ = max(min(1.0, guess[2]+0.2), 0.2)

	ωspace = ωmax - qʰv*h.ωmin
	minω = min(max(qʰv*h.ωmin, guess[1] - 0.2*ωspace), qʰv*h.ωmin + 0.8 * ωspace)
	maxω = max(min(ωmax,       guess[1] + 0.2*ωspace), qʰv*h.ωmin + 0.2 * ωspace)

	ap, bp, cmax, fmax = 0., 0., 0., 0.
	if optim_type == "multivariate"
		guess[1] = max(min(guess[1], ωmax-1e-6), qʰv*h.ωmin+1e-6)
		guess[2] = max(min(guess[2], 1.0-1e-6), 1e-6)

		try
			res = Optim.optimize(
				x -> -value(h, x[1], x[2], itp_vf_s, jϵ, jz, thres, RHS, qʰv, qᵍv, qᵍp, pCv, jdef)
				, guess, [minω, minθ], [maxω, maxθ], Fminbox{LBFGS}())
		catch
			res = Optim.optimize(
				x -> -value(h, x[1], x[2], itp_vf_s, jϵ, jz, thres, RHS, qʰv, qᵍv, qᵍp, pCv, jdef)
				, guess, [minω, minθ], [maxω, maxθ], Fminbox{NelderMead}())
		end

		sp, θp = res.minimizer
		fmax = -res.minimum

		ap, bp, cmax = get_abc(RHS, h.ωmin, qʰv, qᵍv, pCv, sp, θp)
	end
	if optim_type == "sequential"
		function sub_value(h, sp, itp_vf_s, jϵ, jz, thres, RHS, qʰv, qᵍv, qᵍp, pCv, jdef; get_all::Bool=false)

			# minθ = 0.
			# maxθ = 1.

			res = Optim.optimize(
				θ -> -value(h, sp, θ, itp_vf_s, jϵ, jz, thres, RHS, qʰv, qᵍv, qᵍp, pCv, jdef),
				minθ, maxθ, GoldenSection(), rel_tol=h.tol_θ
				)

			if get_all
				θp = res.minimizer
				ap, bp, cmax = get_abc(RHS, h.ωmin, qʰv, qᵍv, pCv, sp, θp)
				return ap, bp, cmax, θp
			else
				fmax = -res.minimum
				return fmax
			end
		end

		res = Optim.optimize(
				sp -> -sub_value(h, sp, itp_vf_s, jϵ, jz, thres, RHS, qʰv, qᵍv, qᵍp, pCv, jdef),
					minω, maxω, GoldenSection(), rel_tol=h.tol_θ
				)
		sp = res.minimizer
		fmax = -res.minimum
		ap, bp, cmax, θp = sub_value(h, sp, itp_vf_s, jϵ, jz, thres, RHS, qʰv, qᵍv, qᵍp, pCv, jdef; get_all=true)		
	else
		curr_min = 1e10
		θp_grid = linspace(0,1,8)
		for θp in θp_grid
			res = Optim.optimize(
				sp -> -value(h, sp, θp, itp_vf_s, jϵ, jz, thres, RHS, qʰv, qᵍv, qᵍp, pCv, jdef),
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

		minimum(μpv) < minimum(h.μgrid) || maximum(μpv) > maximum(h.μgrid)? print_save("\nμ out of bounds at $([jb, jμ, jσ, jw, jζ, jz])"): Void
		minimum(σpv) < minimum(h.σgrid) || maximum(σpv) > maximum(h.σgrid)? print_save("\nσ out of bounds at $([jb, jμ, jσ, jw, jζ, jz])"): Void
		wpv < minimum(h.wgrid) || wpv > maximum(h.wgrid)? print_save("\nw out of bounds at $([jb, jμ, jσ, jw, jζ, jz])"): Void

		# Issuance should never be outside the bounds
		if bpv < minimum(h.bgrid) || bpv > maximum(h.bgrid)
			if isapprox(bpv, h.bgrid[1])
				bpv = h.bgrid[1]
			elseif isapprox(bpv, h.bgrid[end])
				bpv = h.bgrid[end]
			else
				print_save("\nb out of bounds at $([jb, jμ, jσ, jw, jζ, jz])")
			end
		end

		jdef = (h.ζgrid[jζ] != 1.0)

		qᵍp = Array{Float64}(h.Nz, h.Nϵ, 3)
		itp_vf_s = Array{Interpolations.ScaledInterpolation{Float64,2,Interpolations.BSplineInterpolation{Float64,2,Array{Float64,2},Tuple{Interpolations.BSpline{Interpolations.Quadratic{Interpolations.Line}},Interpolations.NoInterp},Interpolations.OnGrid,(1, 0)},Tuple{Interpolations.BSpline{Interpolations.Quadratic{Interpolations.Line}},Interpolations.NoInterp},Interpolations.OnGrid,Tuple{StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}},UnitRange{Int64}}}, 2}(h.Nz, 3)
		for (jzp, zpv) in enumerate(h.zgrid)
			qᵍp[jzp, 1] = itp_qᵍ[bpv, μpv[jzp, 1], σpv[jzp, 1], wpv, 1, jzp]				# Normal times
			qᵍp[jzp, 2] = itp_qᵍ[bpv, μpv[jzp, 2], σpv[jzp, 2], wpv, 2, jzp]				# Remain in default
			qᵍp[jzp, 3] = itp_qᵍ[(1.0 - h.ℏ)*bpv, μpv[jzp, 1], σpv[jzp, 1], wpv, 2, jzp]	# Enter default

			vf_mat = Array{Float64}(h.Nω, h.Nϵ, 3)
			for (jϵp, ϵpv) in enumerate(h.ϵgrid)
				for (jωp, ωpv) in enumerate(h.ωgrid)
					vf_mat[jωp, jϵp, 1] = itp_vf[ωpv, jϵp, bpv, μpv[jzp, 1], σpv[jzp, 1], wpv, 1, jzp]
					vf_mat[jωp, jϵp, 2] = itp_vf[ωpv, jϵp, bpv, μpv[jzp, 2], σpv[jzp, 2], wpv, 2, jzp]
					vf_mat[jωp, jϵp, 3] = itp_vf[ωpv, jϵp, (1.0-h.ℏ)*bpv, μpv[jzp, 1], σpv[jzp, 1], wpv, 2, jzp]
				end
				unscaled = interpolate(vf_mat[:,:,1], (BSpline(Quadratic(Line())), NoInterp()), OnGrid())
				itp_vf_s[jzp, 1] = Interpolations.scale(unscaled, linspace(h.ωgrid[1], h.ωgrid[end], h.Nω), 1:h.Nϵ)

				unscaled = interpolate(vf_mat[:,:,2], (BSpline(Quadratic(Line())), NoInterp()), OnGrid())
				itp_vf_s[jzp, 2] = Interpolations.scale(unscaled, linspace(h.ωgrid[1], h.ωgrid[end], h.Nω), 1:h.Nϵ)

				unscaled = interpolate(vf_mat[:,:,3], (BSpline(Quadratic(Line())), NoInterp()), OnGrid())
				itp_vf_s[jzp, 3] = Interpolations.scale(unscaled, linspace(h.ωgrid[1], h.ωgrid[end], h.Nω), 1:h.Nϵ)
			end
		end

		for (jϵ, ϵv) in enumerate(h.ϵgrid), (jω, ωv) in enumerate(h.ωgrid)

			RHS = ωv + wv * exp(ϵv) - Tv

			ap, bp, cmax, fmax = 0., 0., 0., 0.
			ag, bg = h.ϕa[jω, jϵ, jb, jμ, jσ, jw, jζ, jz], h.ϕb[jω, jϵ, jb, jμ, jσ, jw, jζ, jz]

			ωg = qʰv * ag + qᵍv * bg
			θg = qʰv * (ag - h.ωmin) / (ωg - qʰv*h.ωmin)
			# print_save("a,b,s,θ = $([ag, bg, ωg, θg])")
			ωmax = RHS - 1e-10
			if ωg > ωmax
				ωg = max(ωmax - 1e-2, 0)
			end
			if ωmax < qʰv * h.ωmin
				print_save("\nCan't afford positive consumption at $([jb, jμ, jσ, jw, jζ, jz]) with w*Lᵈ=$(round(wv,2)), T=$(round(Tv,2))")
				ap, bp, cmax = h.ωmin, 0., 1e-10
				fmax = value(h, qʰv*ap, 0., itp_vf, jϵ, jz, thres, RHS, qʰv, qᵍv, qᵍp, pCv, jdef)
			end
			isapprox(θg, 1) && θg > 1? θg = 1.0: Void
			
			if resolve
				# θg = 1.0
				guess = [ωg, θg]

				ap, bp, cmax, fmax = solve_optvalue(h, guess, itp_vf_s, jϵ, jz, thres, RHS, qʰv, qᵍv, qᵍp, pCv, jdef, ωmax)
			else
				ap = h.ϕa[jω, jϵ, jb, jμ, jσ, jw, jζ, jz]
				bp = h.ϕb[jω, jϵ, jb, jμ, jσ, jw, jζ, jz]
				cmax = h.ϕc[jω, jϵ, jb, jμ, jσ, jw, jζ, jz]
				fmax = value(h, ωg, θg, itp_vf_s, jϵ, jz, thres, RHS, qʰv, qᵍv, qᵍp, pCv, jdef)
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
	# Interpolate the value function and debt prices
	itp_vf = make_itps(h, h.vf; agg=false)
	itp_qᵍ = make_itps(h, h.qᵍ; agg=true)

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