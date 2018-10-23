function utility(h::Hank, c::Float64)
	if c > 1e-10
		if h.ψ != 1.0
			u = c^(1.0 - h.ψ) / (1.0 - h.ψ)
		else
			u = log(c)
		end
	else
		u = -1e10
	end
end

function EZ_G(h::Hank, v::Float64)
	if h.EpsteinZin
		if h.γ != 1.0
			return v^(1.0-h.γ)
		else
			return log(v)
		end
	else
		return v
	end
end

function EZ_T(h::Hank, Ev::Float64)
	if h.γ != 1.0
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
			u[jc] = cv.^(-1.0/h.ψ)
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

function get_abec(RHS::Float64, ωmin::Float64, qʰ::Float64, qᵍ::Float64, pC::Float64, sp::Float64, θa::Float64)

	θe = 0.0
	qᵉ = 1.0

	""" Recovers private and public debt purchases and consumption from savings decisions """
	ap = ωmin + θa * (sp - qʰ*ωmin) / qʰ
	bp = (1.0-θa) * (1.0-θe) * (sp - qʰ*ωmin) / qᵍ
	ep = (1.0-θa) *    θe    * (sp - qʰ*ωmin) / qᵉ
	C  = (RHS - sp) / pC

	return ap, bp, ep, C
end

function eval_itp_vf(itp_vf_s::Arr_itp_VF, ωpv::Float64, jϵp::Int64, jzp::Int64, jj::Int64)
	itp_obj = itp_vf_s[jzp, jj]
	v = itp_obj[ωpv, jϵp]
	return v
end

function value(h::Hank, sp::Float64, θa::Float64, itp_vf_s::Arr_itp_VF, jϵ, jz, thres, exp_rep, RHS, qʰ, qᵍ, qᵍp, profits, pC, jdefault)
# function value(h::Hank, sp::Float64, θa::Float64, itp_vf_s::Array{Interpolations.GriddedInterpolation{Float64,2,Float64,Tuple{Interpolations.Gridded{Interpolations.Linear},Interpolations.NoInterp},Tuple{Array{Float64,1},Array{Int64,1}},0}}, jϵ, jz, thres, exp_rep, RHS, qʰ, qᵍ, qᵍp, profits, pC, jdefault)

	ap, bp, ep, C = get_abec(RHS, h.ωmin, qʰ, qᵍ, pC, sp, θa)

	itp_s = true

	# Basis matrix for continuation values
	check, Ev, test, ut = 0., 0., 0, 0.

	if jdefault
		for (jzp, zpv) in enumerate(h.zgrid)
			for (jϵp, ϵpv) in enumerate(h.ϵgrid)
				prob =  h.Pz[jz, jzp] * h.Pϵ[jϵ, jϵp]

				# Reentry
				ζpv = 1
				Rb = h.κ + (1. - h.ρ) * qᵍp[jzp, 1]
				# Re = profits[jzp, 1]
				ωpv = ap + bp * Rb# + ep * Re
				if ωpv < h.ωmin
					Ev += prob * h.θ * 1e-10
				else
					ωpv = min(h.ωmax, ωpv)
					v = eval_itp_vf(itp_vf_s, ωpv, jϵp, jzp, 1)
					Ev += EZ_G(h, v) * prob * h.θ
				end

				# Continue in default
				ζpv = 2
				Rb = (1. - h.ρ) * qᵍp[jzp, 2]
				# Re = profits[jzp, 2]
				ωpv = ap + bp * Rb# + ep * Re
				if ωpv < h.ωmin
					Ev += prob * (1. - h.θ) * 1e-10
				else
					ωpv = min(h.ωmax, ωpv)
					v = eval_itp_vf(itp_vf_s, ωpv, jϵp, jzp, 2)
					Ev += EZ_G(h, v) * prob * (1. - h.θ)
				end
				check += prob
			end
		end
	else
		for (jzp, zpv) in enumerate(h.zgrid)
			for (jϵp, ϵpv) in enumerate(h.ϵgrid)
				prob =  h.Pz[jz, jzp] * h.Pϵ[jϵ, jϵp]
				check += prob

				# Repayment
				ζpv = 1
				Rb = h.κ + (1. - h.ρ) * qᵍp[jzp, 1]
				# Re = profits[jzp, 1]
				ωpv = ap + bp * Rb# + ep * Re
				if ωpv < h.ωmin
					Ev += prob * exp_rep[jzp] * 1e-10
				else
					ωpv = min(h.ωmax, ωpv)
					v = eval_itp_vf(itp_vf_s, ωpv, jϵp, jzp, 1)
					Ev += EZ_G(h, v) * prob * exp_rep[jzp]
				end
				# Default
				ζpv = 2
				Rb = (1.0 - h.ρ) * (1.0 - h.ℏ) * qᵍp[jzp, 2]
				# Re = profits[jzp, 3]
				ωpv = ap + bp * Rb# + ep * Re
				if ωpv < h.ωmin
					Ev += prob * (1. - exp_rep[jzp]) * 1e-10
				else
					ωpv = min(h.ωmax, ωpv)
					v = eval_itp_vf(itp_vf_s, ωpv, jϵp, jzp, 2)
					Ev += EZ_G(h, v) * prob * (1. - exp_rep[jzp])
				end
			end
		end
	end

	isapprox(check, 1) || print_save("\nwrong expectation operator")

	""" CHANGE THIS FOR GHH """
	# Compute value
	if h.EpsteinZin
		if Ev < -1
			throw(error("Something wrong in the choice of ωpv"))
		elseif Ev < 0
			Ev = 1e-6
		end
		Tv = EZ_T(h, Ev)

		if h.ψ != 1.0
			EZ_exp = (h.ψ-1.0)/h.ψ
			C > 1e-10 ? ut = C^(EZ_exp) : ut = 1e-10

			vf = (1.0 - h.β) * ut + h.β * Tv^(EZ_exp)
			vf = vf^(1.0/EZ_exp)
		else
			C > 1e-10 ? ut = C^(1.0-h.β) : ut = 1e-10
			vf = ut * Tv^(h.β) # This is the same as saying that vf = exp( (1.0-h.β)*log(c) + h.β * log(Tv) )
		end
		return vf
	else
		ℓ = 0
		u = utility(h, C - ℓ)
		vf = (1.0 - h.β) * u + h.β * Ev
		return vf
	end
	Void
end

function solve_optvalue(h::Hank, guess::Vector, itp_vf_s, jϵ, jz, thres, exp_rep, RHS, qʰv, qᵍv, qᵍp, profits, pCv, jdef, ωmax)


	optim_type = "multivariate"

	minθ = min(max(0.0, guess[2]-0.2), 0.8)
	maxθ = max(min(1.0, guess[2]+0.2), 0.2)

	if minθ < 0.1
		minθ = 0.1
		maxθ = max(maxθ, 0.3)
	end

	ωspace = ωmax - qʰv*h.ωmin
	minω = min(max(qʰv*h.ωmin, guess[1] - 0.2*ωspace), qʰv*h.ωmin + 0.8 * ωspace)
	maxω = max(min(ωmax,       guess[1] + 0.2*ωspace), qʰv*h.ωmin + 0.2 * ωspace)

	if maxω - minω < 1e-4
		maxω = minω + 1e-4
	end

	# minθ, maxθ = 0., 1.
	# minω, maxω = qʰv*h.ωmin, ωmax

	ap, bp, ep, cmax, fmax = 0., 0., 0., 0., 0.
	if optim_type == "sequential"
		function sub_value(h, sp, itp_vf_s, jϵ, jz, thres, RHS, qʰv, qᵍv, qᵍp, profits, pCv, jdef; get_all::Bool=false)

			# minθ = 0.
			# maxθ = 1.

			res = Optim.optimize(
				θ -> -value(h, sp, θ, itp_vf_s, jϵ, jz, thres, exp_rep, RHS, qʰv, qᵍv, qᵍp, profits, pCv, jdef),
				minθ, maxθ, GoldenSection(), rel_tol=h.tol_θ
				)

			if get_all
				θa = res.minimizer
				ap, bp, ep, cmax = get_abec(RHS, h.ωmin, qʰv, qᵍv, pCv, sp, θa)
				return ap, bp, ep, cmax, θa
			else
				fmax = -res.minimum
				return fmax
			end
		end

		res = Optim.optimize(
				sp -> -sub_value(h, sp, itp_vf_s, jϵ, jz, thres, RHS, qʰv, qᵍv, qᵍp, profits, pCv, jdef),
					minω, maxω, GoldenSection(), rel_tol=h.tol_θ
				)
		sp = res.minimizer
		fmax = -res.minimum
		ap, bp, ep, cmax, θa = sub_value(h, sp, itp_vf_s, jϵ, jz, thres, RHS, qʰv, qᵍv, qᵍp, profits, pCv, jdef; get_all=true)
	elseif optim_type == "multivariate"

		guess[1] = max(min(guess[1], maxω-1e-6), minω+1e-6)
		guess[2] = max(min(guess[2], maxθ-1e-6), minθ+1e-6)

		try
			inner_opt = LBFGS(;linesearch=LineSearches.HagerZhang(linesearchmax=100))
			res = Optim.optimize(
				x -> -value(h, x[1], x[2], itp_vf_s, jϵ, jz, thres, exp_rep, RHS, qʰv, qᵍv, qᵍp, profits, pCv, jdef)
				, [minω, minθ], [maxω, maxθ], guess, Fminbox(inner_opt))
			if Optim.converged(res)
			else
				throw("Main algorithm failed")
			end
		catch
			if minω < guess[1] < maxω && minθ < guess[2] < maxθ
			else
				print_save("\nWARNING: MAYBE PROBLEM WITH BOUNDS")
				print_save("\n[$(minω), $(guess[1]), $(maxω)]")
				print_save("\n[$(minθ), $(guess[2]), $(maxθ)]")
			end
			res = Optim.optimize(
				x -> -value(h, x[1], x[2], itp_vf_s, jϵ, jz, thres, exp_rep, RHS, qʰv, qᵍv, qᵍp, profits, pCv, jdef)
				, [minω, minθ], [maxω, maxθ], guess, Fminbox(NelderMead()))
		end

		sp, θa = res.minimizer
		fmax = -res.minimum

		ap, bp, ep, cmax = get_abec(RHS, h.ωmin, qʰv, qᵍv, pCv, sp, θa)
	else
		curr_min = 1e10
		θa_grid = linspace(0,1,8)
		for θa in θa_grid
			res = Optim.optimize(
				sp -> -value(h, sp, θa, itp_vf_s, jϵ, jz, thres, exp_rep, RHS, qʰv, qᵍv, qᵍp, profits, pCv, jdef),
					qʰv*h.ωmin, ωmax, GoldenSection()
				)
			if res.minimum < curr_min
				sp = res.minimizer
				ap, bp, ep, cmax = get_abec(RHS, h.ωmin, qʰv, qᵍv, pCv, sp, θa)
				fmax = -res.minimum
				curr_min = res.minimum
			end
		end
	end

	return ap, bp, ep, cmax, fmax
end


function opt_value(h::Hank, qʰ_mat, qᵍ_mat, wL_mat, T_mat, pC_mat, Π_mat, itp_qᵍ, itp_vf; resolve::Bool = true, verbose::Bool=true)

	vf = SharedArray{Float64}(size(h.vf))
	ϕa = SharedArray{Float64}(size(h.ϕa))
	ϕb = SharedArray{Float64}(size(h.ϕb))
	ϕe = SharedArray{Float64}(size(h.ϕe))
	ϕc = SharedArray{Float64}(size(h.ϕc))
	warnc0 = SharedArray{Float64}(h.Nb, h.Nμ, h.Nσ, h.Nw, h.Nζ, h.Nz)
	warnc0 *= 0.
	@sync @parallel for js in 1:size(h.Jgrid,1)
		jb = h.Jgrid[js, 1]
		jμ = h.Jgrid[js, 2]
		jσ = h.Jgrid[js, 3]
		jw = h.Jgrid[js, 4]
		jζ = h.Jgrid[js, 5]
		jz = h.Jgrid[js, 6]

		qʰv = qʰ_mat[jb, jμ, jσ, jw, jζ, jz]
		qᵍv = qᵍ_mat[jb, jμ, jσ, jw, jζ, jz]
		wL  = wL_mat[jb, jμ, jσ, jw, jζ, jz]
		Tv  = T_mat[jb, jμ, jσ, jw, jζ, jz]
		pCv = pC_mat[jb, jμ, jσ, jw, jζ, jz]
		profits = Π_mat[jb, jμ, jσ, jw, jζ, jz]

		bpv = h.issuance[js]
		μpv = h.μ′[js,:,:]
		σpv = h.σ′[js,:,:]
		wpv = h.w′[js]
		thres = h.def_thres[js]

		rep_mat = reshape(h.repay, h.Nb, h.Nμ, h.Nσ, h.Nw, h.Nζ, h.Nz, h.Nz)
		exp_rep = rep_mat[jb, jμ, jσ, jw, jζ, jz, :]

		if verbose
			minimum(μpv) < minimum(h.μgrid) || maximum(μpv) > maximum(h.μgrid) ? print_save("\nμ out of bounds at $([jb, jμ, jσ, jw, jζ, jz])") : Void
			minimum(σpv) < minimum(h.σgrid) || maximum(σpv) > maximum(h.σgrid) ? print_save("\nσ out of bounds at $([jb, jμ, jσ, jw, jζ, jz])") : Void
			bpv - minimum(h.bgrid) < -1e-4 || bpv - maximum(h.bgrid) > 1e-4 ? print_save("\nb = $(round(bpv,6)) out of bounds at $([jb, jμ, jσ, jw, jζ, jz])") : Void
			wpv < minimum(h.wgrid) || wpv > maximum(h.wgrid) ? print_save("\nw out of bounds at $([jb, jμ, jσ, jw, jζ, jz])") : Void
		end


		jdef = (h.ζgrid[jζ] != 1.0)

		qᵍp = Array{Float64}(h.Nz, 2)
		itp_vf_s = Arr_itp_VF(h.Nz, 2)
		for (jzp, zpv) in enumerate(h.zgrid)
			qᵍp[jzp, 1] = itp_qᵍ[bpv, μpv[jzp, 1], σpv[jzp, 1], wpv, 1, jzp]
			if jdef
				qᵍp[jzp, 2] = itp_qᵍ[bpv, μpv[jzp, 2], σpv[jzp, 2], wpv, 2, jzp]
			else
				qᵍp[jzp, 2] = itp_qᵍ[(1.0 - h.ℏ)*bpv, μpv[jzp, 2], σpv[jzp, 2], wpv, 2, jzp]
			end

			vf_mat = Array{Float64}(h.Nω, h.Nϵ, 3)
			for (jϵp, ϵpv) in enumerate(h.ϵgrid)
				for (jωp, ωpv) in enumerate(h.ωgrid)
					vf_mat[jωp, jϵp, 1] = itp_vf[ωpv, jϵp, bpv, μpv[jzp, 1], σpv[jzp, 1], wpv, 1, jzp]
					if jdef
						vf_mat[jωp, jϵp, 2] = itp_vf[ωpv, jϵp, bpv, μpv[jzp, 2], σpv[jzp, 2], wpv, 2, jzp]
					else
						vf_mat[jωp, jϵp, 2] = itp_vf[ωpv, jϵp, (1.0-h.ℏ)*bpv, μpv[jzp, 1], σpv[jzp, 1], wpv, 2, jzp]
					end
				end

				knots = (h.ωgrid, 1:h.Nϵ)
				for jj in 1:2
					itp_vf_s[jzp, jj] = interpolate(knots, vf_mat[:,:,jj], (Gridded(Linear()), NoInterp()))
					# unscaled = interpolate(vf_mat[:,:,jj], (BSpline(Quadratic(Line())), NoInterp()), OnGrid())
					# itp_vf_s[jzp, jj] = Interpolations.scale(unscaled, linspace(h.ωgrid[1], h.ωgrid[end], h.Nω), 1:h.Nϵ)
				end
			end
		end

		# @code_warntype v = eval_itp_vf(itp_vf_s, 10.3, 2, 2, 1)
		adjustment = sum(h.λϵ.*exp.(h.ϵgrid))
		for (jϵ, ϵv) in enumerate(h.ϵgrid), (jω, ωv) in enumerate(h.ωgrid)

			RHS = ωv + wL * exp(ϵv) - Tv + profits * exp(ϵv) / adjustment

			ap, bp, ep, cmax, fmax = 0., 0., 0., 0., 0.
			ag, bg = h.ϕa[jω, jϵ, jb, jμ, jσ, jw, jζ, jz], h.ϕb[jω, jϵ, jb, jμ, jσ, jw, jζ, jz]

			ωg = qʰv * ag + qᵍv * bg
			ωmax = RHS - 1e-10
			if ωg > ωmax
				ωg = max(ωmax - 1e-2, 0)
			end
			if ωg - qʰv*h.ωmin > 1e-4
				θg = qʰv * (ag - h.ωmin) / (ωg - qʰv*h.ωmin)
			else
				θg = 1.
			end
			isapprox(θg, 1) && θg > 1 ? θg = 1.0 : Void

			if resolve && ωmax > qʰv * h.ωmin
				guess = [ωg, θg]

				ap, bp, ep, cmax, fmax = solve_optvalue(h, guess, itp_vf_s, jϵ, jz, thres, exp_rep, RHS, qʰv, qᵍv, qᵍp, profits, pCv, jdef, ωmax)
				if cmax < 0
					warn("c = $cmax")
					ap, bp, ep, cmax = h.ωmin, 0., 0., 1e-8
					fmax = 1e-10
				end
			elseif ωmax < qʰv * h.ωmin
				if verbose
					# print_save("\nCan't afford positive consumption at $([jb, jμ, jσ, jw, jζ, jz]) with w*Lᵈ=$(round(wL,2)), T=$(round(Tv,2))")
					warnc0[jb, jμ, jσ, jw, jζ, jz] = 1.
				end
				ap, bp, ep, cmax = h.ωmin, 0., 0., 1e-8
				fmax = 1e-10
			else
				ap = h.ϕa[jω, jϵ, jb, jμ, jσ, jw, jζ, jz]
				bp = h.ϕb[jω, jϵ, jb, jμ, jσ, jw, jζ, jz]
				cmax = h.ϕc[jω, jϵ, jb, jμ, jσ, jw, jζ, jz]
				fmax = value(h, ωg, θg, itp_vf_s, jϵ, jz, thres, exp_rep, RHS, qʰv, qᵍv, qᵍp, profits, pCv, jdef)
			end
			!isnan(fmax) || print_save("\nWARNING: NaN in value function at (ap, bp, c) = ($(round(ap, 2)), $(round(bp, 2)), $(cmax))")

			ϕa[jω, jϵ, jb, jμ, jσ, jw, jζ, jz] = ap
			ϕb[jω, jϵ, jb, jμ, jσ, jw, jζ, jz] = bp
			ϕe[jω, jϵ, jb, jμ, jσ, jw, jζ, jz] = ep
			ϕc[jω, jϵ, jb, jμ, jσ, jw, jζ, jz] = cmax
			vf[jω, jϵ, jb, jμ, jσ, jw, jζ, jz] = fmax
		end
	end

	return vf, ϕa, ϕb, ϕe, ϕc, warnc0
end

function bellman_iteration!(h::Hank, qʰ_mat, qᵍ_mat, wL_mat, T_mat, pC_mat, Π_mat; resolve::Bool=true, verbose::Bool=true)
	t1 = time()
	# Interpolate the value function
	itp_vf = make_itp(h, h.vf; agg=false)
	itp_qᵍ = make_itp(h, h.qᵍ; agg=true)
	# itp_vf, itp_qᵍ = make_itps(h, h.vf, qᵍ_mat)
	# print_save("\ninterp in $(time()-t1)")

	# Compute values
	t1 = time()
	vf, ϕa, ϕb, ϕe, ϕc, warnc0 = opt_value(h, qʰ_mat, qᵍ_mat, wL_mat, T_mat, pC_mat, Π_mat, itp_qᵍ, itp_vf, resolve = resolve, verbose = verbose)
	# print_save("\nopt in $(time()-t1)")

	t1 = time()
	sum(isnan.(vf)) > 0 ? print_save("\n$(sum(isnan.(vf))) NaNs found in vf") : Void
	sum(isnan.(ϕa)) > 0 ? print_save("$(sum(isnan.(ϕa))) NaNs found in ϕa") : Void
	sum(isnan.(ϕb)) > 0 ? print_save("$(sum(isnan.(ϕb))) NaNs found in ϕb") : Void
	sum(isnan.(ϕc)) > 0 ? print_save("$(sum(isnan.(ϕc))) NaNs found in ϕc") : Void

	# Store results in the type
	h.ϕa = ϕa
	h.ϕb = ϕb
	h.ϕe = ϕe
	h.ϕc = ϕc
	h.vf = vf
	# print_save("\nsave in $(time()-t1)")

	return warnc0
end
