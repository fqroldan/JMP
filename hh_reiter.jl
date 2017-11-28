@everywhere begin

using QuantEcon, BasisMatrices, Interpolations, Optim, NLopt, MINPACK, Plots, LaTeXStrings, Distributions
pyplot()

type Hank
	# Utility parameters
	β::Float64
	γ::Float64
	θ::Float64
	χ::Float64
	ρϵ::Float64
	σϵ::Float64
	Ξ::Float64

	# Debt parameters
	ρ::Float64
	κ::Float64

	# NK parameters
	Φπ::Float64
	ΦL::Float64
	η::Float64
	elast::Float64

	# Policy functions
	gc::Vector{Float64}
	gω::Vector{Float64}

	gc_ext::Array{Float64, 8}
	gω_ext::Array{Float64, 8}

	# Coefficients of value and policy functions
	cc::Vector{Float64}
	cω::Vector{Float64}
	cv::Vector{Float64}
	ce::Vector{Float64}

	# Exogenous state
	ρz::Float64
	σz::Float64

	# Grid points
	Nω::Int64
	Nϵ::Int64
	Nb::Int64
	Nμ::Int64
	Nσ::Int64
	Nz::Int64
	Ns::Int64
	Nω_fine::Int64

	# Transition matrices
	Pϵ::Matrix{Float64}
	Pz::Matrix{Float64}

	Ps::Matrix{Float64}

	# Distributions
	λ::Vector{Float64}
	λϵ::Vector{Float64}
	ℏ::Float64
	thr_def::Float64

	# Parameters of the ω grid
	curv::Float64
	order::Int64
	ωmin::Float64
	ωmax::Float64

	ωgrid0::Vector{Float64}
	ωgrid::Vector{Float64}
	ϵgrid::Vector{Float64}
	bgrid::Vector{Float64}
	μgrid::Vector{Float64}
	σgrid::Vector{Float64}
	zgrid::Vector{Float64}
	s::Matrix{Float64}

	# Extra grids for prices
	qgrid::Vector{Float64}
	wgrid::Vector{Float64}

	# Collocation objects
	basis::Basis
	bs::BasisMatrix
	Φ::SparseMatrixCSC
	dΦ::SparseMatrixCSC
	Emat::SparseMatrixCSC
	Φnotω::SparseMatrixCSC

	ωgrid_fine::Vector{Float64}
	snodes::Array{Float64, 2}

	# Forecasting rules
	μ′::Vector{Float64}
	σ′::Vector{Float64}

	# Functions of the state
	Aplus::Vector{Float64}
	Aminus::Vector{Float64}
	debt_repay::Vector{Float64}
	MF_rS::Vector{Float64}
	τ::Float64
	lump_sum::Vector{Float64}
	issuance_B::Vector{Float64}
	spending::Vector{Float64}
	wage::Vector{Float64}
	Ld::Vector{Float64}
	# sdf_vec::Vector{Float64}
	debtprice::Vector{Float64}
	q::Vector{Float64}
	inflation::Vector{Float64}
	Πstar::Float64
	i_star::Float64

	ξg::Vector{Float64}
	ξf::Vector{Float64}
end


function Hank(;	β = (1/1.06)^(1/4),
				γ = 2.,
				θ = 1.,
				χ = 2.,
				τ = 0.3,
				ωmax = 7.5,
				curv = .4,
				order = 3,
				Nω = 6,
				Nω_fine = 1000,
				ρϵ = 0.9,
				σϵ = 0.1,
				Nϵ = 3,
				Nμ = 5,
				Nσ = 5,
				Nb = 5,
				ρz = 0.92,
				σz = 0.007,
				Nz = 5,
				ℏ = .5,
				thr_def = -0.1,
				Nq = 7,
				Nw = 3)
	# Prepare discretized processes	
	if -1 < ρz < 1
		z_chain = tauchen(Nz, ρz, σz, 0, 1)
		Pz = z_chain.p
		zgrid = exp.(z_chain.state_values)
	elseif ρz == 1.
		Pz = eye(Nz)
		zgrid = exp.(collect(linspace(-0.05, 0.05, Nz)))
	else
		throw(error("Parameters of z shock badly defined"))
	end
	ϵ_chain = tauchen(Nϵ, ρϵ, σϵ, 0, 1)
	Pϵ = ϵ_chain.p
	ϵgrid = exp.(ϵ_chain.state_values)

	λϵ = (Pϵ^100)[1,:]

	Ξ = dot(ϵgrid.^(1/χ), λϵ)^χ

	θ = (1-τ) * Ξ

	# Grids for endogenous aggregate states
	bgrid = collect(linspace(0.7, 0.9, Nb))
	μgrid = collect(linspace(0.1, 0.3, Nμ))
	σgrid = collect(linspace(1.0, 2.0, Nσ))

	# Debt parameters
	ρ = 0.05 # Target average maturity of 7 years: ~0.05 at quarterly freq
	κ = ρ

	Φπ = 5.0
	ΦL = 1.0

	η  = 250.0
	elast = 6.0

	# Transitions
	Ps  = Array{Float64}(Nb*Nμ*Nσ*Nz, Nb*Nμ*Nσ*Nz)

	# Prepare grid for cash in hand.
	""" Make sure that the lowest ω point affords positive c at the worst prices """
	ωmin	= -1.0
	ωgrid0	= linspace(0., (ωmax-ωmin)^curv, Nω).^(1/curv)
	ωgrid0	= ωgrid0 + ωmin

	ωgrid_fine	= linspace(0., (ωmax-ωmin)^curv, Nω_fine).^(1/curv)
	ωgrid_fine	= ωgrid_fine + ωmin

	snodes = [kron(ones(Nϵ,), ωgrid_fine) kron(ϵgrid, ones(Nω_fine,))]

	# Define the basis over the state variables
	basis = Basis(SplineParams(ωgrid0, 0, order),
				  LinParams(ϵgrid, 0),
				  LinParams(bgrid, 0),
				  LinParams(μgrid, 0),
				  LinParams(σgrid, 0),
				  LinParams(zgrid, 0))
	s, (ωgrid, ϵgrid, bgrid, μgrid, σgrid, zgrid) = nodes(basis)
	Nω, Ns = size(ωgrid, 1), size(s, 1)

	# Compute the basis matrix and expectations matrix
	bs = BasisMatrix(basis, Direct(), s, [0 0 0 0 0 0])
	Φ = convert(Expanded, bs).vals[1]

	dΦ = BasisMatrix(basis, Expanded(), s, [1 0 0 0 0 0]).vals[1]

	# Save the parts of the interpolation that are not 'ω'
	Φϵ = bs.vals[2]
	Φb = bs.vals[3]
	Φμ = bs.vals[4]
	Φσ = bs.vals[5]
	Φz = bs.vals[6]
	Φnotω = row_kron(Φz, row_kron(Φσ, row_kron(Φμ, row_kron(Φb, Φϵ))))

	Emat = kron(Ps, kron(Pϵ, speye(Nω))) * Φ

	cc = ones(Ns,)
	cω = ones(Ns,)

	gc = ones(Ns,)
	gω = ones(Ns,)

	cv = ones(Ns,)
	ce = ones(Ns,)

	λ = ones(Nω_fine*Nϵ)
	λ = λ/sum(λ)

	qgrid = collect(linspace(.85, 1.1, Nq))
	wgrid = collect(linspace(.9, 1.1, Nw))

	gc_ext = zeros(Nω, Nϵ, Nb, Nμ, Nσ, Nz, Nq, Nw)
	gω_ext = zeros(Nω, Nϵ, Nb, Nμ, Nσ, Nz, Nq, Nw)

	haircut = ℏ * (log.(s[:, 6]) .< thr_def)
	debt_repay = 1-haircut

	r_quart = 0.5/100
	R = (1+r_quart)^(1/4)
	r = ones(Ns,) * (R-1.)

	MF_rS = r

	lump_sum = ones(Ns,) * 0.05

	issuance_B = s[:, 3] + 0.1 * log.(s[:, 6])
	μ′ = Array{Float64}(Nb, Nμ, Nσ, Nz)
	for (jμ, μv) in enumerate(μgrid)
		μ′[:,jμ,:,:] = μv
	end
	μ′ = reshape(μ′, Nb*Nμ*Nσ*Nz)
	σ′ = Array{Float64}(Nb, Nμ, Nσ, Nz)
	for (jσ, σv) in enumerate(σgrid)
		σ′[:,:,jσ,:] = σv
	end
	σ′ = reshape(σ′, Nb*Nμ*Nσ*Nz)

	spending = 0.2 - 0.05 * log.(s[:, 6])

	wage = s[:, 6]
	Ld 	 = 1 * ones(Ns,)

	debtprice = 0.95 * ones(Ns,)
	q = 0.99 * ones(Ns,)

	Πstar = 1.02^(1/4)
	i_star = (1 / 100 + Πstar^4)^(1/4) - 1
	inflation = Πstar * ones(Ns,)

	ξg = zeros(Ns,)
	ξf = zeros(Ns,)

	function compute_grosspositions(μ,σ)
		valplus, valminus, sum_prob = 0.,0.,0.
		for (jϵ, ϵv) in enumerate(ϵgrid)
			for jω = 1:length(ωgrid_fine)-1
				ωv  = ωgrid_fine[jω]
				ω1v = ωgrid_fine[jω+1]
				ωmv = 0.5*(ωv+ω1v)

				prob = pdf(LogNormal(μ, σ), ωmv-ωmin) * λϵ[jϵ] * (ω1v - ωv)

				valplus  += prob * max(ωmv, 0)
				valminus += prob * max(-ωmv, 0)
				sum_prob += prob
			end
		end
		Aup  = valplus  / sum_prob
		Adown = valminus / sum_prob
		
		return Aup, Adown
	end

	Aplus_mat, Aminus_mat = zeros(Nμ,Nσ), zeros(Nμ,Nσ)
	for (jσ, σv) in enumerate(σgrid), (jμ, μv) in enumerate(μgrid)
		Aplus_mat[jμ, jσ], Aminus_mat[jμ, jσ] = compute_grosspositions(μv, σv)
	end

	js = 0
	Aplus, Aminus = zeros(Nb*Nμ*Nσ*Nz), zeros(Nb*Nμ*Nσ*Nz)
	for (jz, zv) in enumerate(zgrid), (jσ,σv) in enumerate(σgrid), (jμ,μv) in enumerate(μgrid), (jb,bv) in enumerate(bgrid)
		js += 1
		Aplus[js], Aminus[js] = Aplus_mat[jμ,jσ], Aminus_mat[jμ,jσ]
	end
	
	return Hank(β, γ, θ, χ, ρϵ, σϵ, Ξ, ρ, κ, Φπ, ΦL, η, elast, gc, gω, gc_ext, gω_ext, cc, cω, cv, ce, ρz, σz, Nω, Nϵ, Nb, Nμ, Nσ, Nz, Ns, Nω_fine, Pϵ, Pz, Ps, λ, λϵ, ℏ, thr_def, curv, order, ωmin, ωmax, ωgrid0, ωgrid, ϵgrid, bgrid, μgrid, σgrid, zgrid, s, qgrid, wgrid, basis, bs, Φ, dΦ, Emat, Φnotω, ωgrid_fine, snodes, μ′, σ′, Aplus, Aminus, debt_repay, MF_rS, τ, lump_sum, issuance_B, spending, wage, Ld, debtprice, q, inflation, Πstar, i_star, ξg, ξf)
end

function _unpackstatefs(h::Hank)
	R = 1 + h.MF_rS
	T = h.lump_sum
	ℓ = h.θ^(-1/h.χ) * (h.s[:,2] .* h.wage .* (1 - h.τ)).^((1+h.χ)/h.χ)
	Rep = h.debt_repay
	q = h.q
	Π = h.inflation

	return R, T, ℓ, Rep, q, Π
end

function utility(h::Hank, c_vec::Vector)
	u = zeros(size(c_vec))
	for (jc, cv) in enumerate(c_vec)
		if cv > 0
			if h.γ == 1.
				u[jc] = log.(cv)
			else
				u[jc] = (cv.^(1-h.γ) - 1)/(1-h.γ)
			end
		else
			u[jc] = -1e10
		end
	end
	if length(c_vec) == 1
		return u[1]
	else
		return u
	end
end

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

function value(h::Hank, gω::Vector{Float64}, s, RHS::Vector{Float64}, ℓ, q; solved::Bool=false)
	h.gc = (RHS - gω) .* q

	Ut = utility(h, h.gc - ℓ/(1+h.χ) )

	Φωp	= BasisMatrix(h.basis[1], Direct(), gω, 0).vals[1]
	Φ_new = row_kron(h.Φnotω, Φωp)

	if solved
		h.cc = h.Φ\h.gc
		return Φ_new
	end

	# Compute value
	vf = Ut + h.β * Φ_new*h.ce 
	return vf
end

function opt_value!(h::Hank, s::Matrix{Float64}, R, T, ℓ, q, Π; newton::Bool=false, resolve::Bool = true)
	rS = (R - 1) .* (s[:,1] .>= 0)
	# Decide savings
	lower_bound = ones(h.Ns,) * h.ωmin				# Borrowing constraint
	upper_bound = ( (1+rS).* s[:,1]./Π + ℓ - T )./q	# Budget constraint; c ≧ 0
	if minimum(upper_bound-lower_bound) < 0
		warn("Budget set empty. $(minimum(upper_bound-lower_bound))")
	end
	vf = zeros(h.gω)
	if resolve
		h.gω, vf = golden_method((gω -> value(h, gω, s, upper_bound, ℓ, q)), lower_bound, upper_bound)
		if minimum(h.gω) > h.ωmin
			Void
		elseif isapprox(minimum(h.gω), h.ωmin)
			h.gω = max.(h.gω, h.ωmin)
			h.gc = (upper_bound - h.gω) .* q
		else
			throw(error("Something wrong with the individual's problem"))
		end

		h.cω = h.Φ\h.gω
	else
		vf = value(h, h.gω, s, upper_bound, ℓ, q)
	end

	# Compute expected value function
	ve = h.Emat * h.cv


	if newton
		# Get basis for the Newton step
		Φ_new = value(h, h.gω, s, upper_bound, ℓ, q, solved = true)
		jac = [h.Φ -h.β*Φ_new; -h.Emat h.Φ]
		return vf, ve, jac
	else
		return vf, ve
	end
end

function bellman_iteration!(h::Hank, R, T, ℓ, q, Π; resolve::Bool=true)
	# Compute values
	vf, ve = opt_value!(h, h.s, R, T, ℓ, q, Π, newton = false, resolve = resolve)

	# Update coefficients
	h.cv = h.Φ\vf
	h.ce = h.Φ\ve

	Void
end

function newton_iteration!(h::Hank, R, T, ℓ, q, Π)
	# Compute values
	vf, ve, jac = opt_value!(h, h.s, R, T, ℓ, q, Π)

	# Update coefficients
	cold = [h.cv; h.ce]
	c = cold - jac\([h.Φ*h.cv - vf; h.Φ*h.ce - ve])

	h.cv = c[1:h.Ns]
	h.ce = c[h.Ns+1:2*h.Ns]

	Void
end

function extend_state_space!(h::Hank, R, T, q, Π)

	gc_ext = SharedArray{Float64}(h.Nω, h.Nϵ, h.Nb, h.Nμ, h.Nσ, h.Nz, length(h.qgrid), length(h.wgrid))
	gω_ext = SharedArray{Float64}(h.Nω, h.Nϵ, h.Nb, h.Nμ, h.Nσ, h.Nz, length(h.qgrid), length(h.wgrid))
	
	Nq, Nw = length(h.qgrid), length(h.wgrid)

	@sync @parallel for jp in 0:Nq*Nw-1
		jq = 1 + jp % Nq
		qv = h.qgrid[jq]
		jw = 1 + Int( (jp - jp % Nq) / Nq)
		wv = h.wgrid[jw]

		# Re-solve for these values of w and q
		ℓ = h.θ^(-1/h.χ) * (h.s[:,2] .* wv .* (1 - h.τ)).^((1+h.χ)/h.χ)
		opt_value!(h, h.s, R, T, ℓ, qv, Π)
			
		gc_ext[:,:,:,:,:,:,jq,jw] = reshape(h.gc, h.Nω, h.Nϵ, h.Nb, h.Nμ, h.Nσ, h.Nz)
		gω_ext[:,:,:,:,:,:,jq,jw] = reshape(h.gω, h.Nω, h.Nϵ, h.Nb, h.Nμ, h.Nσ, h.Nz)
	end

	h.gc_ext = gc_ext
	h.gω_ext = gω_ext

	# Run once more at the 'normalized' values to mess down (?) the type
	ℓ = h.θ^(-1/h.χ) * (h.s[:,2] .* h.wage .* (1 - h.τ)).^((1+h.χ)/h.χ)
	opt_value!(h, h.s, R, T, ℓ, q, Π)

	Void
end

function _unpack_origvars(x, xmax, xmin)	
	f(m::Float64, cmax, cmin) = cmax - (cmax-cmin)/(1+exp(m))
	w  = f(x[1], xmax[1], xmin[1])
	Π  = f(x[2], xmax[2], xmin[2])
	qg = f(x[3], xmax[3], xmin[3])

	return w, Π, qg
end


function mkt_clearing(h::Hank, itp_ξg, itp_ξf, b, μ, σ, z, B′, Aplus, Aminus, rep, Rᵉ, Tᵉ, G, Πᵉ, itp_gc, itp_gω, itp_uc, x, xmax=x, xmin=x; get_others::Bool = false, orig_vars::Bool=true)
	F = zeros(x)
	w, Π, qg = x[1], x[2], x[3]
	if orig_vars == false
		w, Π, qg = _unpack_origvars(x, xmax, xmin)
	end

	L = (w * (1-h.τ)/h.θ * h.Ξ)^(1./h.χ)
	Y = z * L

	Tʳ = G + h.κ * rep / Π * b - qg * (B′ - (1-h.ρ)*b) - h.τ * w*L

	i = (1+h.i_star) * (Π/h.Πstar)^h.Φπ * L^h.ΦL - 1
	q = 1/(1+i)

	ψ = Y * (1 - w/z - 0.5*h.η*(Π/h.Πstar - 1)^2)

	rS = ((Aminus + rep*b*(h.κ + (1-h.ρ)*qg) + ψ) / Aplus) - 1

	valf, valg, valp, valv, sum_prob = 0., 0., 0., 0., 0.
	totcount, incount, outcount = 0, 0, 0
	for (jϵ, ϵv) in enumerate(h.ϵgrid)
		for jω = 1:length(h.ωgrid_fine)-1
			ωv  = h.ωgrid_fine[jω]
			ω1v = h.ωgrid_fine[jω+1]
			ωmv = 0.5*(ωv+ω1v)

			prob = pdf(LogNormal(μ, σ), ωmv-h.ωmin) * h.λϵ[jϵ] * (ω1v - ωv)

			Rʳ = (1+rS*(ωmv>=0))/Π
			rᵉ = (Rᵉ - 1)*(ωmv>=0)
			ω_corrected = (Rʳ*ωmv - Tʳ + Tᵉ)/((1+rᵉ)/Πᵉ)

			gω = itp_gω[ω_corrected, ϵv, b, μ, σ, z, q, w]
			uc = itp_uc[ω_corrected, ϵv, b, μ, σ, z, q, w]
			ξg = itp_ξg[ω_corrected, ϵv, b, μ, σ, z]
			ξf = itp_ξf[ω_corrected, ϵv, b, μ, σ, z]
			if ω_corrected < h.ωgrid[1] || ω_corrected > h.ωgrid[end] || q < h.qgrid[1] || q > h.qgrid[end]
				ext_gω = extrapolate(itp_gω, Interpolations.Flat()) 
				ext_uc = extrapolate(itp_uc, Interpolations.Flat()) 
				ext_ξg = extrapolate(itp_ξg, Interpolations.Flat())
				ext_ξf = extrapolate(itp_ξf, Interpolations.Flat())
				gω = ext_gω[ω_corrected, ϵv, b, μ, σ, z, q, w]

				gω < h.ωmin && ispprox(gω, h.ωmin)? gω = h.ωmin: Void
				ℓ = h.θ^(-1/h.χ) * (ϵv .* w .* (1 - h.τ)).^((1+h.χ)/h.χ)
				BC = ( Rʳ*ωmv - Tʳ + ℓ )./q
				gc = (BC - gω) .* q
				uc = gc^(-h.γ)
				
				ξg = ext_ξg[ω_corrected, ϵv, b, μ, σ, z]
				ξf = ext_ξf[ω_corrected, ϵv, b, μ, σ, z]
				
				if gω < h.ωmin
					warn("gω - h.ωmin = $(gω - h.ωmin) at ω_corr = $(ω_corrected)")
					outcount += 1
				end
			else
				if gω < h.ωmin
					warn("gω - h.ωmin = $(gω - h.ωmin) at ω_corr = $(ω_corrected)")
					incount += 1
				end
			end

			valf += prob * (gω / uc * ξf / Y)
			valg += prob * (gω / uc * ξg)
			valp += prob * (gω)
			valv += prob * (gω)^2
			sum_prob += prob

			totcount += 1
		end
	end
	if incount + outcount > 0 
		throw(error("$incount errors within bounds, $outcount outside. Total $totcount"))
	end

	F[1] = qg - valg / valp
	isnan(F[1])? warn("govt debt pricing error = $(F[1])"): Void
	Rot  = valf / valp
	""" Pensar Rotemberg + Subsidio!!! """
	Rotemberg_RHS = h.elast * (w/z - (h.elast - 1)/h.elast) + h.η * Rot
	
	if 0.25 + Rotemberg_RHS/h.η > 0
		F[2] = Π / h.Πstar - (0.5 + sqrt.(0.25 + Rotemberg_RHS/h.η))
	else
		F[2] = h.η * Π/h.Πstar * (Π/h.Πstar - 1) - Rotemberg_RHS
	end
	isnan(F[2])? warn("rotemberg error = $(F[2])"): Void

	A′ = valp / sum_prob
	varprime = valv / sum_prob - A′^2

	1 + varprime / ((A′-h.ωmin)^2) > 0 || warn("potentially negative variance at w = $w, Π = $Π, qg = $qg, q = $q")

	σ2 = log( 1 + varprime / ((A′-h.ωmin)^2) )
	μ′ = log(A′-h.ωmin) - 0.5 * σ2
	σ′ = sqrt(σ2)

	F[3] = q * A′ - qg * B′
	isnan(F[3])? warn("mf budget constraint error = $(F[3])"): Void

	if get_others
		return q, A′, μ′, σ′, rS, Tʳ, F
	else
		return F
	end
end

function wrap_find_mktclearing(h::Hank, itp_ξg, itp_ξf, b, μ, σ, z, B′, Aplus, Aminus, rep, Rᵉ, Tᵉ, G, Πᵉ, itp_gc, itp_gω, itp_uc, xguess)

	wguess, Πguess, qgguess = xguess[1], xguess[2], xguess[3]
	w, Π, qg = 1e10, 1e10, 1e10

	minw, maxw   = minimum(h.wgrid), maximum(h.wgrid)
	minΠ, maxΠ   = 0.8, 1.1
	minqg, maxqg = 0.6, h.Πstar

	function wrap_mktclear_minpack!(x::Vector, fvec=similar(x))

		xmax = [maxw, maxΠ, maxqg]
		xmin = [minw, minΠ, minqg]

		out = mkt_clearing(h, itp_ξg, itp_ξf, b, μ, σ, z, B′, Aplus, Aminus, rep, Rᵉ, Tᵉ, G, Πᵉ, itp_gc, itp_gω, itp_uc, x, xmax, xmin; orig_vars=false)

		fvec[:] = out[:]
	end

	res = fsolve(wrap_mktclear_minpack!, xguess)

	f(m::Float64, cmax, cmin) = cmax - (cmax-cmin)/(1+exp(m))
	g(c, cmax, cmin) = log( (c - cmin) / (cmax - c) )

	mw, mΠ, mqg = collect(res.:x)

	w  = f(mw,  maxw,  minw)
	Π  = f(mΠ,  maxΠ,  minΠ)
	qg = f(mqg, maxqg, minqg)

	return res.:converged, res.:f, w, Π, qg
end

function find_prices(h::Hank, itp_ξg, itp_ξf, b, μ, σ, z, B′, Aplus, Aminus, rep, Rᵉ, Tᵉ, G, Πᵉ, itp_gc, itp_gω, itp_uc, guess)

	wguess, Πguess, qgguess = guess[1], guess[2], guess[3]	

	flag, minf, w, Π, qg = wrap_find_mktclearing(h, itp_ξg, itp_ξf, b, μ, σ, z, B′, Aplus, Aminus, rep, Rᵉ, Tᵉ, G, Πᵉ, itp_gc, itp_gω, itp_uc, [wguess, Πguess, qgguess])

	curr_min = sum(minf.^2)
	minx = [w, Π, qg]
	if flag == false
		wrap_mktclearing_nlopt(x::Vector, grad::Vector) = sum(mkt_clearing(h, itp_ξg, itp_ξf, b, μ, σ, z, B′, Aplus, Aminus, rep, Rᵉ, Tᵉ, G, Πᵉ, itp_gc, itp_gω, itp_uc, x).^2)

		# alg_list = [:LN_BOBYQA; :LN_COBYLA] :GN_ISRES :GN_DIRECT_L_RAND
		alg_list = [:GN_ISRES; :LN_COBYLA]
		minw, maxw   = minimum(h.wgrid), maximum(h.wgrid)
		minΠ, maxΠ   = 0.8, 1.1
		minqg, maxqg = 0.6, h.Πstar

		ftol = 1e-4
		for jalg in 1:length(alg_list)
			opt = Opt(alg_list[jalg], 3)
			upper_bounds!(opt, [maxw, maxΠ, maxqg])
			lower_bounds!(opt, [maxw, maxΠ, maxqg])
			xtol_rel!(opt, 1e-10)
			# maxtime!(opt, 5)
			# maxeval!(opt, 50)
			min_objective!(opt, wrap_mktclearing_nlopt)
			(minf,minx,ret) = NLopt.optimize(opt, minx)

			if minf < curr_min
				curr_min = minf
				w, Π, qg = minx
			end
		end
	end

	q, A′, μ′, σ′, rS, Tʳ, minf = mkt_clearing(h, itp_ξg, itp_ξf, b, μ, σ, z, B′, Aplus, Aminus, rep, Rᵉ, Tᵉ, G, Πᵉ, itp_gc, itp_gω, itp_uc, [w, Π, qg]; get_others = true)

	return [w, Π, qg, q, μ′, σ′, rS, Tʳ], minf
end

function find_all_prices(h::Hank, itp_ξg, itp_ξf, itp_gc, itp_gω, itp_uc, repay, issuance, Rᵉ_mat, Tᵉ_mat, G_mat, Πᵉ_mat, Aplus_mat, Aminus_mat)
	results = SharedArray{Float64}(h.Nb, h.Nμ, h.Nσ, h.Nz, 8)
	minf	= SharedArray{Float64}(h.Nb, h.Nμ, h.Nσ, h.Nz, 3)

	R, T, ℓ, Rep, q, Π = _unpackstatefs(h)

	Π	= reshape(Π, h.Nω, h.Nϵ, h.Nb, h.Nμ, h.Nσ, h.Nz)[1,1,:,:,:,:]
	w 	= reshape(h.wage, h.Nω, h.Nϵ, h.Nb, h.Nμ, h.Nσ, h.Nz)[1,1,:,:,:,:]
	qg 	= reshape(h.debtprice, h.Nω, h.Nϵ, h.Nb, h.Nμ, h.Nσ, h.Nz)[1,1,:,:,:,:]

	@sync @parallel for jz in 1:length(h.zgrid)
		zv = h.zgrid[jz]
		for (jσ, σv) in enumerate(h.σgrid), (jμ, μv) in enumerate(h.μgrid), (jb, bv) in enumerate(h.bgrid)

			rep 	= repay[jb, jμ, jσ, jz]
			B′  	= issuance[jb, jμ, jσ, jz]
			Rᵉ 		= Rᵉ_mat[jb, jμ, jσ, jz]
			Tᵉ		= Tᵉ_mat[jb, jμ, jσ, jz]
			G		= G_mat[jb, jμ, jσ, jz]
			Πᵉ 		= Πᵉ_mat[jb, jμ, jσ, jz]
			Aplus	= Aplus_mat[jb, jμ, jσ, jz]
			Aminus	= Aminus_mat[jb, jμ, jσ, jz]

			guess = [w[jb,jμ,jσ,jz]; Π[jb,jμ,jσ,jz]; qg[jb,jμ,jσ,jz]]

			results[jb, jμ, jσ, jz, :], minf[jb, jμ, jσ, jz, :] = find_prices(h, itp_ξg, itp_ξf, bv, μv, σv, zv, B′, Aplus, Aminus, rep, Rᵉ, Tᵉ, G, Πᵉ, itp_gc, itp_gω, itp_uc, guess)
		end
	end
							
	return results, minf
end

function upd_P!(h::Hank, B′, μ′, σ′)
	""" Use the linear interpolation to turn laws of motion into transition matrices """
	basis 	= Basis(LinParams(h.bgrid, 0))
	B′		= max.(min.(B′, maximum(h.bgrid)), minimum(h.bgrid))
	Pb 		= BasisMatrix(basis, Expanded(), B′, 0).vals[1]

	basis 	= Basis(LinParams(h.μgrid, 0))
	μ′		= max.(min.(μ′, maximum(h.μgrid)), minimum(h.μgrid))
	Pμ 		= BasisMatrix(basis, Expanded(), μ′, 0).vals[1]

	basis 	= Basis(LinParams(h.σgrid, 0))
	σ′		= max.(min.(σ′, maximum(h.σgrid)), minimum(h.σgrid))
	Pσ 		= BasisMatrix(basis, Expanded(), σ′, 0).vals[1]

	Qz 		= kron(h.Pz, ones(h.Nb*h.Nμ*h.Nσ, 1))

	h.Ps 	= row_kron(row_kron(row_kron(Qz, Pσ), Pμ), Pb)

	for js in 1:size(h.Ps)[1]
		temp = sum(h.Ps[js,:])
		isapprox(temp,1) || warn("∑ P(s'|s) - 1 = $(@sprintf("%.3g", temp-1))")
	end

	h.Emat 	= kron(h.Ps, kron(h.Pϵ, speye(h.Nω))) * h.Φ

	Void
end

function update_state_functions!(h::Hank, upd_η)

	itp_gc = interpolate((h.ωgrid, h.ϵgrid, h.bgrid, h.μgrid, h.σgrid, h.zgrid, h.qgrid, h.wgrid), h.gc_ext, Gridded(Linear()))
	itp_uc = interpolate((h.ωgrid, h.ϵgrid, h.bgrid, h.μgrid, h.σgrid, h.zgrid, h.qgrid, h.wgrid), h.gc_ext.^(-h.γ), Gridded(Linear()))
	itp_gω  = interpolate((h.ωgrid, h.ϵgrid, h.bgrid, h.μgrid, h.σgrid, h.zgrid, h.qgrid, h.wgrid), h.gω_ext, Gridded(Linear()))
	ξg = reshape(h.ξg, h.Nω, h.Nϵ, h.Nb, h.Nμ, h.Nσ, h.Nz)
	ξf = reshape(h.ξf, h.Nω, h.Nϵ, h.Nb, h.Nμ, h.Nσ, h.Nz)
	itp_ξg = interpolate((h.ωgrid, h.ϵgrid, h.bgrid, h.μgrid, h.σgrid, h.zgrid), ξg, Gridded(Linear()))
	itp_ξf = interpolate((h.ωgrid, h.ϵgrid, h.bgrid, h.μgrid, h.σgrid, h.zgrid), ξf, Gridded(Linear()))

	repay 	 	= reshape(h.debt_repay, h.Nω, h.Nϵ, h.Nb, h.Nμ, h.Nσ, h.Nz)[1,1,:,:,:,:]
	issuance 	= reshape(h.issuance_B, h.Nω, h.Nϵ, h.Nb, h.Nμ, h.Nσ, h.Nz)[1,1,:,:,:,:]
	Rᵉ_mat 	 	= reshape(1. + h.MF_rS, h.Nω, h.Nϵ, h.Nb, h.Nμ, h.Nσ, h.Nz)[1,1,:,:,:,:]
	Tᵉ_mat	 	= reshape(h.lump_sum, h.Nω, h.Nϵ, h.Nb, h.Nμ, h.Nσ, h.Nz)[1,1,:,:,:,:]
	G_mat	 	= reshape(h.spending, h.Nω, h.Nϵ, h.Nb, h.Nμ, h.Nσ, h.Nz)[1,1,:,:,:,:]
	Πᵉ_mat	 	= reshape(h.inflation, h.Nω, h.Nϵ, h.Nb, h.Nμ, h.Nσ, h.Nz)[1,1,:,:,:,:]
	Aplus_mat	= reshape(h.Aplus, h.Nb, h.Nμ, h.Nσ, h.Nz)
	Aminus_mat	= reshape(h.Aminus, h.Nb, h.Nμ, h.Nσ, h.Nz)

	results, minf = find_all_prices(h, itp_ξg, itp_ξf, itp_gc, itp_gω, itp_uc, repay, issuance, Rᵉ_mat, Tᵉ_mat, G_mat, Πᵉ_mat, Aplus_mat, Aminus_mat)

	""" Pensar cómo suavizar el update de μ′ y σ′ """
	μ′	= reshape(results[:, :, :, :, 5], h.Nb*h.Nμ*h.Nσ*h.Nz)
	σ′	= reshape(results[:, :, :, :, 6], h.Nb*h.Nμ*h.Nσ*h.Nz)

	h.μ′ = upd_η * μ′ + (1-upd_η) * h.μ′
	h.σ′ = upd_η * σ′ + (1-upd_η) * h.σ′

	out = Array{Float64}(h.Nω, h.Nϵ, h.Nb, h.Nμ, h.Nσ, h.Nz, 8)
	for (jϵ, ϵv) in enumerate(h.ϵgrid), (jω, ωv) in enumerate(h.ωgrid)
		out[jω, jϵ, :, :, :, :, :] = results
	end
	
	h.wage 		= upd_η * reshape(out[:, :, :, :, :, :, 1], h.Ns) + (1-upd_η) * h.wage
	h.inflation = upd_η * reshape(out[:, :, :, :, :, :, 2], h.Ns) + (1-upd_η) * h.inflation
	h.debtprice = upd_η * reshape(out[:, :, :, :, :, :, 3], h.Ns) + (1-upd_η) * h.debtprice
	h.q 		= upd_η * reshape(out[:, :, :, :, :, :, 4], h.Ns) + (1-upd_η) * h.q
	h.MF_rS 	= upd_η * reshape(out[:, :, :, :, :, :, 7], h.Ns) + (1-upd_η) * h.MF_rS
	h.lump_sum 	= upd_η * reshape(out[:, :, :, :, :, :, 8], h.Ns) + (1-upd_η) * h.lump_sum

	return mean(minf[:,:,:,:,1]), mean(minf[:,:,:,:,2]), mean(minf[:,:,:,:,3])
end

function compute_ξ!(h::Hank)
	""" Tener en cuenta que los ξ's hay que calcularlos con el guess viejo """
	Y = h.Ld .* h.s[:,6]
	c = h.gc

	rep = h.debt_repay
	qg 	= h.debtprice
	Π  	= h.inflation
	P 	= kron(h.Ps, kron(h.Pϵ, speye(h.Nω)))

	ret_g = c.^(-h.γ) .* rep .* (h.κ + (1-h.ρ).*qg) ./ Π
	ret_f = c.^(-h.γ) .* Y .* Π/h.Πstar .* (Π/h.Πstar - 1)
	for js in 1:h.Ns
		h.ξg[js] = h.β * dot(P[js,:], ret_g)
		h.ξf[js] = h.β * dot(P[js,:], ret_f)
	end

	Void
end

end # @everywhere

function vfi!(h::Hank; tol::Float64=1e-2, verbose::Bool=true, maxiter::Int64=5000, bellman_iter::Int64=maxiter)
	print_save("\nSolving household problem: ")
	time_init = time()
	t_old = time_init
	iter = 0
	iter_cycle = 0
	dist, dist_s = 10., 10.
	upd_tol = 0.1

	B′ = reshape(h.issuance_B, h.Nω, h.Nϵ, h.Nb, h.Nμ, h.Nσ, h.Nz)[1,1,:,:,:,:]
	B′ = reshape(B′, h.Nb*h.Nμ*h.Nσ*h.Nz)

	upd_P!(h, B′, h.μ′, h.σ′)
	upd_η = 0.1
	dist_statefuncs = [dist]

	while dist > tol && iter < maxiter
		t_old = time()
		iter += 1
		iter_cycle += 1

		R, T, ℓ, Rep, q, Π = _unpackstatefs(h)

		compute_ξ!(h)

		c_old = [h.cv; h.ce]
		if iter_cycle <= bellman_iter
			if iter > 5 && iter % 5 != 0 || iter > 20 && iter % 15 != 0
				bellman_iteration!(h, R, T, ℓ, q, Π; resolve = false)
			else
				bellman_iteration!(h, R, T, ℓ, q, Π; resolve = true)
			end
		else
			newton_iteration!(h, R, T, ℓ, q, Π)
			c_old = [h.cv; h.ce]
			bellman_iteration!(h, R, T, ℓ, q, Π)
		end
		c_new = [h.cv; h.ce]
		
		dist = norm(c_new - c_old) / norm(c_old)
		if verbose
			# plot_hh_policies(h)
			t_new = time()
			print_save("\nd(cv, cv′) = $(@sprintf("%0.3g",dist)) at ‖v‖ = $(@sprintf("%0.3g",norm(h.Φ*h.cv))) after $(time_print(t_new-t_old)) and $iter iterations ")
		end

		if dist < upd_tol
			print_save("\nExtending the state space")
			t1 = time()
			extend_state_space!(h, R, T, q, Π)
			print_save(": done in $(time_print(time()-t1))\n")
			t1 = time()
			print_save("\nUpdating functions of the state")

			err_g, err_R, err_M = update_state_functions!(h, upd_η)
			upd_P!(h, B′, h.μ′, h.σ′)
			print_save(": done in $(time_print(time()-t1)) \nAverage errors in gov debt, PC, MF = ($(@sprintf("%0.3g",mean(err_g))), $(@sprintf("%0.3g",mean(err_R))), $(@sprintf("%0.3g",mean(err_M))))")
			iter_cycle = 0

			new_R, new_T, new_ℓ, new_Rep, new_q, new_Π = _unpackstatefs(h)

			dist_R	 = norm(new_R - R) / norm(R)
			dist_T	 = norm(new_T - T) / norm(T)
			dist_ℓ	 = norm(new_ℓ - ℓ) / norm(ℓ)
			dist_Rep = norm(new_Rep - Rep) / norm(Rep)
			dist_q	 = norm(new_q - q) / norm(q)
			dist_Π	 = norm(new_Π - Π) / norm(Π)

			dist_s = 1/upd_η * max(dist_R, dist_T, dist_ℓ, dist_Rep, dist_q, dist_Π)
			print_save("\nDistance in state functions = $(@sprintf("%0.3g", dist_s)) ")
			upd_tol = update_tolerance(dist, dist_s)
			dist = max(dist, dist_s)
			plot_hh_policies(h)
			plot_state_funcs(h)
			plot_LoM(h, h.μ′, h.σ′)
			print_save("\nNew threshold = $(@sprintf("%0.3g", upd_tol)) ")
			push!(dist_statefuncs, dist_s)
			plot(1:length(dist_statefuncs), dist_statefuncs, xlabel = L"t", yscale=:log10, label="")
			savefig(pwd() * "/../Graphs/convergence_f(S).png")
			dist_s < 1 ? dist_s < 0.5? upd_η = 0.5: upd_η = 0.25: upd_η = 0.1
		end

		dist = max(dist, dist_s)

		if iter % 10 == 0
			plot_hh_policies(h)
		end

		if isnan.(dist)
			error("NaN encountered")
		end

	end

	if dist <= tol
		print_save("\nConverged in $iter iterations. ")
	else
		print_save("\nStopping at distance $(@sprintf("%0.3g",dist)). ")
	end

	if verbose
		plot_hh_policies(h)
		plot_state_funcs(h)
	end

	print_save("\nTotal time: $(time_print(time()-time_init))\n")
	Void
end

function decomp(x::Float64)

	pot = floor(log10(x))
	bas = floor(x / 10.0^pot)

	return pot, bas
end

function maketol(bas, pot)
	step_tol, min_tol = 10.0^pot, 5*10.0^pot
	if bas <= 5
		step_tol, min_tol = 5*10.0^(pot-1), 10.0^pot
	end

	return step_tol, min_tol
end

function update_tolerance(upd_tol::Float64, dist_s::Float64)

	pot, bas = decomp(upd_tol)
	pot_s, bas_s = decomp(dist_s)
	pot = max(pot, pot_s-3)

	step_tol, min_tol = maketol(bas, pot)
	upd_tol = max(upd_tol - step_tol, min_tol)

	return upd_tol
end

function plot_hh_policies(h::Hank)
	jshow = (h.s[:,3].==median(h.bgrid)) .* (h.s[:,4].==median(h.μgrid)) .* (h.s[:,5].==median(h.σgrid)) .* (h.s[:,6].==median(h.zgrid))

	leg = Array{LaTeXStrings.LaTeXString}(1, h.Nϵ)
	for jϵ in 1:h.Nϵ
		leg[jϵ] = latexstring("\\epsilon = $(round(h.ϵgrid[jϵ],2))")
	end

	vf = h.Φ * h.cv

	jshow_b, jshow_μ, jshow_σ, jshow_z, jshow_q, jshow_w = ceil(Int64, h.Nb/2), ceil(Int64, h.Nμ/2), ceil(Int64, h.Nσ/2), ceil(Int64, h.Nz/2), findfirst(h.qgrid.>=1), findfirst(h.wgrid.>=1)

	# pc = plot(h.ωgrid, reshape(h.gc[jshow],h.Nω,h.Nϵ), lw = 2, title = "Consumption", label = leg, legend = :bottomright)
	pc = plot(h.ωgrid, reshape(h.gc_ext[:,:,jshow_b,jshow_μ,jshow_σ,jshow_z,jshow_q,jshow_w],h.Nω, h.Nϵ), title = "Consumption", label = leg, legend = :bottomright)
	# pc = plot!(h.ωgrid, reshape(h.gc[jshow],h.Nω,h.Nϵ), title = "Consumption", label = "")
	pω = plot(h.ωgrid, reshape(h.gω[jshow],h.Nω,h.Nϵ), title = "Savings", label = "")
	pv = plot(h.ωgrid, reshape(vf[jshow],h.Nω,h.Nϵ), title = "Value function", label = "")

	l = @layout([a; b c])

	plot(pc, pω, pv, layout=l, lw = 1.5, xlabel = L"\omega_t", size = (540,720))
	#plot!(bg_outside = RGBA(0.99,0.99,0.99, 0.))
	# plot!(right_margin=10px, titlefont=font(11,"Palatino"), guidefont=font(8,"Palatino"), tickfont=font(7,"Palatino"), titlefont=font(12,"Palatino"))
	savefig(pwd() * "/../Graphs/hh.pdf")
	savefig(pwd() * "/../Graphs/hh.png")

	jshow_ω, jshow_ϵ = ceil(Int64, h.Nω/2), ceil(Int64, h.Nϵ/2)

	leg = Array{LaTeXStrings.LaTeXString}(1, length(h.wgrid))
	for jw in 1:length(h.wgrid)
		leg[jw] = latexstring("w = $(round(h.wgrid[jw],2))")
	end
	pc = plot(h.qgrid, h.gc_ext[jshow_ω, jshow_ϵ,jshow_b,jshow_μ,jshow_σ,jshow_z,:,:], title = "Consumption", label = leg, legend = :bottomright)
	pω = plot(h.qgrid, h.gω_ext[jshow_ω, jshow_ϵ,jshow_b,jshow_μ,jshow_σ,jshow_z,:,:], title = "Savings", label = "")
	plot(pc, pω, layout=(2,1), lw = 1.5, xlabel = L"q_t", size = (540,720))

	savefig(pwd() * "/../Graphs/hh_qw.png")
	savefig(pwd() * "/../Graphs/hh_qw.pdf")

	return Void
end

function plot_state_funcs(h::Hank)
	R, T, ℓ, Rep, q, Π = _unpackstatefs(h)

	Aplus_mat	= reshape(h.Aplus, h.Nb, h.Nμ, h.Nσ, h.Nz)
	Aminus_mat	= reshape(h.Aminus, h.Nb, h.Nμ, h.Nσ, h.Nz)

	P = kron(h.Ps, kron(h.Pϵ, speye(h.Nω)))

	EβR = P * (R ./ Π) ./ q
	EβR = h.β * reshape(EβR, h.Nω, h.Nϵ, h.Nb, h.Nμ, h.Nσ, h.Nz)[1,1,:,:,:,:]

	R	= reshape(R, h.Nω, h.Nϵ, h.Nb, h.Nμ, h.Nσ, h.Nz)[1,1,:,:,:,:]
	T	= reshape(T, h.Nω, h.Nϵ, h.Nb, h.Nμ, h.Nσ, h.Nz)[1,1,:,:,:,:]
	q	= reshape(q, h.Nω, h.Nϵ, h.Nb, h.Nμ, h.Nσ, h.Nz)[1,1,:,:,:,:]
	Π	= reshape(Π, h.Nω, h.Nϵ, h.Nb, h.Nμ, h.Nσ, h.Nz)[1,1,:,:,:,:]
	w 	= reshape(h.wage, h.Nω, h.Nϵ, h.Nb, h.Nμ, h.Nσ, h.Nz)[1,1,:,:,:,:]
	qg 	= reshape(h.debtprice, h.Nω, h.Nϵ, h.Nb, h.Nμ, h.Nσ, h.Nz)[1,1,:,:,:,:]


	Z = zeros(h.Nb, h.Nμ, h.Nσ, h.Nz)
	for (jz, zv) in enumerate(h.zgrid)
		Z[:,:,:,jz] = zv
	end	
	L = (w * (1-h.τ)/h.θ * h.Ξ).^(1/h.χ)
	Y = Z .* L

	ψ = Y .* (1 - w./Z - 0.5*h.η*(Π/h.Πstar - 1).^2)

	Π 	= (Π.^4 - 1)*100

	i = ((1./q).^4 - 1) * 100 # Annualized percent nominal rate

	j = 3

	# l = @layout([a b; c d; e f; g h i])

	pEβR= plot(h.bgrid, vec(EβR[:,j,j,j]), title=L"\beta E[R^S]", label = "")
	pqg = plot(h.bgrid, vec(qg[:,j,j,j]), title=L"q^g", label = "")
	pΠ = plot(h.bgrid, vec(Π[:,j,j,j]), title=L"Π", label = "")
	pAp = plot(h.bgrid, vec(Aplus_mat[:,j,j,j]), title=L"A^+", label = "")
	pAm = plot(h.bgrid, vec(Aminus_mat[:,j,j,j]), title=L"A^-", label = "")
	p_i = plot(h.bgrid, vec(i[:,j,j,j]), title=L"i", label = "")
	pT = plot(h.bgrid, vec(T[:,j,j,j]), title=L"T", label = "")
	pw = plot(h.bgrid, vec(w[:,j,j,j]), title=L"w", label = "")
	pψ = plot(h.bgrid, vec(ψ[:,j,j,j]), title=L"\psi", label = "")
	
	plot(pEβR, pqg, pΠ, pAp, pAm, p_i, pT, pw, pψ, xlabel = L"B_t", layout = (3,3), lw = 1.5)
	savefig(pwd() * "/../Graphs/fs_b.png")

	# l = @layout([a b; c d; e f g])
	pEβR= plot(h.μgrid, vec(EβR[j,:,j,j]), title=L"\beta E[R^S]", label = "")
	pqg = plot(h.μgrid, vec(qg[j,:,j,j]), title=L"q^g", label = "")
	pΠ = plot(h.μgrid, vec(Π[j,:,j,j]), title=L"Π", label = "")
	pAp = plot(h.μgrid, vec(Aplus_mat[j,:,j,j]), title=L"A^+", label = "")
	pAm = plot(h.μgrid, vec(Aminus_mat[j,:,j,j]), title=L"A^-", label = "")
	p_i = plot(h.μgrid, vec(i[j,:,j,j]), title=L"i", label = "")
	pT = plot(h.μgrid, vec(T[j,:,j,j]), title=L"T", label = "")
	pw = plot(h.μgrid, vec(w[j,:,j,j]), title=L"w", label = "")
	pψ = plot(h.μgrid, vec(ψ[j,:,j,j]), title=L"\psi", label = "")
	
	plot(pEβR, pqg, pΠ, pAp, pAm, p_i, pT, pw, pψ, xlabel = L"\mu_t", layout = (3,3), lw = 1.5)
	savefig(pwd() * "/../Graphs/fs_mu.png")

	# l = @layout([a b; c d; e f g])
	pEβR= plot(h.σgrid, vec(EβR[j,j,:,j]), title=L"\beta E[R^S]", label = "")
	pqg = plot(h.σgrid, vec(qg[j,j,:,j]), title=L"q^g", label = "")
	pΠ = plot(h.σgrid, vec(Π[j,j,:,j]), title=L"Π", label = "")
	pAp = plot(h.σgrid, vec(Aplus_mat[j,j,:,j]), title=L"A^+", label = "")
	pAm = plot(h.σgrid, vec(Aminus_mat[j,j,:,j]), title=L"A^-", label = "")
	p_i = plot(h.σgrid, vec(i[j,j,:,j]), title=L"i", label = "")
	pT = plot(h.σgrid, vec(T[j,j,:,j]), title=L"T", label = "")
	pw = plot(h.σgrid, vec(w[j,j,:,j]), title=L"w", label = "")
	pψ = plot(h.σgrid, vec(ψ[j,j,:,j]), title=L"\psi", label = "")
	
	plot(pEβR, pqg, pΠ, pAp, pAm, p_i, pT, pw, pψ, xlabel = L"\sigma_t", layout = (3,3), lw = 1.5)
	savefig(pwd() * "/../Graphs/fs_sigma.png")

	# l = @layout([a b; c d; e f g])
	pEβR= plot(h.zgrid, vec(EβR[j,j,j,:]), title=L"\beta E[R^S]", label = "")
	pqg = plot(h.zgrid, vec(qg[j,j,j,:]), title=L"q^g", label = "")
	pΠ = plot(h.zgrid, vec(Π[j,j,j,:]), title=L"Π", label = "")
	pAp = plot(h.zgrid, vec(Aplus_mat[j,j,j,:]), title=L"A^+", label = "")
	pAm = plot(h.zgrid, vec(Aminus_mat[j,j,j,:]), title=L"A^-", label = "")
	p_i = plot(h.zgrid, vec(i[j,j,j,:]), title=L"i", label = "")
	pT = plot(h.zgrid, vec(T[j,j,j,:]), title=L"T", label = "")
	pw = plot(h.zgrid, vec(w[j,j,j,:]), title=L"w", label = "")
	pψ = plot(h.zgrid, vec(ψ[j,j,j,:]), title=L"\psi", label = "")
	
	plot(pEβR, pqg, pΠ, pAp, pAm, p_i, pT, pw, pψ, xlabel = L"z_t", layout = (3,3), lw = 1.5)
	savefig(pwd() * "/../Graphs/fs_z.png")

	Void
end

function plot_LoM(h::Hank, μ′, σ′)

	μ′ = reshape(μ′, h.Nb, h.Nμ, h.Nσ, h.Nz)
	σ′ = reshape(σ′, h.Nb, h.Nμ, h.Nσ, h.Nz)

	j 	= ceil(Int, length(h.bgrid)/2)
	pμb = plot(h.bgrid, vec(μ′[:,j,j,j]), title=L"\mu", xlabel = L"B_t", label = "")
	pσb = plot(h.bgrid, vec(σ′[:,j,j,j]), title=L"\sigma", xlabel = L"B_t", label = "")

	j 	= ceil(Int, length(h.μgrid)/2)
	pμμ = plot(h.μgrid, vec(μ′[j,:,j,j]), title=L"\mu", xlabel = L"\mu_t", label = "")
	pσμ = plot(h.μgrid, vec(σ′[j,:,j,j]), title=L"\sigma", xlabel = L"\mu_t", label = "")

	j 	= ceil(Int, length(h.σgrid)/2)
	pμσ = plot(h.σgrid, vec(μ′[j,j,:,j]), title=L"\mu", xlabel = L"\sigma_t", label = "")
	pσσ = plot(h.σgrid, vec(σ′[j,j,:,j]), title=L"\sigma", xlabel = L"\sigma_t", label = "")

	j 	= ceil(Int, length(h.zgrid)/2)
	pμz = plot(h.zgrid, vec(μ′[j,j,j,:]), title=L"\mu", xlabel = L"z_t", label = "")
	pσz = plot(h.zgrid, vec(σ′[j,j,j,:]), title=L"\sigma", xlabel = L"z_t", label = "")
	
	plot(pμb, pσb, pμμ, pσμ, pμσ, pσσ, pμz, pσz, layout = (4,2), lw = 1.5, size = (600,800))
	savefig(pwd() * "/../Graphs/LoMs.png")
	Void
end


function print_save(s::String)
	print(s)
	output = readstring(pwd()*"/../output.txt")
	write(pwd()*"/../output.txt", output * s)

	Void
end
