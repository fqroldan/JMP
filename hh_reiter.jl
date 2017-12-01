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
	ga::Vector{Float64}

	gc_ext::Array{Float64, 9}
	ga_ext::Array{Float64, 9}

	# Coefficients of value and policy functions
	cc::Vector{Float64}
	ca::Vector{Float64}
	cv::Vector{Float64}
	ce::Vector{Float64}

	# Exogenous state
	ρz::Float64
	σz::Float64

	# Grid points
	Na::Int64
	Nϵ::Int64
	Nb::Int64
	Nμ::Int64
	Nσ::Int64
	Nz::Int64
	Ns::Int64
	Na_fine::Int64

	# Transition matrices
	Pϵ::Matrix{Float64}
	Pz::Matrix{Float64}

	Ps::Matrix{Float64}

	# Distributions
	λ::Vector{Float64}
	λϵ::Vector{Float64}
	ℏ::Float64
	thr_def::Float64

	# Parameters of the a grid
	curv::Float64
	order::Int64
	amin::Float64
	amax::Float64

	agrid0::Vector{Float64}
	agrid::Vector{Float64}
	ϵgrid::Vector{Float64}
	bgrid::Vector{Float64}
	μgrid::Vector{Float64}
	σgrid::Vector{Float64}
	zgrid::Vector{Float64}
	s::Matrix{Float64}

	# Extra grids for prices
	qˢgrid::Vector{Float64}
	qᵇgrid::Vector{Float64}
	wgrid::Vector{Float64}

	# Collocation objects
	basis::Basis
	bs::BasisMatrix
	Φ::SparseMatrixCSC
	dΦ::SparseMatrixCSC
	Emat::SparseMatrixCSC
	Φnota::SparseMatrixCSC

	agrid_fine::Vector{Float64}
	snodes::Array{Float64, 2}

	# Forecasting rules
	μ′::Vector{Float64}
	σ′::Vector{Float64}

	# Functions of the state
	A⁺::Vector{Float64}
	A⁻::Vector{Float64}
	debt_repay::Vector{Float64}
	MF_rS::Vector{Float64}
	τ::Float64
	lump_sum::Vector{Float64}
	issuance_B::Vector{Float64}
	spending::Vector{Float64}
	wage::Vector{Float64}
	Ld::Vector{Float64}
	# sdf_vec::Vector{Float64}
	qᵍ::Vector{Float64}
	qˢ::Vector{Float64}
	qᵇ::Vector{Float64}
	inflation::Vector{Float64}
	Πstar::Float64
	i_star::Float64

	# Mutual fund quantities
	ξg::Vector{Float64}
	ξf::Vector{Float64}
	ξp::Vector{Float64}
	sdf::String
end


function Hank(;	β = (1/1.04)^(1/4),
				γ = 2.,
				θ = 1.,
				χ = 2.,
				τ = 0.3,
				amax = 7.5,
				curv = .4,
				order = 3,
				Na = 8,
				Na_fine = 1000,
				ρϵ = 0.95,
				σϵ = 0.005,
				Nϵ = 3,
				Nμ = 5,
				Nσ = 5,
				Nb = 5,
				ρz = 0.92,
				σz = 0.005,
				Nz = 5,
				ℏ = .5,
				thr_def = -0.1,
				Nq = 7,
				Nw = 3,
				sdf = "risk_neutral"		# possible values: risk_neutral, agg_C, avg
				)
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

	# Grids for endogenous aggregate states
	bgrid = collect(linspace(0.1, 0.4, Nb))
	μgrid = collect(linspace(0.0, 0.3, Nμ))
	σgrid = collect(linspace(1.0, 2.0, Nσ))

	Φπ = 2.0
	ΦL = 0.0

	η  = 250.0 
	elast = 6.0

	# Transitions
	Ps  = Array{Float64}(Nb*Nμ*Nσ*Nz, Nb*Nμ*Nσ*Nz)

	# Prepare grid for cash in hand.
	""" Make sure that the lowest a point affords positive c at the worst prices """
	amin	= -1.0
	agrid0	= linspace(0., (amax-amin)^curv, Na).^(1/curv)
	agrid0	= agrid0 + amin

	agrid_fine	= linspace(0., (amax-amin)^curv, Na_fine).^(1/curv)
	agrid_fine	= agrid_fine + amin

	snodes = [kron(ones(Nϵ,), agrid_fine) kron(ϵgrid, ones(Na_fine,))]

	# Define the basis over the state variables
	basis = Basis(SplineParams(agrid0, 0, order),
				  LinParams(ϵgrid, 0),
				  LinParams(bgrid, 0),
				  LinParams(μgrid, 0),
				  LinParams(σgrid, 0),
				  LinParams(zgrid, 0))
	s, (agrid, ϵgrid, bgrid, μgrid, σgrid, zgrid) = nodes(basis)
	Na, Ns = size(agrid, 1), size(s, 1)

	# Compute the basis matrix and expectations matrix
	bs = BasisMatrix(basis, Direct(), s, [0 0 0 0 0 0])
	Φ = convert(Expanded, bs).vals[1]

	dΦ = BasisMatrix(basis, Expanded(), s, [1 0 0 0 0 0]).vals[1]

	# Save the parts of the interpolation that are not 'a'
	Φϵ = bs.vals[2]
	Φb = bs.vals[3]
	Φμ = bs.vals[4]
	Φσ = bs.vals[5]
	Φz = bs.vals[6]
	Φnota = row_kron(Φz, row_kron(Φσ, row_kron(Φμ, row_kron(Φb, Φϵ))))

	Emat = kron(Ps, kron(Pϵ, speye(Na))) * Φ

	cc = ones(Ns,)
	ca = ones(Ns,)

	gc = ones(Ns,)
	ga = ones(Ns,)

	cv = ones(Ns,)
	ce = ones(Ns,)

	λ = ones(Na_fine*Nϵ)
	λ = λ/sum(λ)
	
	qmin, qmax = (1.05)^(-0.25), (0.90)^(-0.25)
	qˢgrid = collect(linspace(qmin, qmax, Nq))
	qᵇgrid = collect(linspace(qmin, qmax, Nq))
	wgrid = collect(linspace(.9, 1.1, Nw)) * (elast-1)/elast

	gc_ext = zeros(Na, Nϵ, Nb, Nμ, Nσ, Nz, Nq, Nq, Nw)
	ga_ext = zeros(Na, Nϵ, Nb, Nμ, Nσ, Nz, Nq, Nq, Nw)

	haircut = ℏ * (log.(s[:, 6]) .< thr_def)
	debt_repay = 1-haircut

	lump_sum = ones(Ns,) * 0.05

	issuance_B = s[:, 3] + 0.1 * log.(s[:, 6])
	μ′ = Array{Float64}(Nb, Nμ, Nσ, Nz)
	for (jμ, μv) in enumerate(μgrid)
		μ′[:,jμ,:,:] = μv + ( mean(μgrid) - μv )/2
	end
	μ′ = reshape(μ′, Nb*Nμ*Nσ*Nz)
	σ′ = Array{Float64}(Nb, Nμ, Nσ, Nz)
	for (jσ, σv) in enumerate(σgrid)
		σ′[:,:,jσ,:] = σv + ( mean(σgrid) - σv )/2
	end
	σ′ = reshape(σ′, Nb*Nμ*Nσ*Nz)

	spending = 0.2 - 0.05 * log.(s[:, 6])

	wage = s[:, 6] * (elast-1)/elast
	θ = (elast-1)/elast * (1-τ) * Ξ
	Ld 	 = 1 * ones(Ns,)

	Πstar = 1.02^(0.25)
	i_star = (1 / 100 + Πstar^4)^(1/4) - 1
	inflation = Πstar * ones(Ns,)
	qˢ = (1/(1+i_star)) * ones(Ns,)
	qᵇ = (1/(1+i_star)) * ones(Ns,)
	qᵍ = 1.0 * ones(Ns,)

	R = 1+i_star
	r = ones(Ns,) * (R-1.)

	MF_rS = r

	# Debt parameters
	ρ = 0.05 # Target average maturity of 7 years: ~0.05 at quarterly freq
	κ = ρ + i_star

	ξg = zeros(Ns,)
	ξf = zeros(Ns,)
	ξp = zeros(Ns,)

	function compute_grosspositions(μ,σ)
		val⁺, val⁻, sum_prob = 0.,0.,0.
		for (jϵ, ϵv) in enumerate(ϵgrid)
			for ja = 1:length(agrid_fine)-1
				av  = agrid_fine[ja]
				a1v = agrid_fine[ja+1]
				amv = 0.5*(av+a1v)

				prob = pdf(LogNormal(μ, σ), amv-amin) * λϵ[jϵ] * (a1v - av)

				amv > 0? val⁺ += prob * amv: val⁻ += prob * abs(amv)
				sum_prob += prob
			end
		end
		a⁺ = val⁺ / sum_prob
		a⁻ = val⁻ / sum_prob
		
		return a⁺, a⁻
	end

	A⁺_mat, A⁻_mat = zeros(Nμ,Nσ), zeros(Nμ,Nσ)
	for (jσ, σv) in enumerate(σgrid), (jμ, μv) in enumerate(μgrid)
		A⁺_mat[jμ, jσ], A⁻_mat[jμ, jσ] = compute_grosspositions(μv, σv)
	end

	js = 0
	A⁺, A⁻ = zeros(Nb*Nμ*Nσ*Nz), zeros(Nb*Nμ*Nσ*Nz)
	for (jz, zv) in enumerate(zgrid), (jσ,σv) in enumerate(σgrid), (jμ,μv) in enumerate(μgrid), (jb,bv) in enumerate(bgrid)
		js += 1
		A⁺[js], A⁻[js] = A⁺_mat[jμ,jσ], A⁻_mat[jμ,jσ]
	end
	
	return Hank(β, γ, θ, χ, ρϵ, σϵ, Ξ, ρ, κ, Φπ, ΦL, η, elast, gc, ga, gc_ext, ga_ext, cc, ca, cv, ce, ρz, σz, Na, Nϵ, Nb, Nμ, Nσ, Nz, Ns, Na_fine, Pϵ, Pz, Ps, λ, λϵ, ℏ, thr_def, curv, order, amin, amax, agrid0, agrid, ϵgrid, bgrid, μgrid, σgrid, zgrid, s, qˢgrid, qᵇgrid, wgrid, basis, bs, Φ, dΦ, Emat, Φnota, agrid_fine, snodes, μ′, σ′, A⁺, A⁻, debt_repay, MF_rS, τ, lump_sum, issuance_B, spending, wage, Ld, qᵍ, qˢ, qᵇ, inflation, Πstar, i_star, ξg, ξf, ξp, sdf)
end

function _unpackstatefs(h::Hank)
	R = 1 + h.MF_rS
	T = h.lump_sum
	ℓ = h.θ^(-1/h.χ) * (h.s[:,2] .* h.wage .* (1 - h.τ)).^((1+h.χ)/h.χ)
	Rep = h.debt_repay
	qˢ = h.qˢ
	qᵇ = h.qᵇ
	Π = h.inflation

	return R, T, ℓ, Rep, qˢ, qᵇ, Π
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

function value(h::Hank, gc::Vector{Float64}, ga::Vector{Float64}, s, RHS::Vector{Float64}, ℓ, qˢ, qᵇ; solved::Bool=false)

	q = qˢ .* (ga.>0) + qᵇ .* (ga.<=0)

	gc[:] = RHS - ga .* q

	Ut = utility(h, gc - ℓ/(1+h.χ) )

	Φap	= BasisMatrix(h.basis[1], Direct(), ga, 0).vals[1]
	Φ_new = row_kron(h.Φnota, Φap)

	if solved
		h.cc = h.Φ\gc
		return Φ_new
	end

	# Compute value
	vf = Ut + h.β * Φ_new*h.ce 
	return vf
end

function opt_value!(h::Hank, s::Matrix{Float64}, R, T, ℓ, qˢ, qᵇ, Π; newton::Bool=false, resolve::Bool = true)
	rS = (R - 1) .* (s[:,1] .>= 0)
	# Decide savings
	lower_bound = ones(h.Ns,) * h.amin					# Borrowing constraint
	BC = (1+rS).* s[:,1]./Π + ℓ - T 
	q = qˢ .* (BC.>0) + qᵇ .* (BC.<=0)
	upper_bound = BC./q	# Budget constraint; c ≧ 0
	
	# if minimum(upper_bound-lower_bound) < 0
	# 	warn("Budget set empty. $(minimum(upper_bound-lower_bound))")
	# end
	
	vf = zeros(h.ga)
	gc = copy(h.gc)
	ga = copy(h.ga)
	if resolve
		ga_debt, vf_debt = golden_method((ga -> value(h, gc, ga, s, BC, ℓ, qˢ, qᵇ)), lower_bound, min.(0,upper_bound))
		ga_save, vf_save = golden_method((ga -> value(h, gc, ga, s, BC, ℓ, qˢ, qᵇ)), zeros(BC), upper_bound)
		
		ga = ga_debt + (ga_save - ga_debt) .* (vf_save .>= vf_debt) .* (BC .> 0)
		vf = vf_debt + (vf_save - vf_debt) .* (vf_save .>= vf_debt) .* (BC .> 0)

		if minimum(ga) > h.amin
			Void
		elseif isapprox(minimum(ga), h.amin)
			ga = max.(ga, h.amin)
		else
			throw(error("Something wrong with the individual's problem"))
		end
		q = qˢ .* (ga.>0) + qᵇ .* (ga.<=0)
		gc = BC - ga .* q
	else
		vf = value(h, gc, h.ga, s, BC, ℓ, qˢ, qᵇ)
	end

	# Compute expected value function
	ve = h.Emat * h.cv

	if newton
		# Get basis for the Newton step
		Φ_new = value(h, gc, ga, s, BC, ℓ, qˢ, qᵇ, solved = true)
		jac = [h.Φ -h.β*Φ_new; -h.Emat h.Φ]
		return vf, ve, jac
	else
		return vf, ve, gc, ga
	end
end

function bellman_iteration!(h::Hank, R, T, ℓ, qˢ, qᵇ, Π; resolve::Bool=true)
	# Compute values
	vf, ve, gc, ga = opt_value!(h, h.s, R, T, ℓ, qˢ, qᵇ, Π, resolve = resolve)

	h.gc = gc
	h.ga = ga

	# Update coefficients
	h.ca = h.Φ\ga
	h.cc = h.Φ\gc
	h.cv = h.Φ\vf
	h.ce = h.Φ\ve

	Void
end

function newton_iteration!(h::Hank, R, T, ℓ, qˢ, qᵇ, Π)
	# Compute values
	vf, ve, jac = opt_value!(h, h.s, R, T, ℓ, qˢ, qᵇ, Π)

	# Update coefficients
	cold = [h.cv; h.ce]
	c = cold - jac\([h.Φ*h.cv - vf; h.Φ*h.ce - ve])

	h.cv = c[1:h.Ns]
	h.ce = c[h.Ns+1:2*h.Ns]

	Void
end

function extend_state_space!(h::Hank, R, T, qˢ, qᵇ, Π)

	gc_ext = SharedArray{Float64}(size(h.gc_ext))
	ga_ext = SharedArray{Float64}(size(h.ga_ext))
	
	Nqˢ, Nqᵇ, Nw = length(h.qˢgrid), length(h.qᵇgrid), length(h.wgrid)

	qwgrid = [kron(ones(Int,Nqᵇ*Nw),1:Nqˢ) kron(kron(ones(Int,Nw),1:Nqᵇ),ones(Int,Nqˢ)) kron(1:Nw,ones(Int,Nqᵇ*Nqˢ))]
	Nqw = size(qwgrid)[1]

	@sync @parallel for jp in 1:Nqw
		jqˢ = qwgrid[jp, 1]
		jqᵇ = qwgrid[jp, 2]
		jw  = qwgrid[jp, 3]
		
		qˢv = h.qˢgrid[jqˢ]
		qᵇv = h.qᵇgrid[jqᵇ]
		wv  = h.wgrid[jw]
		# Re-solve for these values of w and q
		ℓ = h.θ^(-1/h.χ) * (h.s[:,2] .* wv .* (1 - h.τ)).^((1+h.χ)/h.χ)
		_, _, gc, ga = opt_value!(h, h.s, R, T, ℓ, qˢv, qᵇv, Π)
			
		gc_ext[:,:,:,:,:,:,jqˢ,jqᵇ,jw] = reshape(gc, h.Na, h.Nϵ, h.Nb, h.Nμ, h.Nσ, h.Nz)
		ga_ext[:,:,:,:,:,:,jqˢ,jqᵇ,jw] = reshape(ga, h.Na, h.Nϵ, h.Nb, h.Nμ, h.Nσ, h.Nz)
	end

	h.gc_ext = gc_ext
	h.ga_ext = ga_ext

	Void
end

transform_vars(m::Float64, cmax, cmin) = cmax - (cmax-cmin)/(1+exp(m))

function _unpack_origvars(x, xmax, xmin)
	y = zeros(x)
	for (jx, xv) in enumerate(x)
		y[jx] = transform_vars(xv, xmax[jx], xmin[jx])
	end

	return y
end


function mkt_clearing(h::Hank, itp_ξg, itp_ξf, itp_ξp, b, μ, σ, z, B′, A⁺, A⁻, rep, Rᵉ, Tᵉ, G, Πᵉ, itp_ga, x, xmax=x, xmin=x; get_others::Bool = false, orig_vars::Bool=true)
	F = zeros(x)
	w, Π, qᵍ, qˢ, qᵇ = collect(x)
	if orig_vars == false
		w, Π, qᵍ, qˢ, qᵇ = _unpack_origvars(x, xmax, xmin)
	end

	L = (w * (1-h.τ)/h.θ * h.Ξ)^(1./h.χ)
	Y = z * L

	Tʳ = G + h.κ * rep / Π * b - qᵍ * (B′ - rep / Π * (1-h.ρ)*b) - h.τ * w*L

	ψ = Y * (1 - w/z - 0.5*h.η*(Π/h.Πstar - 1)^2)

	rS = ((A⁻ + rep*b*(h.κ + (1-h.ρ)*qᵍ) + Π*ψ) / A⁺) - 1

	valf, valg, valnorm, valv, valp, val⁺, val⁻, sum_prob = 0., 0., 0., 0., 0., 0., 0., 0.
	for (jϵ, ϵv) in enumerate(h.ϵgrid)
		for ja = 1:length(h.agrid_fine)-1
			av  = h.agrid_fine[ja]
			a1v = h.agrid_fine[ja+1]
			amv = 0.5*(av+a1v)

			prob = pdf(LogNormal(μ, σ), amv-h.amin) * h.λϵ[jϵ] * (a1v - av)

			Rʳ = (1+rS*(amv>=0))/Π
			rᵉ = (Rᵉ - 1)*(amv>=0)
			a_corrected = (Rʳ*amv - Tʳ + Tᵉ)/((1+rᵉ)/Πᵉ)

			itp_obj_ga = itp_ga
			if a_corrected < h.agrid[1] || a_corrected > h.agrid[end] 
				itp_obj_ga = extrapolate(itp_ga, Interpolations.Linear())
			end
			if qᵇ < h.qᵇgrid[1] || qᵇ > h.qᵇgrid[end]
				itp_obj_ga = extrapolate(itp_ga, Interpolations.Flat())
			end
			ga = itp_obj_ga[a_corrected, ϵv, b, μ, σ, z, qˢ, qᵇ, w]
			ga < h.amin && isapprox(ga, h.amin)? ga = h.amin: Void

			ξg = itp_ξg[ga, ϵv, b, μ, σ, z]
			ξf = itp_ξf[ga, ϵv, b, μ, σ, z]
			ξp = itp_ξp[ga, ϵv, b, μ, σ, z]
			ℓ = h.θ^(-1/h.χ) * (ϵv .* w .* (1 - h.τ)).^((1+h.χ)/h.χ)
			BC = ( Rʳ*amv - Tʳ + ℓ )
			ga > 0? q = qˢ: q = qᵇ
			gc = BC - ga * q
			uc = gc^(-h.γ)

			if h.sdf == "risk_neutral"
				uc = 1.0
			end

			if ga > 0
				valf += prob * (ga / uc * ξf / Y)
				valg += prob * (ga / uc * ξg)
				valp += prob * (ga / uc * ξp)
				val⁺ += prob * ga
			else
				val⁻ += prob * abs(ga)
			end
			valnorm += prob * (ga)
			valv += prob * (ga)^2
			
			sum_prob += prob
		end
	end

	F[1] = qᵍ - valg / val⁺
	isnan(F[1])? warn("govt debt pricing error = $(F[1])"): Void

	Rot  = valf / val⁺
	""" Pensar Rotemberg + Subsidio!!! """
	Rotemberg_RHS = h.elast * (w/z - (h.elast - 1)/h.elast) + h.η * Rot
	
	F[2] = Π / h.Πstar - (0.5 + sqrt.(0.25 + Rotemberg_RHS/h.η))
	if 0.25 + Rotemberg_RHS/h.η > 0
	else
		F[2] = h.η * Π/h.Πstar * (Π/h.Πstar - 1) - Rotemberg_RHS
	end
	isnan(F[2])? warn("rotemberg error = $(F[2])"): Void

	A′ = valnorm / sum_prob
	varprime = valv / sum_prob - A′^2

	1 + varprime / ((A′-h.amin)^2) > 0 || warn("potentially negative variance at w = $w, Π = $Π, qᵍ = $qᵍ, q = $q")

	σ2 = log( 1 + varprime / ((A′-h.amin)^2) )
	μ′ = log(A′-h.amin) - 0.5 * σ2
	σ′ = sqrt(σ2)

	savings 	= val⁺ / sum_prob
	borrowing 	= val⁻ / sum_prob

	F[3] = qˢ * savings - (qᵇ * borrowing + qᵍ * B′)
	isnan(F[3])? warn("mf budget constraint error = $(F[3])"): Void


	F[4] = qᵇ - valp / val⁺
	isnan(F[4])? warn("Euler equation error = $(F[4])"): Void

	F[5] = 1/qᵇ - (1+h.i_star) * ((Π)/h.Πstar)^h.Φπ * L^h.ΦL

	# weights = [1;1e-2;1e-2;1]
	weights = ones(F)

	if get_others
		return A′, μ′, σ′, rS, Tʳ, F
	else
		return F .* weights
	end
end

function wrap_find_mktclearing(h::Hank, itp_ξg, itp_ξf, itp_ξp, b, μ, σ, z, B′, A⁺, A⁻, rep, Rᵉ, Tᵉ, G, Πᵉ, itp_ga, xguess, xmax, xmin)

	# wguess, Πguess, qᵍguess, qˢguess = collect(xguess)
	# w, Π, qᵍ, qˢ = collect(xguess)

	function wrap_mktclear_minpack!(x::Vector, fvec=similar(x))

		out = mkt_clearing(h, itp_ξg, itp_ξf, itp_ξp, b, μ, σ, z, B′, A⁺, A⁻, rep, Rᵉ, Tᵉ, G, Πᵉ, itp_ga, x, xmax, xmin; orig_vars=false)

		fvec[:] = out[:]
	end

	res = fsolve(wrap_mktclear_minpack!, xguess)
	if res.:converged == false
		res2 = fsolve(wrap_mktclear_minpack!, xguess, method=:lmdif)

		if res2.:converged || sum(res2.:f.^2) < sum(res.:f.^2)
			res = res2
		end
	end

	w, Π, qᵍ, qˢ, qᵇ = _unpack_origvars(res.:x, xmax, xmin)

	return res.:converged, res.:f, [w; Π; qᵍ; qˢ; qᵇ]
end

function find_prices(h::Hank, itp_ξg, itp_ξf, itp_ξp, b, μ, σ, z, B′, A⁺, A⁻, rep, Rᵉ, Tᵉ, G, Πᵉ, itp_ga, guess, xmax, xmin)

	flag, minf, curr_xmin = wrap_find_mktclearing(h, itp_ξg, itp_ξf, itp_ξp, b, μ, σ, z, B′, A⁺, A⁻, rep, Rᵉ, Tᵉ, G, Πᵉ, itp_ga, guess, xmax, xmin)

	curr_min = sum(minf.^2)
	minx = copy(curr_xmin)
	# if flag == false
	# 	wrap_mktclearing_nlopt(x::Vector, grad::Vector) = sum(mkt_clearing(h, itp_ξg, itp_ξf, itp_ξp, b, μ, σ, z, B′, A⁺, A⁻, rep, Rᵉ, Tᵉ, G, Πᵉ, itp_ga, x).^2)

	# 	# alg_list = [:LN_BOBYQA; :LN_COBYLA] :GN_ISRES :GN_DIRECT_L_RAND
	# 	alg_list = [:GN_ISRES; :LN_COBYLA]

	# 	for jalg in 1:length(alg_list)
	# 		opt = Opt(alg_list[jalg], length(guess))
	# 		upper_bounds!(opt, xmax)
	# 		lower_bounds!(opt, xmin)
	# 		xtol_rel!(opt, 1e-10)
	# 		maxtime!(opt, 2)
	# 		# maxeval!(opt, 50)
	# 		min_objective!(opt, wrap_mktclearing_nlopt)
	# 		(minf,minx,ret) = NLopt.optimize(opt, minx)

	# 		if minf < curr_min
	# 			curr_min  = minf
	# 			curr_xmin = minx
	# 		end
	# 	end
	# end

	A′, μ′, σ′, rS, Tʳ, minf = mkt_clearing(h, itp_ξg, itp_ξf, itp_ξp, b, μ, σ, z, B′, A⁺, A⁻, rep, Rᵉ, Tᵉ, G, Πᵉ, itp_ga, curr_xmin; get_others = true)

	w, Π, qᵍ, qˢ, qᵇ = curr_xmin

	return [w, Π, qᵍ, qˢ, qᵇ, μ′, σ′, rS, Tʳ], minf
end

function find_all_prices(h::Hank, itp_ξg, itp_ξf, itp_ξp, itp_ga, repay, issuance, Rᵉ_mat, Tᵉ_mat, G_mat, Πᵉ_mat, A⁺_mat, A⁻_mat)
	results = SharedArray{Float64}(h.Nb, h.Nμ, h.Nσ, h.Nz, 9)
	minf	= SharedArray{Float64}(h.Nb, h.Nμ, h.Nσ, h.Nz, 5)

	R, T, ℓ, Rep, qˢ, qᵇ, Π = _unpackstatefs(h)

	Π 	= reshape(Π, 		h.Na, h.Nϵ, h.Nb, h.Nμ, h.Nσ, h.Nz)[1,1,:,:,:,:]
	w 	= reshape(h.wage, 	h.Na, h.Nϵ, h.Nb, h.Nμ, h.Nσ, h.Nz)[1,1,:,:,:,:]
	qᵍ 	= reshape(h.qᵍ, 	h.Na, h.Nϵ, h.Nb, h.Nμ, h.Nσ, h.Nz)[1,1,:,:,:,:]
	qˢ 	= reshape(qˢ, 		h.Na, h.Nϵ, h.Nb, h.Nμ, h.Nσ, h.Nz)[1,1,:,:,:,:]
	qᵇ 	= reshape(qᵇ, 		h.Na, h.Nϵ, h.Nb, h.Nμ, h.Nσ, h.Nz)[1,1,:,:,:,:]

	exogrid = [kron(ones(Int,h.Nμ*h.Nσ*h.Nz),1:h.Nb) kron(kron(ones(Int,h.Nσ*h.Nz),1:h.Nμ),ones(Int,h.Nb)) kron(kron(ones(Int,h.Nz),1:h.Nσ),ones(Int,h.Nμ*h.Nb)) kron(1:h.Nz,ones(Int,h.Nσ*h.Nμ*h.Nb))]
	N_exo = size(exogrid)[1]

	@sync @parallel for jp in 1:N_exo
		jb = exogrid[jp, 1]
		jμ = exogrid[jp, 2]
		jσ = exogrid[jp, 3]
		jz = exogrid[jp, 4]

		bv = h.bgrid[jb]
		μv = h.μgrid[jμ]
		σv = h.σgrid[jσ]
		zv = h.zgrid[jz]

		rep = repay[jb, jμ, jσ, jz]
		B′  = issuance[jb, jμ, jσ, jz]
		Rᵉ 	= Rᵉ_mat[jb, jμ, jσ, jz]
		Tᵉ	= Tᵉ_mat[jb, jμ, jσ, jz]
		G	= G_mat[jb, jμ, jσ, jz]
		Πᵉ 	= Πᵉ_mat[jb, jμ, jσ, jz]
		A⁺	= A⁺_mat[jb, jμ, jσ, jz]
		A⁻	= A⁻_mat[jb, jμ, jσ, jz]

		guess = [w[jb,jμ,jσ,jz]; Π[jb,jμ,jσ,jz]; qᵍ[jb,jμ,jσ,jz]; qˢ[jb,jμ,jσ,jz]; qᵇ[jb,jμ,jσ,jz]]

		minw, maxw 	= minimum(h.wgrid), maximum(h.wgrid)
		minΠ, maxΠ 	= (h.Πstar-0.1)^(0.25), (h.Πstar+0.1)^(0.25)
		minqᵍ, maxqᵍ	= 0.9, h.Πstar
		minqˢ, maxqˢ	= minimum(h.qˢgrid), maximum(h.qˢgrid)
		minqᵇ, maxqᵇ	= minimum(h.qᵇgrid), maximum(h.qᵇgrid)

		xmin = [minw; minΠ; minqᵍ; minqˢ; minqᵇ]
		xmax = [maxw; maxΠ; maxqᵍ; maxqˢ; maxqᵇ]

		results[jb, jμ, jσ, jz, :], minf[jb, jμ, jσ, jz, :] = find_prices(h, itp_ξg, itp_ξf, itp_ξp, bv, μv, σv, zv, B′, A⁺, A⁻, rep, Rᵉ, Tᵉ, G, Πᵉ, itp_ga, guess, xmax, xmin)
	end
							
	return results, minf
end

function upd_P!(h::Hank, B′, μ′, σ′)
	""" Use linear interpolation to turn laws of motion into transition matrices """
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

	h.Emat 	= kron(h.Ps, kron(h.Pϵ, speye(h.Na))) * h.Φ

	Void
end

function update_state_functions!(h::Hank, upd_η)

	itp_ga  = interpolate((h.agrid, h.ϵgrid, h.bgrid, h.μgrid, h.σgrid, h.zgrid, h.qˢgrid, h.qᵇgrid, h.wgrid), h.ga_ext, Gridded(Linear()))
	ξg 		= reshape(h.ξg, h.Na, h.Nϵ, h.Nb, h.Nμ, h.Nσ, h.Nz)
	ξf 		= reshape(h.ξf, h.Na, h.Nϵ, h.Nb, h.Nμ, h.Nσ, h.Nz)
	ξp 		= reshape(h.ξp, h.Na, h.Nϵ, h.Nb, h.Nμ, h.Nσ, h.Nz)
	itp_ξg 	= interpolate((h.agrid, h.ϵgrid, h.bgrid, h.μgrid, h.σgrid, h.zgrid), ξg, Gridded(Linear()))
	itp_ξf 	= interpolate((h.agrid, h.ϵgrid, h.bgrid, h.μgrid, h.σgrid, h.zgrid), ξf, Gridded(Linear()))
	itp_ξp 	= interpolate((h.agrid, h.ϵgrid, h.bgrid, h.μgrid, h.σgrid, h.zgrid), ξp, Gridded(Linear()))

	repay 	 = reshape(h.debt_repay, h.Na, h.Nϵ, h.Nb, h.Nμ, h.Nσ, h.Nz)[1,1,:,:,:,:]
	issuance = reshape(h.issuance_B, h.Na, h.Nϵ, h.Nb, h.Nμ, h.Nσ, h.Nz)[1,1,:,:,:,:]
	Rᵉ_mat 	 = reshape(1. + h.MF_rS, h.Na, h.Nϵ, h.Nb, h.Nμ, h.Nσ, h.Nz)[1,1,:,:,:,:]
	Tᵉ_mat	 = reshape(h.lump_sum, h.Na, h.Nϵ, h.Nb, h.Nμ, h.Nσ, h.Nz)[1,1,:,:,:,:]
	G_mat	 = reshape(h.spending, h.Na, h.Nϵ, h.Nb, h.Nμ, h.Nσ, h.Nz)[1,1,:,:,:,:]
	Πᵉ_mat	 = reshape(h.inflation, h.Na, h.Nϵ, h.Nb, h.Nμ, h.Nσ, h.Nz)[1,1,:,:,:,:]
	A⁺_mat	 = reshape(h.A⁺, h.Nb, h.Nμ, h.Nσ, h.Nz)
	A⁻_mat	 = reshape(h.A⁻, h.Nb, h.Nμ, h.Nσ, h.Nz)

	results, minf = find_all_prices(h, itp_ξg, itp_ξf, itp_ξp, itp_ga, repay, issuance, Rᵉ_mat, Tᵉ_mat, G_mat, Πᵉ_mat, A⁺_mat, A⁻_mat)

	""" Pensar cómo suavizar el update de μ′ y σ′ """
	μ′	= reshape(results[:, :, :, :, 6], h.Nb*h.Nμ*h.Nσ*h.Nz)
	σ′	= reshape(results[:, :, :, :, 7], h.Nb*h.Nμ*h.Nσ*h.Nz)

	h.μ′ = upd_η * μ′ + (1-upd_η) * h.μ′
	h.σ′ = upd_η * σ′ + (1-upd_η) * h.σ′

	out = Array{Float64}(h.Na, h.Nϵ, h.Nb, h.Nμ, h.Nσ, h.Nz, size(results)[5])
	for (jϵ, ϵv) in enumerate(h.ϵgrid), (ja, av) in enumerate(h.agrid)
		out[ja, jϵ, :, :, :, :, :] = results
	end
	
	h.wage 		= upd_η * reshape(out[:, :, :, :, :, :, 1], h.Ns) + (1-upd_η) * h.wage
	h.inflation = upd_η * reshape(out[:, :, :, :, :, :, 2], h.Ns) + (1-upd_η) * h.inflation
	h.qᵍ 		= upd_η * reshape(out[:, :, :, :, :, :, 3], h.Ns) + (1-upd_η) * h.qᵍ
	h.qˢ 		= upd_η * reshape(out[:, :, :, :, :, :, 4], h.Ns) + (1-upd_η) * h.qˢ
	h.qᵇ 		= upd_η * reshape(out[:, :, :, :, :, :, 5], h.Ns) + (1-upd_η) * h.qᵇ
	h.MF_rS 	= upd_η * reshape(out[:, :, :, :, :, :, 8], h.Ns) + (1-upd_η) * h.MF_rS
	h.lump_sum 	= upd_η * reshape(out[:, :, :, :, :, :, 9], h.Ns) + (1-upd_η) * h.lump_sum

	meanf = zeros(size(minf)[end])
	for jf in 1:size(minf)[end]
		meanf[jf] = mean(minf[:,:,:,:,jf])
	end

	return meanf
end

function compute_ξ!(h::Hank)
	rep = h.debt_repay
	qᵍ 	= h.qᵍ
	Π  	= h.inflation
	P 	= kron(h.Ps, kron(h.Pϵ, speye(h.Na)))
	w   = h.wage
	L   = (w * (1-h.τ)/h.θ * h.Ξ).^(1/h.χ)
	Z 	= h.s[:,6]

	Y 	= L .* Z
	uc 	= h.gc.^(-h.γ)

	if h.sdf == "risk_neutral"
		uc = ones(uc)
	end

	ret_g = uc .* rep .* (h.κ + (1-h.ρ).*qᵍ) ./ Π
	ret_f = uc .* Y .* Π/h.Πstar .* (Π/h.Πstar - 1)
	ret_p = uc .* 1.0 ./ Π
	
	h.ξg = h.β * P * ret_g
	h.ξf = h.β * P * ret_f
	h.ξp = h.β * P * ret_p

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

	B′ = reshape(h.issuance_B, h.Na, h.Nϵ, h.Nb, h.Nμ, h.Nσ, h.Nz)[1,1,:,:,:,:]
	B′ = reshape(B′, h.Nb*h.Nμ*h.Nσ*h.Nz)

	upd_P!(h, B′, h.μ′, h.σ′)
	upd_η = 0.1
	dist_statefuncs = [dist]

	while dist > tol && iter < maxiter
		t_old = time()
		iter += 1
		iter_cycle += 1

		R, T, ℓ, Rep, qˢ, qᵇ, Π = _unpackstatefs(h)

		c_old = [h.cv; h.ce]
		if iter_cycle <= bellman_iter
			if iter <= 5 || iter % 7 == 0 || iter_cycle == 1
				bellman_iteration!(h, R, T, ℓ, qˢ, qᵇ, Π; resolve = true)
			else
				bellman_iteration!(h, R, T, ℓ, qˢ, qᵇ, Π; resolve = false)
			end
		else
			newton_iteration!(h, R, T, ℓ, qˢ, qᵇ, Π)
			c_old = [h.cv; h.ce]
			bellman_iteration!(h, R, T, ℓ, qˢ, qᵇ, Π)
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
			extend_state_space!(h, R, T, qˢ, qᵇ, Π)
			print_save(": done in $(time_print(time()-t1))\n")
			t1 = time()
			print_save("\nUpdating functions of the state")

			compute_ξ!(h)
			err_g, err_R, err_M, err_E, err_T = update_state_functions!(h, upd_η)
			upd_P!(h, B′, h.μ′, h.σ′)
			print_save(": done in $(time_print(time()-t1)) \nAverage errors in GD, PC, MF, RF, TR = ($(@sprintf("%0.3g",mean(err_g))), $(@sprintf("%0.3g",mean(err_R))), $(@sprintf("%0.3g",mean(err_M))), $(@sprintf("%0.3g",mean(err_E))), $(@sprintf("%0.3g",mean(err_T))))")
			iter_cycle = 0

			new_R, new_T, new_ℓ, new_Rep, new_qˢ, new_qᵇ, new_Π = _unpackstatefs(h)

			dist_R	 = norm(new_R - R) / norm(R)
			dist_T	 = norm(new_T - T) / norm(T)
			dist_ℓ	 = norm(new_ℓ - ℓ) / norm(ℓ)
			dist_Rep = norm(new_Rep - Rep) / norm(Rep)
			dist_qˢ	 = norm(new_qˢ - qˢ) / norm(qˢ)
			dist_qᵇ	 = norm(new_qᵇ - qᵇ) / norm(qᵇ)
			dist_Π	 = norm(new_Π - Π) / norm(Π)

			dist_s = 1/upd_η * max(dist_R, dist_T, dist_ℓ, dist_Rep, dist_qˢ, dist_qᵇ, dist_Π)
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

	jshow_b, jshow_μ, jshow_σ, jshow_z, jshow_qˢ, jshow_qᵇ, jshow_w = ceil(Int64, h.Nb/2), ceil(Int64, h.Nμ/2), ceil(Int64, h.Nσ/2), ceil(Int64, h.Nz/2), findfirst(h.qˢgrid.>=1), findfirst(h.qᵇgrid.>=1), findfirst(h.wgrid.>=1)

	# pc = plot(h.agrid, reshape(h.gc[jshow],h.Na,h.Nϵ), lw = 2, title = "Consumption", label = leg, legend = :bottomright)
	pc = plot(h.agrid, reshape(h.gc[jshow],h.Na,h.Nϵ), title = "Consumption", label = leg, legend = :bottomright)
	# pc = plot!(h.agrid, reshape(h.gc[jshow],h.Na,h.Nϵ), title = "Consumption", label = "")
	pa = plot(h.agrid, reshape(h.ga[jshow],h.Na,h.Nϵ), title = "Savings", label = "")
	pv = plot(h.agrid, reshape(vf[jshow],h.Na,h.Nϵ), title = "Value function", label = "")

	l = @layout([a; b c])

	plot(pc, pa, pv, layout=l, lw = 1.5, xlabel = L"a_t", size = (540,720))
	#plot!(bg_outside = RGBA(0.99,0.99,0.99, 0.))
	# plot!(right_margin=10px, titlefont=font(11,"Palatino"), guidefont=font(8,"Palatino"), tickfont=font(7,"Palatino"), titlefont=font(12,"Palatino"))
	savefig(pwd() * "/../Graphs/hh.pdf")
	savefig(pwd() * "/../Graphs/hh.png")

	Nq = length(h.qˢgrid)
	jshow_a, jshow_ϵ, jshow_qˢ, jshow_qᵇ = ceil(Int64, h.Na/2), ceil(Int64, h.Nϵ/2), ceil(Int64, Nq/2), ceil(Int64, Nq/2)

	leg = Array{LaTeXStrings.LaTeXString}(1, length(h.wgrid))
	for jw in 1:length(h.wgrid)
		leg[jw] = latexstring("w = $(round(h.wgrid[jw],2))")
	end
	pc = plot(h.qˢgrid, h.gc_ext[jshow_a, jshow_ϵ,jshow_b,jshow_μ,jshow_σ,jshow_z,:,jshow_qᵇ,:], title = "Consumption", label = leg, legend = :bottomright)
	pa = plot(h.qˢgrid, h.ga_ext[jshow_a, jshow_ϵ,jshow_b,jshow_μ,jshow_σ,jshow_z,:,jshow_qᵇ,:], title = "Savings", label = "")
	plot(pc, pa, layout=(2,1), lw = 1.5, xlabel = L"q^s_t", size = (540,720))

	savefig(pwd() * "/../Graphs/hh_qw.png")
	savefig(pwd() * "/../Graphs/hh_qw.pdf")

	return Void
end

function plot_state_funcs(h::Hank)
	R, T, ℓ, Rep, qˢ, qᵇ, Π = _unpackstatefs(h)

	A⁺_mat	= reshape(h.A⁺, h.Nb, h.Nμ, h.Nσ, h.Nz)
	A⁻_mat	= reshape(h.A⁻, h.Nb, h.Nμ, h.Nσ, h.Nz)

	P = kron(h.Ps, kron(h.Pϵ, speye(h.Na)))

	EβR = (P * (R ./ Π)) ./ qˢ
	EβR = h.β * reshape(EβR, h.Na, h.Nϵ, h.Nb, h.Nμ, h.Nσ, h.Nz)[1,1,:,:,:,:]

	R	= reshape(R, 		h.Na, h.Nϵ, h.Nb, h.Nμ, h.Nσ, h.Nz)[1,1,:,:,:,:]
	T	= reshape(T, 		h.Na, h.Nϵ, h.Nb, h.Nμ, h.Nσ, h.Nz)[1,1,:,:,:,:]
	qˢ	= reshape(qˢ, 		h.Na, h.Nϵ, h.Nb, h.Nμ, h.Nσ, h.Nz)[1,1,:,:,:,:]
	qᵇ	= reshape(qᵇ, 		h.Na, h.Nϵ, h.Nb, h.Nμ, h.Nσ, h.Nz)[1,1,:,:,:,:]
	Π	= reshape(Π, 		h.Na, h.Nϵ, h.Nb, h.Nμ, h.Nσ, h.Nz)[1,1,:,:,:,:]
	w 	= reshape(h.wage, 	h.Na, h.Nϵ, h.Nb, h.Nμ, h.Nσ, h.Nz)[1,1,:,:,:,:]
	qᵍ 	= reshape(h.qᵍ, 	h.Na, h.Nϵ, h.Nb, h.Nμ, h.Nσ, h.Nz)[1,1,:,:,:,:]

	Πdev = (Π)/h.Πstar

	Z = zeros(h.Nb, h.Nμ, h.Nσ, h.Nz)
	for (jz, zv) in enumerate(h.zgrid)
		Z[:,:,:,jz] = zv
	end	
	L = (w * (1-h.τ)/h.θ * h.Ξ).^(1/h.χ)
	Y = Z .* L

	ψ = Y .* (1 - w./Z - 0.5*h.η*(Π/h.Πstar - 1).^2)

	Π = (Π.^4 - 1)*100

	i = ((1./qᵇ).^4 - 1) * 100 # Annualized percent nominal rate
	iˢ = ((1./qˢ).^4 - 1) * 100
	
	j = 3

	# l = @layout([a b; c d; e f; g h i])

	pEβR = plot(h.bgrid, vec(EβR[:,j,j,j]), title=L"\beta E[R^S]", label = "")
	pqᵍ	 = plot(h.bgrid, vec(qᵍ[:,j,j,j]), 	title=L"q^g", label = "")
	pΠ	 = plot(h.bgrid, vec(Π[:,j,j,j]), 	title=L"Π", label = "")
	pAp	 = plot(h.bgrid, vec(A⁺_mat[:,j,j,j]), title=L"A^+", label = "")
	pAm	 = plot(h.bgrid, vec(A⁻_mat[:,j,j,j]), title=L"A^-", label = "")
	p_i	 = plot(h.bgrid, [vec(iˢ[:,j,j,j]) vec(i[:,j,j,j])], title=L"i", label = [L"i^s" L"i^b"])
	pT	 = plot(h.bgrid, vec(T[:,j,j,j]), 	title=L"T", label = "")
	pL	 = plot(h.bgrid, [vec(L[:,j,j,j]) vec(Πdev[:,j,j,j])], title=L"L", label = [L"L" L"\tilde{\pi}"])
	pw 	 = plot(h.bgrid, vec(w[:,j,j,j]), title=L"w", label = "")
	pψ = plot(h.bgrid, vec(ψ[:,j,j,j]), title=L"\psi", label = "")
	
	plot(pEβR, pqᵍ, pΠ, pAp, pAm, p_i, pT, pL, pψ, xlabel = L"B_t", layout = (3,3), lw = 1.5)
	savefig(pwd() * "/../Graphs/fs_b.png")

	# l = @layout([a b; c d; e f g])
	pEβR = plot(h.μgrid, vec(EβR[j,:,j,j]), title=L"\beta E[R^S]", label = "")
	pqᵍ	 = plot(h.μgrid, vec(qᵍ[j,:,j,j]), title=L"q^g", label = "")
	pΠ	 = plot(h.μgrid, vec(Π[j,:,j,j]), title=L"Π", label = "")
	pAp	 = plot(h.μgrid, vec(A⁺_mat[j,:,j,j]), title=L"A^+", label = "")
	pAm	 = plot(h.μgrid, vec(A⁻_mat[j,:,j,j]), title=L"A^-", label = "")
	p_i	 = plot(h.μgrid, [vec(iˢ[j,:,j,j]) vec(i[j,:,j,j])], title=L"i", label = [L"i^s" L"i^b"])
	pT	 = plot(h.μgrid, vec(T[j,:,j,j]), title=L"T", label = "")
	pL	 = plot(h.μgrid, [vec(L[j,:,j,j]) vec(Πdev[j,:,j,j])], title=L"L", label = [L"L" L"\tilde{\pi}"])
	pw = plot(h.μgrid, vec(w[j,:,j,j]), title=L"w", label = "")
	pψ = plot(h.μgrid, vec(ψ[j,:,j,j]), title=L"\psi", label = "")
	
	plot(pEβR, pqᵍ, pΠ, pAp, pAm, p_i, pT, pL, pψ, xlabel = L"\mu_t", layout = (3,3), lw = 1.5)
	savefig(pwd() * "/../Graphs/fs_mu.png")

	# l = @layout([a b; c d; e f g])
	pEβR = plot(h.σgrid, vec(EβR[j,j,:,j]), title=L"\beta E[R^S]", label = "")
	pqᵍ	 = plot(h.σgrid, vec(qᵍ[j,j,:,j]), title=L"q^g", label = "")
	pΠ	 = plot(h.σgrid, vec(Π[j,j,:,j]), title=L"Π", label = "")
	pAp	 = plot(h.σgrid, vec(A⁺_mat[j,j,:,j]), title=L"A^+", label = "")
	pAm	 = plot(h.σgrid, vec(A⁻_mat[j,j,:,j]), title=L"A^-", label = "")
	p_i	 = plot(h.σgrid, [vec(iˢ[j,j,:,j]) vec(i[j,j,:,j])], title=L"i", label = [L"i^s" L"i^b"])
	pT	 = plot(h.σgrid, vec(T[j,j,:,j]), title=L"T", label = "")
	pL	 = plot(h.σgrid, [vec(L[j,j,:,j]) vec(Πdev[j,j,:,j])], title=L"L", label = [L"L" L"\tilde{\pi}"])
	pw = plot(h.σgrid, vec(w[j,j,:,j]), title=L"w", label = "")
	pψ = plot(h.σgrid, vec(ψ[j,j,:,j]), title=L"\psi", label = "")
	
	plot(pEβR, pqᵍ, pΠ, pAp, pAm, p_i, pT, pL, pψ, xlabel = L"\sigma_t", layout = (3,3), lw = 1.5)
	savefig(pwd() * "/../Graphs/fs_sigma.png")

	# l = @layout([a b; c d; e f g])
	pEβR = plot(h.zgrid, vec(EβR[j,j,j,:]), title=L"\beta E[R^S]", label = "")
	pqᵍ	 = plot(h.zgrid, vec(qᵍ[j,j,j,:]), title=L"q^g", label = "")
	pΠ	 = plot(h.zgrid, vec(Π[j,j,j,:]), title=L"Π", label = "")
	pAp	 = plot(h.zgrid, vec(A⁺_mat[j,j,j,:]), title=L"A^+", label = "")
	pAm	 = plot(h.zgrid, vec(A⁻_mat[j,j,j,:]), title=L"A^-", label = "")
	p_i	 = plot(h.zgrid, [vec(iˢ[j,j,j,:]) vec(i[j,j,j,:])], title=L"i", label = [L"i^s" L"i^b"])
	pT	 = plot(h.zgrid, vec(T[j,j,j,:]), title=L"T", label = "")
	pL	 = plot(h.zgrid, [vec(L[j,j,j,:]) vec(Πdev[j,j,j,:])], title=L"L", label = [L"L" L"\tilde{\pi}"])
	pw = plot(h.zgrid, vec(w[j,j,j,:]), title=L"w", label = "")
	pψ = plot(h.zgrid, vec(ψ[j,j,j,:]), title=L"\psi", label = "")
	
	plot(pEβR, pqᵍ, pΠ, pAp, pAm, p_i, pT, pL, pψ, xlabel = L"z_t", layout = (3,3), lw = 1.5)
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
