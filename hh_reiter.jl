@everywhere begin

using QuantEcon, BasisMatrices, Interpolations, Optim, NLopt, MINPACK, Plots, LaTeXStrings, Distributions
pyplot()

type Hank
	# Utility parameters
	β::Float64
	γ::Float64
	γw::Float64

	# Debt parameters
	ρ::Float64
	κ::Float64
	r_star::Float64

	# Policy functions
	ϕa::Array{Float64, 7}
	ϕb::Array{Float64, 7}
	ϕc::Array{Float64, 7}

	ϕa_ext::Array{Float64, 10}
	ϕb_ext::Array{Float64, 10}
	ϕc_ext::Array{Float64, 10}

	vf::Array{Float64, 7}

	# Exogenous states
	ρϵ::Float64
	σϵ::Float64
	ρz::Float64
	σz::Float64

	# Grid points
	Nω::Int64
	Nϵ::Int64
	Nb::Int64
	Nμ::Int64
	Nσ::Int64
	Nz::Int64
	Nw::Int64
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

	# Parameters of the a grid
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
	wgrid::Vector{Float64}
	s::Matrix{Float64}

	Jgrid::Matrix{Int64}
	Sgrid::Matrix{Float64}

	# Extra grids for prices
	qʰgrid::Vector{Float64}
	qᵍgrid::Vector{Float64}
	wngrid::Vector{Float64}

	# Collocation objects
	basis::Basis
	bs::BasisMatrix
	Φ::SparseMatrixCSC
	# Emat::SparseMatrixCSC

	ωgrid_fine::Vector{Float64}
	snodes::Array{Float64, 2}

	# Forecasting rules
	μ′::Vector{Float64}
	σ′::Vector{Float64}
	w′::Vector{Float64}

	# Functions of the state
	A⁺::Vector{Float64}
	A⁻::Vector{Float64}
	repay::Vector{Float64}
	τ::Float64
	T::Vector{Float64}
	issuance::Vector{Float64}
	spending::Vector{Float64}
	wage::Vector{Float64}
	Ld::Vector{Float64}
	# sdf_vec::Vector{Float64}
	qʰ::Vector{Float64}
	qᵍ::Vector{Float64}

	# Mutual fund quantities
	ξg::Vector{Float64}
	ξf::Vector{Float64}
	ξp::Vector{Float64}
	sdf::String
end


function Hank(;	β = (1/1.04)^(1/4),
				γ = 2.,
				γw = 0.99,
				τ = 0.3,
				r_star = 0.02,
				ωmax = 7.5,
				curv = .4,
				order = 3,
				Nω = 8,
				Nω_fine = 1000,
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
				thr_def = -0.01,
				Nq = 7,
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

	wgrid = zgrid
	Nw = Nz

	λϵ = (Pϵ^100)[1,:]

	# Grids for endogenous aggreϕate states
	bgrid = collect(linspace(0.1, 0.4, Nb))
	μgrid = collect(linspace(0.0, 0.3, Nμ))
	σgrid = collect(linspace(1.0, 2.0, Nσ))

	Φπ = 2.0
	ΦL = 0.0

	η  = 250.0 
	elast = 6.0

	# Transitions
	Ps = Array{Float64}(Nb*Nμ*Nσ*Nz*Nw, Nb*Nμ*Nσ*Nz*Nw)

	# Prepare grid for cash in hand.
	""" Make sure that the lowest a point affords positive c at the worst prices """
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
				  LinParams(zgrid, 0),
				  LinParams(wgrid, 0))
	s, (ωgrid, ϵgrid, bgrid, μgrid, σgrid, zgrid, wgrid) = nodes(basis)
	Nω, Ns = size(ωgrid, 1), size(s, 1)

	agg_basis = Basis(LinParams(1:Nb, 0),
					  LinParams(1:Nμ, 0),
					  LinParams(1:Nσ, 0),
					  LinParams(1:Nz, 0),
					  LinParams(1:Nw, 0))
	Jgrid, _ = nodes(agg_basis)

	agg_basis = Basis(LinParams(bgrid, 0),
					  LinParams(μgrid, 0),
					  LinParams(σgrid, 0),
					  LinParams(zgrid, 0),
					  LinParams(wgrid, 0))
	Sgrid, _ = nodes(agg_basis)

	# Compute the basis matrix and expectations matrix
	bs = BasisMatrix(basis, Direct(), s, [0 0 0 0 0 0 0])
	Φ = convert(Expanded, bs).vals[1]

	ϕa = zeros(Nω, Nϵ, Nb, Nμ, Nσ, Nz, Nw)
	ϕb = zeros(Nω, Nϵ, Nb, Nμ, Nσ, Nz, Nw)
	ϕc = zeros(Nω, Nϵ, Nb, Nμ, Nσ, Nz, Nw)

	vf = ones(Nω, Nϵ, Nb, Nμ, Nσ, Nz, Nw)

	λ = ones(Nω_fine*Nϵ)
	λ = λ/sum(λ)
	
	qmin, qmax = (1.05)^(-0.25), (0.90)^(-0.25)
	qʰgrid = collect(linspace(qmin, qmax, Nq))
	qᵍgrid = collect(linspace(qmin, qmax, Nq))
	wngrid = wgrid

	ϕa_ext = zeros(Nω, Nϵ, Nb, Nμ, Nσ, Nz, Nw, Nq, Nq, Nw)
	ϕb_ext = zeros(Nω, Nϵ, Nb, Nμ, Nσ, Nz, Nw, Nq, Nq, Nw)
	ϕc_ext = zeros(Nω, Nϵ, Nb, Nμ, Nσ, Nz, Nw, Nq, Nq, Nw)

	μ′ = Array{Float64}(Nb, Nμ, Nσ, Nz, Nw)
	for (jμ, μv) in enumerate(μgrid)
		μ′[:,jμ,:,:,:] = μv + ( mean(μgrid) - μv )/2
	end
	μ′ = reshape(μ′, Nb*Nμ*Nσ*Nz*Nw)
	σ′ = Array{Float64}(Nb, Nμ, Nσ, Nz, Nw)
	for (jσ, σv) in enumerate(σgrid)
		σ′[:,:,jσ,:,:] = σv + ( mean(σgrid) - σv )/2
	end
	σ′ = reshape(σ′, Nb*Nμ*Nσ*Nz*Nw)
	w′ = Array{Float64}(Nb, Nμ, Nσ, Nz, Nw)
	for (jw, wv) in enumerate(wgrid)
		w′[:,:,:,:,jw] = wv
	end
	w′ = reshape(w′, Nb*Nμ*Nσ*Nz*Nw)

	# Debt parameters
	ρ = 0.05 # Target average maturity of 7 years: ~0.05 at quarterly freq
	κ = ρ + r_star
	
	# State functions
	Ld = ones(Nb*Nμ*Nσ*Nz*Nw)
	T  = ones(Nb*Nμ*Nσ*Nz*Nw) * 0.05
	qʰ = ones(Nb*Nμ*Nσ*Nz*Nw) * (1.0+r_star)^(-1)
	qᵍ = ones(Nb*Nμ*Nσ*Nz*Nw)

	repay = 	Array{Float64}(Nb, Nμ, Nσ, Nz, Nw)
	wage = 		Array{Float64}(Nb, Nμ, Nσ, Nz, Nw)
	spending = 	Array{Float64}(Nb, Nμ, Nσ, Nz, Nw)
	issuance = 	Array{Float64}(Nb, Nμ, Nσ, Nz, Nw)
	for (jz, zv) in enumerate(zgrid)
		repay[:,:,:,jz,:]	= 1.0 - ℏ * (log(zv) .< thr_def)
		wage[:,:,:,jz,:]	= zv
		spending[:,:,:,jz,:]= 0.2 - 0.05 * log.(zv)
		for (jb, bv) in enumerate(bgrid)
			issuance[jb,:,:,jz,:] = bv + 0.1 * log.(zv)
		end
	end
	repay	 = reshape(repay, 	 Nb*Nμ*Nσ*Nz*Nw)
	wage	 = reshape(wage, 	 Nb*Nμ*Nσ*Nz*Nw)
	spending = reshape(spending, Nb*Nμ*Nσ*Nz*Nw)
	issuance = reshape(issuance, Nb*Nμ*Nσ*Nz*Nw)

	
	ξg = zeros(Ns,)
	ξf = zeros(Ns,)
	ξp = zeros(Ns,)

	function compute_grosspositions(μ,σ)
		val⁺, val⁻, sum_prob = 0.,0.,0.
		for (jϵ, ϵv) in enumerate(ϵgrid)
			for jω = 1:length(ωgrid_fine)-1
				ωv  = ωgrid_fine[jω]
				ω1v = ωgrid_fine[jω+1]
				ωmv = 0.5*(ωv+ω1v)

				prob = pdf(LogNormal(μ, σ), ωmv-ωmin) * λϵ[jϵ] * (ω1v - ωv)

				ωmv > 0? val⁺ += prob * ωmv: val⁻ += prob * abs(ωmv)
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
	
	return Hank(β, γ, γw, ρ, κ, r_star, ϕa, ϕb, ϕc, ϕa_ext, ϕb_ext, ϕc_ext, vf, ρϵ, σϵ, ρz, σz, Nω, Nϵ, Nb, Nμ, Nσ, Nz, Nw, Ns, Nω_fine, Pϵ, Pz, Ps, λ, λϵ, ℏ, thr_def, curv, order, ωmin, ωmax, ωgrid0, ωgrid, ϵgrid, bgrid, μgrid, σgrid, zgrid, wgrid, s, Jgrid, Sgrid, qʰgrid, qᵍgrid, wngrid, basis, bs, Φ, ωgrid_fine, snodes, μ′, σ′, w′, A⁺, A⁻, repay, τ, T, issuance, spending, wage, Ld, qʰ, qᵍ, ξg, ξf, ξp, sdf)
end

function _unpackstatefs(h::Hank)

	w = h.Ld .* h.wage .* (1-h.τ)
	R = h.repay .* (h.κ + (1-h.ρ)*h.qᵍ)

	qʰ_mat = reshape(h.qʰ, 	h.Nb, h.Nμ, h.Nσ, h.Nz, h.Nw)
	qᵍ_mat = reshape(h.qᵍ, 	h.Nb, h.Nμ, h.Nσ, h.Nz, h.Nw)
	w_mat  = reshape(w, 	h.Nb, h.Nμ, h.Nσ, h.Nz, h.Nw)
	T_mat  = reshape(h.T, 	h.Nb, h.Nμ, h.Nσ, h.Nz, h.Nw)
	R_mat  = reshape(R, 	h.Nb, h.Nμ, h.Nσ, h.Nz, h.Nw)

	return qʰ_mat, qᵍ_mat, w_mat, T_mat, R_mat
end

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

function get_c(RHS::Float64, ϕa::Float64, ϕb::Float64, qʰ::Float64, qᵍ::Float64)
	return RHS - ϕa * qʰ - ϕb * qᵍ
end
function get_c(RHS::Array{Float64}, ϕa::Array{Float64}, ϕb::Array{Float64}, qʰ::Float64, qᵍ::Float64)
	return RHS - ϕa * qʰ - ϕb * qᵍ
end

function value(h::Hank, ϕa, ϕb, itp_vf, jϵ, js, RHS, qʰ, qᵍ, R_mat)

	C = get_c(RHS, ϕa, ϕb, qʰ, qᵍ)

	u = utility(h, C)

	# Basis matrix for continuation values
	jsp, Ev = 0, 0.
	for (jwp, wpv) in enumerate(h.wgrid), (jzp, zpv) in enumerate(h.zgrid), (jσp, σpv) in enumerate(h.σgrid), (jμp, μpv) in enumerate(h.μgrid), (jbp, bpv) in enumerate(h.bgrid)
		jsp += 1
		for (jϵp, ϵpv) in enumerate(h.ϵgrid)	
			ωpv = ϕa + ϕb * R_mat[jbp, jμp, jσp, jzp, jwp]
			Ev += itp_vf[ωpv, ϵpv, bpv, μpv, σpv, zpv, wpv] * h.Ps[js, jsp] * h.Pϵ[jϵ, jϵp]
		end
	end

	# Compute value
	vf = u + h.β * Ev
	return vf
end

function opt_value!(h::Hank, qʰ_mat, qᵍ_mat, w_mat, T_mat, R_mat, itp_vf; resolve::Bool = true)
	
	vf = SharedArray{Float64}(size(h.vf))
	ϕa = SharedArray{Float64}(size(h.ϕa))
	ϕb = SharedArray{Float64}(size(h.ϕb))
	ϕc = SharedArray{Float64}(size(h.ϕc))
	@sync @parallel for js in 1:size(h.Sgrid,1)
		jb = h.Jgrid[js, 1]
		jμ = h.Jgrid[js, 2]
		jσ = h.Jgrid[js, 3]
		jz = h.Jgrid[js, 4]
		jw = h.Jgrid[js, 5]
	# js = 0
	# for (jw, wv) in enumerate(h.wgrid), (jz, zv) in enumerate(h.zgrid), (jσ, σv) in enumerate(h.σgrid), (jμ, μv) in enumerate(h.μgrid), (jb, bv) in enumerate(h.bgrid)
		# js += 1
		qʰv = qʰ_mat[jb, jμ, jσ, jz, jw]
		qᵍv = qᵍ_mat[jb, jμ, jσ, jz, jw]
		wv  = w_mat[jb, jμ, jσ, jz, jw]
		Tv  = T_mat[jb, jμ, jσ, jz, jw]

		for (jϵ, ϵv) in enumerate(h.ϵgrid), (jω, ωv) in enumerate(h.ωgrid)

			RHS = ωv + wv * ϵv - Tv
			ag, bg = ϕa[jω, jϵ, jb, jμ, jσ, jz, jw], ϕb[jω, jϵ, jb, jμ, jσ, jz, jw]
			
			if qʰv * ag + qᵍv * bg > RHS
				ag, bg = RHS/qʰv - 1e-2, 0.0
			end
			
			if resolve
				guess = [ag, bg]

				wrap_value(x::Vector, grad::Vector) = value(h, x[1], x[2], itp_vf, jϵ, js, RHS, qʰv, qᵍv, R_mat)

				# println(length(guess))
				opt = Opt(:LD_LBFGS, length(guess))

				upper_bounds!(opt, [RHS/qʰv, max.(RHS/qᵍv, 0.0)])
				lower_bounds!(opt, [h.ωmin, 0.0])
				max_objective!(opt, wrap_value)
				(fmax, xmax, ret) = NLopt.optimize(opt, guess)
				ret == :SUCCESS || println(ret)

				cmax = get_c(RHS, xmax[1], xmax[2], qʰv, qᵍv)
				vf[jω, jϵ, jb, jμ, jσ, jz, jw] = fmax
				ϕa[jω, jϵ, jb, jμ, jσ, jz, jw] = xmax[1]
				ϕb[jω, jϵ, jb, jμ, jσ, jz, jw] = xmax[2]
				ϕc[jω, jϵ, jb, jμ, jσ, jz, jw] = cmax
				cmax < 0? println("c = $cmax"): Void
			else
				ap = h.ϕa[jω, jϵ, jb, jμ, jσ, jz, jw]
				bp = h.ϕb[jω, jϵ, jb, jμ, jσ, jz, jw]
				vf = value(h, ap, bp, itp_vf, jϵ, js, RHS, qʰ, qᵍ, R_mat)
			end
		end
	end

	return vf, ϕa, ϕb, ϕc
end

function bellman_iteration!(h::Hank, qʰ_mat, qᵍ_mat, w_mat, T_mat, R_mat; resolve::Bool=true)
	# Interpolate the value function
	itp_vf = interpolate((h.ωgrid, h.ϵgrid, h.bgrid, h.μgrid, h.σgrid, h.zgrid, h.wgrid), h.vf, Gridded(Linear()))

	# Compute values
	vf, ϕa, ϕb, ϕc = opt_value!(h, qʰ_mat, qᵍ_mat, w_mat, T_mat, R_mat, itp_vf, resolve = resolve)

	# Store results in the type
	h.ϕa = ϕa
	h.ϕb = ϕb
	h.ϕc = ϕc
	h.vf = vf

	Void
end


function extend_state_space!(h::Hank, R, T, q, Π)

	ϕc_ext = SharedArray{Float64}(size(h.ϕc_ext))
	ϕa_ext = SharedArray{Float64}(size(h.ϕa_ext))
	
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
		_, _, ϕc, ϕa = opt_value!(h, h.s, R, T, ℓ, qˢv, qᵇv, Π)
			
		ϕc_ext[:,:,:,:,:,:,jqˢ,jqᵇ,jw] = reshape(ϕc, h.Na, h.Nϵ, h.Nb, h.Nμ, h.Nσ, h.Nz)
		ϕa_ext[:,:,:,:,:,:,jqˢ,jqᵇ,jw] = reshape(ϕa, h.Na, h.Nϵ, h.Nb, h.Nμ, h.Nσ, h.Nz)
	end

	h.ϕc_ext = ϕc_ext
	h.ϕa_ext = ϕa_ext

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


function mkt_clearing(h::Hank, itp_ξg, itp_ξf, itp_ξp, b, μ, σ, z, B′, A⁺, A⁻, rep, Rᵉ, Tᵉ, G, Πᵉ, itp_ϕa, x, xmax=x, xmin=x; get_others::Bool = false, orig_vars::Bool=true)
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

			itp_obj_ϕa = itp_ϕa
			if a_corrected < h.agrid[1] || a_corrected > h.agrid[end] 
				itp_obj_ϕa = extrapolate(itp_ϕa, Interpolations.Linear())
			end
			if qᵇ < h.qᵇgrid[1] || qᵇ > h.qᵇgrid[end]
				itp_obj_ϕa = extrapolate(itp_ϕa, Interpolations.Flat())
			end
			ϕa = itp_obj_ϕa[a_corrected, ϵv, b, μ, σ, z, qˢ, qᵇ, w]
			ϕa < h.amin && isapprox(ϕa, h.amin)? ϕa = h.amin: Void

			ξg = itp_ξg[ϕa, ϵv, b, μ, σ, z]
			ξf = itp_ξf[ϕa, ϵv, b, μ, σ, z]
			ξp = itp_ξp[ϕa, ϵv, b, μ, σ, z]
			ℓ = h.θ^(-1/h.χ) * (ϵv .* w .* (1 - h.τ)).^((1+h.χ)/h.χ)
			BC = ( Rʳ*amv - Tʳ + ℓ )
			ϕa > 0? q = qˢ: q = qᵇ
			ϕc = BC - ϕa * q
			uc = ϕc^(-h.γ)

			if h.sdf == "risk_neutral"
				uc = 1.0
			end

			if ϕa > 0
				valf += prob * (ϕa / uc * ξf / Y)
				valg += prob * (ϕa / uc * ξg)
				valp += prob * (ϕa / uc * ξp)
				val⁺ += prob * ϕa
			else
				val⁻ += prob * abs(ϕa)
			end
			valnorm += prob * (ϕa)
			valv += prob * (ϕa)^2
			
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

	1 + varprime / ((A′-h.amin)^2) > 0 || warn("potentially neϕative variance at w = $w, Π = $Π, qᵍ = $qᵍ, q = $q")

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

function wrap_find_mktclearing(h::Hank, itp_ξg, itp_ξf, itp_ξp, b, μ, σ, z, B′, A⁺, A⁻, rep, Rᵉ, Tᵉ, G, Πᵉ, itp_ϕa, xguess, xmax, xmin)

	# wguess, Πguess, qᵍguess, qˢguess = collect(xguess)
	# w, Π, qᵍ, qˢ = collect(xguess)

	function wrap_mktclear_minpack!(x::Vector, fvec=similar(x))

		out = mkt_clearing(h, itp_ξg, itp_ξf, itp_ξp, b, μ, σ, z, B′, A⁺, A⁻, rep, Rᵉ, Tᵉ, G, Πᵉ, itp_ϕa, x, xmax, xmin; orig_vars=false)

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

function find_prices(h::Hank, itp_ξg, itp_ξf, itp_ξp, b, μ, σ, z, B′, A⁺, A⁻, rep, Rᵉ, Tᵉ, G, Πᵉ, itp_ϕa, guess, xmax, xmin)

	flag, minf, curr_xmin = wrap_find_mktclearing(h, itp_ξg, itp_ξf, itp_ξp, b, μ, σ, z, B′, A⁺, A⁻, rep, Rᵉ, Tᵉ, G, Πᵉ, itp_ϕa, guess, xmax, xmin)

	curr_min = sum(minf.^2)
	minx = copy(curr_xmin)
	# if flag == false
	# 	wrap_mktclearing_nlopt(x::Vector, grad::Vector) = sum(mkt_clearing(h, itp_ξg, itp_ξf, itp_ξp, b, μ, σ, z, B′, A⁺, A⁻, rep, Rᵉ, Tᵉ, G, Πᵉ, itp_ϕa, x).^2)

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

	A′, μ′, σ′, rS, Tʳ, minf = mkt_clearing(h, itp_ξg, itp_ξf, itp_ξp, b, μ, σ, z, B′, A⁺, A⁻, rep, Rᵉ, Tᵉ, G, Πᵉ, itp_ϕa, curr_xmin; get_others = true)

	w, Π, qᵍ, qˢ, qᵇ = curr_xmin

	return [w, Π, qᵍ, qˢ, qᵇ, μ′, σ′, rS, Tʳ], minf
end

function find_all_prices(h::Hank, itp_ξg, itp_ξf, itp_ξp, itp_ϕa, repay, issuance, Rᵉ_mat, Tᵉ_mat, G_mat, Πᵉ_mat, A⁺_mat, A⁻_mat)
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

		results[jb, jμ, jσ, jz, :], minf[jb, jμ, jσ, jz, :] = find_prices(h, itp_ξg, itp_ξf, itp_ξp, bv, μv, σv, zv, B′, A⁺, A⁻, rep, Rᵉ, Tᵉ, G, Πᵉ, itp_ϕa, guess, xmax, xmin)
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

	basis 	= Basis(LinParams(h.wgrid, 0))
	w′		= h.w′
	Pw 		= BasisMatrix(basis, Expanded(), w′, 0).vals[1]

	Q 		= kron(h.Pz, ones(h.Nw*h.Nσ*h.Nμ*h.Nb, 1))

	h.Ps 	= row_kron(row_kron(row_kron(row_kron(Pw, Q), Pσ), Pμ), Pb)

	err_count = 0
	for js in 1:size(h.Ps)[1]
		temp = sum(h.Ps[js,:])
		if ~isapprox(temp,1)
			warn("∑ P(s'|s) - 1 = $(@sprintf("%.3g", temp-1))")
			err_count += 1
		end
	end
	err_count = sum(err_count)

	err_count == 0 || warn("$err_count errors encountered")
	# h.Emat 	= kron(h.Ps, kron(h.Pϵ, speye(h.Nω))) * h.Φ

	Void
end

function update_state_functions!(h::Hank, upd_η)

	itp_ϕa  = interpolate((h.agrid, h.ϵgrid, h.bgrid, h.μgrid, h.σgrid, h.zgrid, h.qˢgrid, h.qᵇgrid, h.wgrid), h.ϕa_ext, Gridded(Linear()))
	ξg 		= reshape(h.ξg, h.Na, h.Nϵ, h.Nb, h.Nμ, h.Nσ, h.Nz)
	ξf 		= reshape(h.ξf, h.Na, h.Nϵ, h.Nb, h.Nμ, h.Nσ, h.Nz)
	ξp 		= reshape(h.ξp, h.Na, h.Nϵ, h.Nb, h.Nμ, h.Nσ, h.Nz)
	itp_ξg 	= interpolate((h.agrid, h.ϵgrid, h.bgrid, h.μgrid, h.σgrid, h.zgrid), ξg, Gridded(Linear()))
	itp_ξf 	= interpolate((h.agrid, h.ϵgrid, h.bgrid, h.μgrid, h.σgrid, h.zgrid), ξf, Gridded(Linear()))
	itp_ξp 	= interpolate((h.agrid, h.ϵgrid, h.bgrid, h.μgrid, h.σgrid, h.zgrid), ξp, Gridded(Linear()))

	repay 	 = reshape(h.repay, h.Na, h.Nϵ, h.Nb, h.Nμ, h.Nσ, h.Nz)[1,1,:,:,:,:]
	issuance = reshape(h.issuance, h.Na, h.Nϵ, h.Nb, h.Nμ, h.Nσ, h.Nz)[1,1,:,:,:,:]
	Rᵉ_mat 	 = reshape(1. + h.MF_rS, h.Na, h.Nϵ, h.Nb, h.Nμ, h.Nσ, h.Nz)[1,1,:,:,:,:]
	Tᵉ_mat	 = reshape(h.T, h.Na, h.Nϵ, h.Nb, h.Nμ, h.Nσ, h.Nz)[1,1,:,:,:,:]
	G_mat	 = reshape(h.spending, h.Na, h.Nϵ, h.Nb, h.Nμ, h.Nσ, h.Nz)[1,1,:,:,:,:]
	Πᵉ_mat	 = reshape(h.inflation, h.Na, h.Nϵ, h.Nb, h.Nμ, h.Nσ, h.Nz)[1,1,:,:,:,:]
	A⁺_mat	 = reshape(h.A⁺, h.Nb, h.Nμ, h.Nσ, h.Nz)
	A⁻_mat	 = reshape(h.A⁻, h.Nb, h.Nμ, h.Nσ, h.Nz)

	results, minf = find_all_prices(h, itp_ξg, itp_ξf, itp_ξp, itp_ϕa, repay, issuance, Rᵉ_mat, Tᵉ_mat, G_mat, Πᵉ_mat, A⁺_mat, A⁻_mat)

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
	h.T 	= upd_η * reshape(out[:, :, :, :, :, :, 9], h.Ns) + (1-upd_η) * h.T

	meanf = zeros(size(minf)[end])
	for jf in 1:size(minf)[end]
		meanf[jf] = mean(minf[:,:,:,:,jf])
	end

	return meanf
end

function compute_ξ!(h::Hank)
	rep = h.repay
	qᵍ 	= h.qᵍ
	Π  	= h.inflation
	P 	= kron(h.Ps, kron(h.Pϵ, speye(h.Na)))
	w   = h.wage
	L   = (w * (1-h.τ)/h.θ * h.Ξ).^(1/h.χ)
	Z 	= h.s[:,6]

	Y 	= L .* Z
	uc 	= h.ϕc.^(-h.γ)

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
	upd_tol = 0.01

	B′ = h.issuance

	upd_P!(h, B′, h.μ′, h.σ′)
	upd_η = 0.1
	dist_statefuncs = [dist]

	while dist > tol && iter < maxiter
		t_old = time()
		iter += 1
		iter_cycle += 1

		qʰ_mat, qᵍ_mat, w_mat, T_mat, R_mat = _unpackstatefs(h)

		v_old = copy(h.vf)
		if iter <= 5 || iter % 7 == 0 || iter_cycle == 1
			bellman_iteration!(h, qʰ_mat, qᵍ_mat, w_mat, T_mat, R_mat; resolve=true)
		else
			bellman_iteration!(h, qʰ_mat, qᵍ_mat, w_mat, T_mat, R_mat; resolve=false)
		end
		v_new = copy(h.vf)
		
		dist = sqrt.(sum( (v_new - v_old).^2 )) / sqrt.(sum(v_old.^2))
		norm_v = sqrt.(sum(v_old.^2))
		if verbose
			# plot_hh_policies(h)
			t_new = time()
			print_save("\nd(cv, cv′) = $(@sprintf("%0.3g",dist)) at ‖v‖ = $(@sprintf("%0.3g",norm_v)) after $(time_print(t_new-t_old)) and $iter iterations ")
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

	# pc = plot(h.agrid, reshape(h.ϕc[jshow],h.Na,h.Nϵ), lw = 2, title = "Consumption", label = leg, legend = :bottomright)
	pc = plot(h.agrid, reshape(h.ϕc[jshow],h.Na,h.Nϵ), title = "Consumption", label = leg, legend = :bottomright)
	# pc = plot!(h.agrid, reshape(h.ϕc[jshow],h.Na,h.Nϵ), title = "Consumption", label = "")
	pa = plot(h.agrid, reshape(h.ϕa[jshow],h.Na,h.Nϵ), title = "Savings", label = "")
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
	pc = plot(h.qˢgrid, h.ϕc_ext[jshow_a, jshow_ϵ,jshow_b,jshow_μ,jshow_σ,jshow_z,:,jshow_qᵇ,:], title = "Consumption", label = leg, legend = :bottomright)
	pa = plot(h.qˢgrid, h.ϕa_ext[jshow_a, jshow_ϵ,jshow_b,jshow_μ,jshow_σ,jshow_z,:,jshow_qᵇ,:], title = "Savings", label = "")
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
