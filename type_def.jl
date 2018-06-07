using QuantEcon, BasisMatrices

type Hank
	# Utility parameters
	β::Float64
	γ::Float64
	ψ::Float64
	EpsteinZin::Bool
	γw::Float64
	θL::Float64
	χ::Float64
	Ξ::Float64
	# prop_transf::Bool

	# Debt parameters
	ρ::Float64
	κ::Float64
	r_star::Float64

	# Open Economy and production parameters
	η::Float64
	ϖ::Float64
	α_T::Float64
	α_N::Float64

	# Policy functions
	ϕa::Array{Float64, 8}
	ϕb::Array{Float64, 8}
	ϕe::Array{Float64, 8}
	ϕc::Array{Float64, 8}

	ϕa_ext::Array{Float64, 9}
	ϕb_ext::Array{Float64, 9}
	ϕe_ext::Array{Float64, 9}
	ϕc_ext::Array{Float64, 9}

	vf::Array{Float64, 8}

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
	Nw::Int64
	Nζ::Int64
	Nz::Int64
	Ns::Int64
	Nω_fine::Int64

	# Transition matrices
	Pϵ::Matrix{Float64}
	Pz::Matrix{Float64}

	# Ps::Matrix{Float64}

	# Distributions
	λ::Vector{Float64}
	λϵ::Vector{Float64}

	# Default parameters
	ℏ::Float64
	θ::Float64
	Δ::Float64

	# Parameters of the a grid
	# curv::Float64
	# order::Int64
	ωmin::Float64
	ωmax::Float64

	ωgrid0::Vector{Float64}
	ωgrid::Vector{Float64}
	ϵgrid::Vector{Float64}
	bgrid::Vector{Float64}
	μgrid::Vector{Float64}
	σgrid::Vector{Float64}
	wgrid::Vector{Float64}
	ζgrid::Vector{Float64}
	zgrid::Vector{Float64}
	s::Matrix{Float64}

	Jgrid::Matrix{Int64}

	# Extra grids for prices
	pngrid::Vector{Float64}

	# Collocation objects
	basis::Basis
	bs::BasisMatrix
	Φ::SparseMatrixCSC
	# Emat::SparseMatrixCSC

	ωgrid_fine::Vector{Float64}
	snodes::Array{Float64, 2}

	# Forecasting rules
	μ′::Array{Float64, 3} # μ′[js, jzp, 1] if reentry (or no default today)
	σ′::Array{Float64, 3} # σ′[js, jzp, 2] if default today and no reentry
	w′::Vector{Float64}

	# Functions of the state
	# A⁺::Vector{Float64}
	# A⁻::Vector{Float64}
	repay::Vector{Float64}
	τ::Float64
	T::Vector{Float64}
	issuance::Vector{Float64}
	def_thres::Vector{Float64}
	output::Vector{Float64}
	profits::Vector{Float64}
	spending::Vector{Float64}
	wage::Vector{Float64}
	Ld::Vector{Float64}
	# sdf_vec::Vector{Float64}
	qʰ::Vector{Float64}
	qᵍ::Vector{Float64}
	pN::Vector{Float64}


	# Options
	upd_tol::Float64
end

type Path
	data::Matrix{Float64}
	n::Dict{Symbol,Int64}
end
function Path(; T::Int64 = 1)
	n = Dict(
		:B => 1,
		:μ => 2,
		:σ => 3,
		:w => 4,
		:ζ => 5,
		:z => 6,
		:π => 7,
		:Y => 8,
		:L => 9,
		:ψ => 10,
		:P => 11,
		:A => 12,
		:Bh => 13,
		:Bf => 14,
		:Pe => 15,
		:Wr => 16,
		:Wd => 17
		)
	data = Matrix{Float64}(T, length(n))
	return Path(data, n)
end

getfrompath(p::Path, t::Int64, sym::Symbol) = p.data[t, p.n[sym]]
getfrompath(p::Path, t::AbstractArray, sym::Symbol) = p.data[t, p.n[sym]]
getfrompath(p::Path, sym::Symbol) = p.data[:, p.n[sym]]

function fill_path!(p::Path, t::Int64, d::Dict{Symbol, Float64}=Dict(:VOID=>-Inf);
					B::Float64=-Inf,
					μ::Float64=-Inf,
					σ::Float64=-Inf,
					w::Float64=-Inf,
					ζ::Float64=-Inf,
					z::Float64=-Inf,
					π::Float64=-Inf,
					Y::Float64=-Inf,
					L::Float64=-Inf,
					ψ::Float64=-Inf,
					P::Float64=-Inf,
					A::Float64=-Inf,
					Bh::Float64=-Inf,
					Bf::Float64=-Inf,
					Wr::Float64=-Inf,
					Wd::Float64=-Inf)
	0 < t <= size(p.data, 1) || throw("t out of bounds")
	if d != Dict(:VOID=>-Inf)
		for (jd, dv) in enumerate(d)
			sym, val = dv[1], dv[2]
			p.data[t, p.n[sym]] = val
		end
	else
		B != -Inf? p.data[t, p.n[:B]] = B: Void
		μ != -Inf? p.data[t, p.n[:μ]] = μ: Void
		σ != -Inf? p.data[t, p.n[:σ]] = σ: Void
		w != -Inf? p.data[t, p.n[:w]] = w: Void
		ζ != -Inf? p.data[t, p.n[:ζ]] = ζ: Void
		z != -Inf? p.data[t, p.n[:z]] = z: Void
		π != -Inf? p.data[t, p.n[:π]] = π: Void
		Y != -Inf? p.data[t, p.n[:Y]] = Y: Void
		L != -Inf? p.data[t, p.n[:L]] = L: Void
		ψ != -Inf? p.data[t, p.n[:ψ]] = ψ: Void
		P != -Inf? p.data[t, p.n[:P]] = P: Void
		A != -Inf? p.data[t, p.n[:A]] = A: Void
		Bh != -Inf? p.data[t, p.n[:Bh]] = Bh: Void
		Bf != -Inf? p.data[t, p.n[:Bf]] = Bf: Void
	end
	Void
end

function trim_path!(p::Path, T_burnin::Int64)
	p.data = p.data[T_burnin+1:end, :]
	Void
end
function series(p::Path, sym::Symbol)
	T = size(p.data, 1)

	y = zeros(T)
	for jt in 1:T
		y[jt] = getfrompath(p, jt, sym)
	end

	return y
end
