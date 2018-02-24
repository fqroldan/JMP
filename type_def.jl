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
	ϕc::Array{Float64, 8}

	ϕa_ext::Array{Float64, 9}
	ϕb_ext::Array{Float64, 9}
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