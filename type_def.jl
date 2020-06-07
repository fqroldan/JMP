using QuantEcon, Distributions

mutable struct SOEdef{Ktot, Kshocks}
	pars::Dict{Symbol, Float64}
	opt::Dict{Symbol, Bool}
	gr::Dict{Symbol, Vector{Float64}}
	prob::Dict{Symbol, Matrix{Float64}}

	ϕ::Dict{Symbol, Array{Float64, Ktot}}
	v::Dict{Symbol, Array{Float64, Ktot}}

	eq::Dict{Symbol, Vector{Float64}}
	gov::Dict{Symbol, Vector{Float64}}
	LoM::Dict{Symbol, Array{Vector{Float64}, Kshocks}}
end
abstract type AbstractPath
end
mutable struct Path{T} <: AbstractPath
	data::Dict{Symbol, Vector{Float64}}
end
function Path(; T::Int64 = 1)
	data = Dict( key => Vector{Float64}(undef, T) for key in [:B,:μ,:σ,:w,:ζ,:z,:π,:Y,:L,:ψ,:P,:A,:Bh,:Bf,:Pe,:Wr,:Wd,:qg,:G,:mean,:var,:CoY,:C,:CoYd,:T,:NX,:ξ,:p25,:p90,:avgω])
	return Path{T}(data)
end

periods(pv::Vector{T}) where T <: AbstractPath = sum([periods(pp) for pp in pv])
periods(p::Path{T}) where T = T

function check_periods(p::Path, t::Int64)
	0 < t <= periods(p) || throw("t out of bounds")
	nothing
end

getfrompath(p::Path, t::Int64, sym::Symbol) = p.data[sym][t]
getfrompath(p::Path, t::AbstractArray, sym::Symbol) = [p.data[sym][tv] for tv in t]
getfrompath(p::Path, t::Int) = Dict(key => p.data[key][t] for key in keys(p.data))
getfrompath(p::Path, sym::Symbol) = p.data[sym]
series(p::Path, sym::Symbol) = getfrompath(p,sym)

function fill_path!(p::Path, t::Int64, d::Dict=Dict())
	check_periods(p,t)
	missing_keys = 0
	for (key, val) in d
		if haskey(p.data, key)
			p.data[key][t] = val
		else
			missing_keys += 1
		end
	end
	
	if missing_keys > 0
		print_save("WARNING: $missing_keys missing keys")
	end
	nothing
end

function trim_path(p::Path{T}, t0::Int64) where T
	check_periods(p,t0)
	
	return Path{T-t0}(Dict(key => val[t0+1:end] for (key, val) in p.data))
end

function tauchen_fun(ny::Int64, ρ::Float64, σe::Float64; m=3, mu=0.0)
	σy = σe/sqrt((1-ρ^2))
	λ1 = -m*σy; λn = m*σy
	λgrid = range(λ1,λn,length=ny) .+ mu
	ww = λgrid[2] - λgrid[1]

	distrib = Normal(0,1)

	Π = Array{Float64}(undef, ny,ny)

	for ii=1:ny
		Π[ii,1] = cdf(distrib,(λ1+ww/2-(1-ρ)*mu-ρ*λgrid[ii])/σe) # For j=1
		Π[ii,end] = 1.0-cdf(distrib,(λn-ww/2-(1-ρ)*mu-ρ*λgrid[ii])/σe) # For j=ny
		for jj=2:ny-1
			Π[ii,jj] = cdf(distrib,(λgrid[jj]+ww/2-(1-ρ)*mu-ρ*λgrid[ii])/σe) - cdf(distrib,(λgrid[jj]-ww/2-(1-ρ)*mu-ρ*λgrid[ii])/σe)
		end
	end
	return λgrid,Π
end
function quarterlize_AR1(ρ, σ)
	ρ4 = ρ^0.25
	σ4 = sqrt(  σ^2 / ( 1 + ρ4^2 + ρ4^4 + ρ4^6 )  )
	return ρ4, σ4
end

function SOEdef(;
	β = (1.0/1.09)^0.25,
	IES = 1.0,
	RRA = 10,
	τ = 0.2,
	r_star = 1.04^0.25 - 1.0,
	tax = 0.02,
	ωmax = 20,
	wbar = 0.89,
	curv = .4,
	income_process = "Mendoza-D'Erasmo",
	EpsteinZin = true,
	order = 3,
	Nω_fine = 1000,
	Nω = 7,
	Nϵ = 5,
	Nμ = 4,
	Nσ = 4,
	Nb = 7,
	Nξ = 2,
	Nz = 7,
	ρz = 0.97,
	σz = 0.003,
	ρξ = 0.95,
	σξ = 0.0025,
	ℏ = 0.45,
	Δ = 0.1,
	θ = .04167,
	Np = 5,
	upd_tol = 5e-3,
	nodef = false,
	noΔ = false,
	rep_agent = false,
	nob = false
	)

	ψ = IES
	γ = 0.
	if EpsteinZin == true
		γ = RRA
	end

	if noΔ
		Δ = 0.0
	end

	σmin, σmax = 0.01, 1.25
	if rep_agent
		σϵ = 0.0001
		Nϵ = 2
		Nσ = 2
		σmin, σmax = 0.01, 0.02
	end

	# Debt parameters
	ρ = 0.05 # Target average maturity of 7 years: ~0.05 at quarterly freq
	κ = ρ + r_star

	## Prepare discretized processes
	# Aggregate risk
	zgrid, Pz = tauchen_fun(Nz, ρz, σz, m=1.5)

	meanξ = tax
	ξgrid, Pξ = tauchen_fun(Nξ, ρξ, σξ, m=0.5, mu=meanξ)

	ρϵ, σϵ = 0., 0.
	if income_process == "Floden-Lindé"
		ρϵ = 0.9136		# Floden-Lindé for US
		σϵ = 0.0426		# Floden-Lindé for US
	elseif income_process == "Mendoza-D'Erasmo"
		ρϵ = 0.85		# Mendoza-D'Erasmo for Spain
		σϵ = 0.2498		# Mendoza-D'Erasmo for Spain
	else
		throw(error("Must specify an income process"))
	end
	ρϵ, σϵ = quarterlize_AR1(ρϵ, σϵ)

	ϵ_chain = tauchen(Nϵ, ρϵ, σϵ, 0, 1)

	Pϵ = ϵ_chain.p
	ϵgrid = ϵ_chain.state_values

	pngrid = range(0.5, 1.1, length=Np)
	ζgrid = 0:1
	Nζ = length(ζgrid)

	α_T = 0.67
	α_N = 0.67

	ϑ = 0.88 # Straight from Anzoategui

	μ_anzo = 0.74 # Taken straight from Anzoategui, from Stockman and Tesar (1995)
	ω_anzo = 0.8  # Taken from Anzoategui, targets SS output share of nontradables at 88%

	η = μ_anzo
	ϖ = ω_anzo^(1.0/μ_anzo)

	# Grids for endogenous aggregate states
	Bmax  = 6.5
	Bbar  = Bmax * 0.5
	bgrid = range(0.0, Bmax, length=Nb)
	μgrid = range(-2.5, 0.75, length=Nμ)
	σgrid = range(σmin, σmax, length=Nσ)

	# Prepare grid for cash in hand.
	ωmin	= -0.5
	ωgrid0	= range(0.0, (ωmax-ωmin)^curv, length=Nω).^(1/curv)
	ωgrid0	= ωgrid0 .+ ωmin
	ωgrid 	= ωgrid0

	ωgrid_fine	= range(0.0, (ωmax-ωmin)^curv, length=Nω_fine).^(1/curv)
	ωgrid_fine	= ωgrid_fine .+ ωmin


	Jagg = gridmake(1:Nb, 1:Nμ, 1:Nσ, 1:Nξ, 1:Nζ, 1:Nz)

	μ = Array{Float64}(undef, Nb, Nμ, Nσ, Nξ, Nζ, Nz)
	for (jμ, μv) in enumerate(μgrid)
		μ[:,jμ,:,:,:,:,:] .= μv + 0.5 * ( mean(μgrid) - μv )
	end
	μ = reshape(μ, Nb*Nμ*Nσ*Nξ*Nζ*Nz)
	σ = Array{Float64}(undef, Nb, Nμ, Nσ, Nξ, Nζ, Nz)
	for (jσ, σv) in enumerate(σgrid)
		σ[:,:,jσ,:,:,:,:] .= σv + 0.5 * ( mean(σgrid) - σv )
	end
	σ = reshape(σ, Nb*Nμ*Nσ*Nξ*Nζ*Nz)
	μ′ = Array{Vector{Float64}}(undef, Nb*Nμ*Nσ*Nξ*Nζ*Nz, Nξ, Nz)
	σ′ = Array{Vector{Float64}}(undef, Nb*Nμ*Nσ*Nξ*Nζ*Nz, Nξ, Nz)
	for js = 1:size(Jagg,1), jξ = 1:Nξ, jz = 1:Nz
		μ′[js, jξ, jz] = [μ[js] for jj in 1:2]
		σ′[js, jξ, jz] = [σ[js] for jj in 1:2]
	end
	w′ = Array{Float64}(undef, Nb, Nμ, Nσ, Nξ, Nζ, Nz)
	for (jw, wv) in enumerate(ξgrid)
		w′[:,:,:,jw,:,:] .= wv
	end
	w′ = reshape(w′, Nb*Nμ*Nσ*Nξ*Nζ*Nz)

	Ld = ones(Nb*Nμ*Nσ*Nξ*Nζ*Nz)
	T  = ones(Nb*Nμ*Nσ*Nξ*Nζ*Nz) * 0.05
	qʰ = ones(Nb*Nμ*Nσ*Nξ*Nζ*Nz) / (1.0+r_star)
	qᵍ = ones(Nb*Nμ*Nσ*Nξ*Nζ*Nz)
	spread = zeros(Nb*Nμ*Nσ*Nξ*Nζ*Nz)

	pN 		  = Array{Float64}(undef, Nb, Nμ, Nσ, Nξ, Nζ, Nz)
	wage 	  = Array{Float64}(undef, Nb, Nμ, Nσ, Nξ, Nζ, Nz)
	spending  =	Array{Float64}(undef, Nb, Nμ, Nσ, Nξ, Nζ, Nz)
	issuance  =	Array{Float64}(undef, Nb, Nμ, Nσ, Nξ, Nζ, Nz)
	output	  = Array{Float64}(undef, Nb, Nμ, Nσ, Nξ, Nζ, Nz)
	repay 	  = Array{Float64}(undef, Nb, Nμ, Nσ, Nξ, Nζ, Nz, Nξ, Nz)
	for (jz, zv) in enumerate(zgrid)
		pN[:,:,:,:,:,jz] .= mean(pngrid) - 0.1 * zv
		output[:,:,:,:,:,jz] .= exp(zv)
		spending[:,:,:,:,:,jz] .= 0.1 - 0.25 * zv
		for (jb, bv) in enumerate(bgrid)
			issuance[jb,:,:,:,1,jz] .= bv - 0.25 * zv + 0.1 * (Bbar-bv)
			issuance[jb,:,:,:,2,jz] .= bv
		end
		for (jζ, ζv) in enumerate(ζgrid)
			def = (ζv != 1.0)
			wage[:,:,:,:,jζ,jz] .= max(exp(zv) * (1.0 - Δ * def), wbar)
		end
		repay[:,:,:,:,:,:,:,jz] .= 1.0 - (zv <= zgrid[1])
	end
	pN	 		= reshape(pN, 	 	 Nb*Nμ*Nσ*Nξ*Nζ*Nz)
	repay	 	= reshape(repay, 	 Nb*Nμ*Nσ*Nξ*Nζ*Nz*Nz*Nξ)
	if nodef
		repay = ones(repay)
	end
	wage	 	= reshape(wage, 	 Nb*Nμ*Nσ*Nξ*Nζ*Nz)
	spending 	= reshape(spending,	 Nb*Nμ*Nσ*Nξ*Nζ*Nz)
	issuance 	= min.(max.(reshape(issuance,  Nb*Nμ*Nσ*Nξ*Nζ*Nz), minimum(bgrid)), maximum(bgrid))
	output 		= reshape(output, Nb*Nμ*Nσ*Nξ*Nζ*Nz)
	profits 	= output - wage .* Ld
	welfare   	= zeros(Nb*Nμ*Nσ*Nξ*Nζ*Nz)
	C 			= zeros(Nb*Nμ*Nσ*Nξ*Nζ*Nz)

	pars = Dict(:β=>β, :γ=>γ, :ψ=>ψ, :wbar=>wbar, :ρ=>ρ, :κ=>κ, :r_star=>r_star, :τ=>τ, :η=>η, :ϖ=>ϖ, :α_T=>α_T, :α_N=>α_N, :ϑ=>ϑ, :ρϵ=>ρϵ, :σϵ=>σϵ, :ρξ=>ρξ, :σξ=>σξ, :ρz=>ρz, :σz=>σz, :ℏ=>ℏ, :θ=>θ, :Δ=>Δ, :ωmin=>ωmin, :ωmax=>ωmax, :meanξ=>meanξ)
	opt = Dict(:EpsteinZin=>EpsteinZin, :nodef=>nodef, :noΔ=>noΔ, :rep_agent=>rep_agent, :nob=>nob)
	gr = Dict(:ω=>ωgrid, :ϵ=>ϵgrid, :b=>bgrid, :μ=>μgrid, :σ=>σgrid, :ξ=>ξgrid, :ζ=>ζgrid, :z=>zgrid, :pN=>pngrid, :ωf => ωgrid_fine)
	prob = Dict(:ϵ => Pϵ, :ξ => Pξ, :z => Pz)

	ϕ = Dict(sym => zeros(Nω, Nϵ, Nb, Nμ, Nσ, Nξ, Nζ, Nz) for sym in [:a, :b, :c, :s, :θ])
	ϕ[:c] .+= 1
	v = Dict(sym => ones(Nω, Nϵ, Nb, Nμ, Nσ, Nξ, Nζ, Nz) * 1e-2 for sym in [:v, :w])

	eq = Dict(:Ld=>Ld, :T=>T,:qʰ=>qʰ,:qᵍ=>qᵍ,:spread=>spread,:pN=>pN,:wage=>wage,:spending=>spending,:issuance=>issuance,:output=>output,:welfare=>welfare,:C=>C)
	gov = Dict(:repay => repay)
	LoM = Dict(:μ => μ′, :σ => σ′)

	Ktot = length(size(ϕ[:c]))
	Kshocks = length(size(μ′))

	return SOEdef{Ktot, Kshocks}(pars, opt, gr, prob, ϕ, v, eq, gov, LoM)
end

function update_probs!(sd::SOEdef)
	update_prob_z!(sd)
	update_prob_ξ!(sd)
	nothing
end

function update_prob_z!(sd::SOEdef, ρ::Real=sd.pars[:ρz], σ::Real=sd.pars[:σz])
	Ns = N(sd, :z)
	
	zgrid, Pz = tauchen_fun(Ns, ρ, σ, m=1.5)

	sd.gr[:z] = zgrid
	sd.prob[:z] = Pz

	nothing
end

function update_prob_ξ!(sd::SOEdef, ρ::Real=sd.pars[:ρξ], σ::Real=sd.pars[:σξ], meanξ::Real=sd.pars[:meanξ])
	Ns = N(sd, :ξ)

	ξgrid, Pξ = tauchen_fun(Ns, ρ, σ, m=0.5)
	ξgrid = ξgrid .+ meanξ

	sd.gr[:ξ] = ξgrid
	sd.prob[:ξ] = Pξ

	nothing
end


ergodic_ϵ(sd::SOEdef) = first(stationary_distributions(MarkovChain(sd.prob[:ϵ], sd.gr[:ϵ])))
N(sd::SOEdef, sym) = length(sd.gr[sym])

agg_grid(sd::SOEdef) = gridmake(1:N(sd,:b), 1:N(sd,:μ), 1:N(sd,:σ), 1:N(sd,:ξ), 1:N(sd,:ζ), 1:N(sd,:z))

price_index(sd::SOEdef, pN) = price_index(sd.pars, pN)
price_index(p::Dict{Symbol,Float64}, pN) = (p[:ϖ] * pN.^(1.0-p[:η]) .+ (1.0-p[:ϖ])).^(1.0/(1.0-p[:η]))

function make_logN(meanX, varX)
	""" Takes mean and variance and returns μ and σ parameters for logNormal dist"""
	Eσ2 = 1 + varX / ( meanX^2 )

	if Eσ2 > 1
		σ2 = log( Eσ2 )
	else
		# print_save("\n1 + vω / Eω² = $(Eσ2)")
		σ2 = 1e-6
	end

	μ = log(meanX) - 0.5 * σ2
	σ = sqrt(σ2)
	return μ, σ
end

function unmake_logN(μ, σ)
	""" Takes parameters and returns mean and variance """

	m = exp.(μ + 0.5*σ.^2)
	v = exp.(σ.^2 .- 1) .* exp.(2*μ + σ.^2)

	return m, v
end
