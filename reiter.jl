using QuantEcon, BasisMatrices, Interpolations, Optim, MINPACK, LaTeXStrings, Distributions, JLD, Cubature

include("hh_pb.jl")

function Hank(;	β = (1.0/1.05)^0.25,
				IES = 1.0,
				RRA = 5.,
				γw = 1.0,#0.99,#^0.25,
				τ = 0.2,
				r_star = 1.04^0.25 - 1.0,
				tax = 0.1,
				ωmax = 20.,
				wbar = 1.05,
				curv = .4,
				income_process = "Floden-Lindé",
				EpsteinZin = true,
				order = 3,
				Nω_fine = 2500,
				Nω = 7,
				Nϵ = 5,
				Nμ = 4,
				Nσ = 4,
				Nb = 9,
				Nw = 2,
				Nz = 13,
				ρz = 0.9,
				σz = 0.025,
				ℏ = 0.4,
				Δ = 0.1,
				θ = .125,
				Np = 5,
				upd_tol = 5e-3,
				nodef::Bool = false,
				rep_agent::Bool = false
				)
	ψ = IES
	γ = 0.
	if EpsteinZin == true
		γ = RRA
	end

	σmin, σmax = 0.1, 1.25
	if rep_agent
		σϵ = 0.0001
		Nϵ = 2
		Nσ = 2
		σmin, σmax = 0.01, 0.02
	end

	## Prepare discretized processes
	function quarterlize_AR1(ρ, σ)
		ρ4 = ρ^0.25
		σ4 = sqrt(  σ^2 / ( 1 + ρ4^2 + ρ4^4 + ρ4^6 )  )
		return ρ4, σ4
	end
	# Aggregate risk
	z_chain = tauchen(Nz, ρz, σz, 0, 2)
	# z_chain = rouwenhorst(Nz, ρz, σz, 0)
	Pz = z_chain.p
	# zgrid = linspace(minimum(z_chain.state_values), maximum(z_chain.state_values), Nz)
	zgrid = z_chain.state_values

	# Idiosyncratic risk
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
	# ϵ_chain = rouwenhorst(Nϵ, ρϵ, σϵ, 0)
	Pϵ = ϵ_chain.p
	ϵgrid = ϵ_chain.state_values

	wgrid = [0.5; wbar] # linspace(0.75, 1.5, Nw)
	Nw = 2
	pngrid = linspace(0.5, 1.1, Np)
	ζgrid = 1:2
	Nζ = length(ζgrid)

	λϵ = stationary_distributions(ϵ_chain)[1]

	χ = 2.0
	Ξ = dot(exp.(ϵgrid).^(1.0/χ), λϵ)^χ
	θL = (1.0-τ) * Ξ

	α_T = 0.6
	α_N = 0.6

	ϑ = 0.88 # Straight from Anzoategui

	# α_T = 0.75
	# α_N = 0.75

	μ_anzo = 0.74 # Taken straight from Anzoategui, from Stockman and Tesar (1995)
	ω_anzo = 0.8  # Taken from Anzoategui, targets SS output share of nontradables at 88%

	η = μ_anzo
	ϖ = ω_anzo^(1.0/μ_anzo)

	# ϖ = 0.7 * ϖ	

	# Grids for endogenous aggregate states
	Bmax  = 4.0
	Bbar  = Bmax * 0.5
	bgrid = linspace(0.0, Bmax, Nb)
	μgrid = linspace(-2.5, 0.75, Nμ)
	σgrid = linspace(σmin, σmax, Nσ)

	# Prepare grid for cash in hand.
	ωmin	= -0.5
	ωgrid0	= linspace(0.0, (ωmax-ωmin)^curv, Nω).^(1/curv)
	# ωgrid0	= linspace(0.0, (ωmax-ωmin), Nω)
	ωgrid0	= ωgrid0 + ωmin
	ωgrid 	= ωgrid0

	ωgrid_fine	= linspace(0., (ωmax-ωmin), Nω_fine)
	ωgrid_fine	= ωgrid_fine + ωmin

	snodes = [kron(ones(Nϵ,), ωgrid_fine) kron(ϵgrid, ones(Nω_fine,))]

	# Define the basis over the state variables
	# basis = Basis(SplineParams(ωgrid0, 0, order),
	basis = Basis(LinParams(ωgrid, 0),
				  LinParams(ϵgrid, 0),
				  LinParams(bgrid, 0),
				  LinParams(μgrid, 0),
				  LinParams(σgrid, 0),
				  LinParams(wgrid, 0),
				  LinParams(ζgrid, 0),
				  LinParams(zgrid, 0))
	s, _ = nodes(basis)
	Nω, Ns = size(ωgrid, 1), size(s, 1)

	Jgrid = gridmake(1:Nb, 1:Nμ, 1:Nσ, 1:Nw, 1:Nζ, 1:Nz)

	# Compute the basis matrix and expectations matrix
	bs = BasisMatrix(basis, Direct(), s, [0 0 0 0 0 0 0 0])
	Φ = convert(Expanded, bs).vals[1]

	λ = ones(Nω_fine*Nϵ)
	λ = λ/sum(λ)

	ϕa_ext = zeros(Nω, Nϵ, Nb, Nμ, Nσ, Nw, Nζ, Nz, Np)
	ϕb_ext = zeros(Nω, Nϵ, Nb, Nμ, Nσ, Nw, Nζ, Nz, Np)
	ϕe_ext = zeros(Nω, Nϵ, Nb, Nμ, Nσ, Nw, Nζ, Nz, Np)
	ϕc_ext = zeros(Nω, Nϵ, Nb, Nμ, Nσ, Nw, Nζ, Nz, Np)

	μ = Array{Float64}(Nb, Nμ, Nσ, Nw, Nζ, Nz)
	for (jμ, μv) in enumerate(μgrid)
		μ[:,jμ,:,:,:,:,:] = μv + 0.5 * ( mean(μgrid) - μv )
	end
	μ = reshape(μ, Nb*Nμ*Nσ*Nw*Nζ*Nz)
	σ = Array{Float64}(Nb, Nμ, Nσ, Nw, Nζ, Nz)
	for (jσ, σv) in enumerate(σgrid)
		σ[:,:,jσ,:,:,:,:] = σv + 0.5 * ( mean(σgrid) - σv )
	end
	σ = reshape(σ, Nb*Nμ*Nσ*Nw*Nζ*Nz)
	μ′ = Array{Float64, 3}(Nb*Nμ*Nσ*Nw*Nζ*Nz, Nz, 2)
	σ′ = Array{Float64, 3}(Nb*Nμ*Nσ*Nw*Nζ*Nz, Nz, 2)
	for j = 1:Nz, jj = 1:2
		μ′[:, j, jj] = μ
		σ′[:, j, jj] = σ
	end
	w′ = Array{Float64}(Nb, Nμ, Nσ, Nw, Nζ, Nz)
	for (jw, wv) in enumerate(wgrid)
		w′[:,:,:,jw,:,:] = wv
	end
	w′ = reshape(w′, Nb*Nμ*Nσ*Nw*Nζ*Nz)

	# Debt parameters
	ρ = 0.05 # Target average maturity of 7 years: ~0.05 at quarterly freq
	κ = ρ + r_star

	# State functions
	Ld = ones(Nb*Nμ*Nσ*Nw*Nζ*Nz)
	T  = ones(Nb*Nμ*Nσ*Nw*Nζ*Nz) * 0.05
	qʰ = ones(Nb*Nμ*Nσ*Nw*Nζ*Nz) / (1.0+r_star)
	qᵍ = zeros(Nb*Nμ*Nσ*Nw*Nζ*Nz)
	spread = zeros(Nb*Nμ*Nσ*Nw*Nζ*Nz)

	pN 		  = Array{Float64}(Nb, Nμ, Nσ, Nw, Nζ, Nz)
	repay 	  = Array{Float64}(Nb, Nμ, Nσ, Nw, Nζ, Nz, Nz)
	wage 	  = Array{Float64}(Nb, Nμ, Nσ, Nw, Nζ, Nz)
	spending  =	Array{Float64}(Nb, Nμ, Nσ, Nw, Nζ, Nz)
	issuance  =	Array{Float64}(Nb, Nμ, Nσ, Nw, Nζ, Nz)
	def_thres =	Array{Float64}(Nb, Nμ, Nσ, Nw, Nζ, Nz)
	output	  = Array{Float64}(Nb, Nμ, Nσ, Nw, Nζ, Nz)
	for (jz, zv) in enumerate(zgrid)
		pN[:,:,:,:,:,jz] = mean(pngrid) - 0.1 * zv
		output[:,:,:,:,:,jz] = exp(zv)
		spending[:,:,:,:,:,jz] = 0.1 - 0.25 * zv
		for (jb, bv) in enumerate(bgrid)
			issuance[jb,:,:,:,1,jz] = bv - 0.25 * zv + 0.1 * (Bbar-bv)
			issuance[jb,:,:,:,2,jz] = bv
		end
		for (jζ, ζv) in enumerate(ζgrid)
			def = (ζv != 1.0)
			for (jw, wv) in enumerate(wgrid)
				wage[:,:,:,jw,jζ,jz] = max(exp(zv) * (1.0 - Δ * def), γw*wv)
			end
		end
		repay[:,:,:,:,:,:,jz] = 1.0 - (zv <= zgrid[1])
		def_thres[:,:,:,:,:,jz] = zgrid[1]
		# def_thres[:,:,:,:,:,jz] = -Inf
	end
	pN	 		= reshape(pN, 	 	 Nb*Nμ*Nσ*Nw*Nζ*Nz)
	repay	 	= reshape(repay, 	 Nb*Nμ*Nσ*Nw*Nζ*Nz*Nz)
	wage	 	= reshape(wage, 	 Nb*Nμ*Nσ*Nw*Nζ*Nz)
	spending 	= reshape(spending,	 Nb*Nμ*Nσ*Nw*Nζ*Nz)
	issuance 	= min.(max.(reshape(issuance,  Nb*Nμ*Nσ*Nw*Nζ*Nz), minimum(bgrid)), maximum(bgrid))
	def_thres 	= reshape(def_thres, Nb*Nμ*Nσ*Nw*Nζ*Nz)
	output 		= reshape(output, Nb*Nμ*Nσ*Nw*Nζ*Nz)
	
	welfare   	= zeros(Nb, Nμ, Nσ, Nw, Nζ, Nz)
	welfare   	= reshape(welfare, Nb*Nμ*Nσ*Nw*Nζ*Nz)
	
	profits 	= output - wage .* Ld


	ϕa = zeros(Nω, Nϵ, Nb, Nμ, Nσ, Nw, Nζ, Nz)
	ϕb = zeros(Nω, Nϵ, Nb, Nμ, Nσ, Nw, Nζ, Nz)
	ϕe = zeros(Nω, Nϵ, Nb, Nμ, Nσ, Nw, Nζ, Nz)
	ϕc = zeros(Nω, Nϵ, Nb, Nμ, Nσ, Nw, Nζ, Nz)

	vf = Array{Float64}(Nω, Nϵ, Nb, Nμ, Nσ, Nw, Nζ, Nz)
	for js in 1:size(Jgrid,1)
		jb = Jgrid[js, 1]
		jμ = Jgrid[js, 2]
		jσ = Jgrid[js, 3]
		jw = Jgrid[js, 4]
		jζ = Jgrid[js, 5]
		jz = Jgrid[js, 6]

		pNv = pN[js]
		pC = (ϖ * pNv.^(1.0-η) + (1.0-ϖ)).^(1.0/(1.0-η))

		Π = profits[js]

		wv = exp(zgrid[jz])
		for (jϵ, ϵv) in enumerate(ϵgrid), (jω, ωv) in enumerate(ωgrid)

			Y = exp(ϵv) * (wv * (1.0-τ) + Π) + (ωv-ωmin)
			Y = max.(Y, 1e-4)
			a = Y * 0.5 ./ (1./(1.0+r_star))
			c = Y * 0.5 ./ pC
			ϕa[jω, jϵ, jb, jμ, jσ, jw, jζ, jz] = a
			ϕb[jω, jϵ, jb, jμ, jσ, jw, jζ, jz] = 0.
			ϕc[jω, jϵ, jb, jμ, jσ, jw, jζ, jz] = c
			ut = log(c)
			γ == 1? Void: ut = c^(1.0-γ) / (1.0-γ)
			if EpsteinZin
				vf[jω, jϵ, jb, jμ, jσ, jw, jζ, jz] = c
			else
				vf[jω, jϵ, jb, jμ, jσ, jw, jζ, jz] = ut / (1.0 - β)
			end
		end
	end

	if nodef
		repay = ones(repay)
	end

	outer_dists = [1.]

	return Hank(β, γ, ψ, EpsteinZin, γw, θL, χ, Ξ, ρ, κ, r_star, tax, η, ϖ, α_T, α_N, ϑ, ϕa, ϕb, ϕe, ϕc, ϕa_ext, ϕb_ext, ϕe_ext, ϕc_ext, vf, ρϵ, σϵ, ρz, σz, Nω, Nϵ, Nb, Nμ, Nσ, Nw, Nζ, Nz, Ns, Nω_fine, Pϵ, Pz, λ, λϵ, ℏ, θ, Δ, #curv, order,
		ωmin, ωmax, ωgrid0, ωgrid, ϵgrid, bgrid, μgrid, σgrid, wgrid, ζgrid, zgrid, s, Jgrid, pngrid, basis, bs, #Φ,
		ωgrid_fine, snodes, μ′, σ′, w′, repay, welfare, τ, T, issuance, def_thres, output, profits, spending, wage, Ld, qʰ, qᵍ, spread, pN, outer_dists, upd_tol)
end

function iterate_qᵍ!(h::Hank; verbose::Bool=false)
	""" Uses the government's repayment function to set the price of debt """
	dist, iter = 10.0, 0
	tol, maxiter = 1e-12, 1500

	init_t = time()

	coupon = h.κ * (1.0 - h.tax)
	qᵍ_mat = reshape(h.qᵍ, h.Nb, h.Nμ, h.Nσ, h.Nw, h.Nζ, h.Nz)
	rep_mat = reshape(h.repay, h.Nb, h.Nμ, h.Nσ, h.Nw, h.Nζ, h.Nz, h.Nz)

	qᵍ = ones(qᵍ_mat)
	spread = zeros(qᵍ)
	while dist > tol && iter < maxiter
		old_q  = copy(qᵍ)
		knots  = (h.bgrid, h.μgrid, h.σgrid, h.wgrid, h.ζgrid, 1:h.Nz)
		itp_qᵍ = interpolate(knots, qᵍ, (Gridded(Linear()), Gridded(Linear()), Gridded(Linear()), Gridded(Linear()), NoInterp(), NoInterp()))

		for js in 1:size(h.Jgrid,1)
			jb = h.Jgrid[js, 1]
			jμ = h.Jgrid[js, 2]
			jσ = h.Jgrid[js, 3]
			jw = h.Jgrid[js, 4]
			jζ = h.Jgrid[js, 5]
			jz = h.Jgrid[js, 6]

			ζv, zv = h.ζgrid[jζ], h.zgrid[jz]

			exp_rep = rep_mat[jb, jμ, jσ, jw, jζ, jz, :]

			jdefault = (ζv != 1.0)

			bpv = h.issuance[js]
			wpv = h.w′[js]
			thres = h.def_thres[js]

			E_rep, check = 0.0, 0.0
			E_fullrep = 0.
			if jdefault == false
				for (jzp, zpv) in enumerate(h.zgrid)
					ζpv = 1.0
					μpv = h.μ′[js, jzp, 1]
					σpv = h.σ′[js, jzp, 1]
					E_rep += h.Pz[jz, jzp] * (coupon + (1.0-h.ρ) * itp_qᵍ[bpv, μpv, σpv, wpv, ζpv, jzp]) * exp_rep[jzp]
					E_fullrep += h.Pz[jz, jzp] * (coupon + (1.-h.ρ) * itp_qᵍ[bpv, μpv, σpv, wpv, ζpv, jzp]) * exp_rep[jzp]
					ζpv = 2.0
					μpv = h.μ′[js, jzp, 2]
					σpv = h.σ′[js, jzp, 2]
					E_rep += h.Pz[jz, jzp] * (1.0-h.ℏ) * (1.0-h.ρ) * itp_qᵍ[(1.0 - h.ℏ)*bpv, μpv, σpv, wpv, ζpv, jzp] * (1.-exp_rep[jzp])
					E_fullrep += h.Pz[jz, jzp] * (coupon + (1.-h.ρ) * itp_qᵍ[bpv, μpv, σpv, wpv, ζpv, jzp]) * (1.-exp_rep[jzp])
					check += h.Pz[jz, jzp] * (1.-exp_rep[jzp])
					check += h.Pz[jz, jzp] * exp_rep[jzp]
				end
			else
				for (jzp, zpv) in enumerate(h.zgrid)
					ζ_reent = 1.0
					μpv = h.μ′[js, jzp, 1]
					σpv = h.σ′[js, jzp, 1]
					E_rep += h.Pz[jz, jzp] * (coupon + (1.0-h.ρ) * itp_qᵍ[bpv, μpv, σpv, wpv, ζ_reent, jzp]) * h.θ
					E_fullrep += h.Pz[jz, jzp] * (coupon + (1.-h.ρ) * itp_qᵍ[bpv, μpv, σpv, wpv, ζ_reent, jzp]) * h.θ
					check += h.Pz[jz, jzp] * h.θ
					ζ_cont = 2.0
					μpv = h.μ′[js, jzp, 2]
					σpv = h.σ′[js, jzp, 2]
					E_rep += h.Pz[jz, jzp] * (1.0-h.ρ) * itp_qᵍ[bpv, μpv, σpv, wpv, ζ_cont, jzp] * (1.0 - h.θ)
					E_fullrep += h.Pz[jz, jzp] * (coupon + (1.-h.ρ) * itp_qᵍ[bpv, μpv, σpv, wpv, ζ_cont, jzp]) * (1.-h.θ)
					check += h.Pz[jz, jzp] * (1.0 - h.θ)
				end
			end

			isapprox(check, 1.0) || print_save("WARNING: wrong transitions in update_qᵍ!")
			qᵍ[jb, jμ, jσ, jw, jζ, jz] = E_rep / (1.0 + h.r_star)
			spread[jb, jμ, jσ, jw, jζ, jz] = E_fullrep / qᵍ[jb, jμ, jσ, jw, jζ, jz] - (1+h.r_star)
		end
		iter += 1
		dist = sum( (qᵍ - old_q).^2 ) / sum(old_q.^2)
	end

	h.qᵍ = reshape(qᵍ, h.Nb*h.Nμ*h.Nσ*h.Nw*h.Nζ*h.Nz)
	h.spread = reshape(spread, h.Nb*h.Nμ*h.Nσ*h.Nw*h.Nζ*h.Nz)

	h.spread = 1.0./h.qᵍ - 1.0

	if verbose
		end_t = time()
		if dist <= tol
			print_save("Updated prices after $iter iterations in $(time_print(end_t-init_t))")
		else
			warn("Iteration on qᵍ aborted at distance $(@sprintf("%.3g",dist)) after $(time_print(end_t-init_t))")
		end
	end

	Void
end

function compute_netexports(h::Hank)
	itp_ϕc = make_itp(h, h.ϕc; agg=false)

	aggC_T = zeros(size(h.Jgrid, 1))
	supply_T = zeros(size(h.Jgrid, 1))
	for js in 1:size(h.Jgrid, 1)
		bv = h.bgrid[h.Jgrid[js, 1]]
		μv = h.μgrid[h.Jgrid[js, 2]]
		σv = h.σgrid[h.Jgrid[js, 3]]
		wv = h.wgrid[h.Jgrid[js, 4]]
		jζ = h.Jgrid[js, 5]
		jz = h.Jgrid[js, 6]
		ζv = h.ζgrid[jζ]
		zv = h.zgrid[jz]

		pNv = h.pN[js]
		pC = price_index(h, pNv)

		aggC = integrate_itp(h, bv, μv, σv, wv, jζ, jz, itp_ϕc)
		aggC_T[js] = aggC  * h.ϖ * (1./pC)^(-h.η)

		wt = h.wage[js]
		Ld_N, Ld_T  = labor_demand(h, wt, zv, ζv, pNv; get_both=true)
		supply_T[js] = TFP_T(zv, h.Δ, ζv) * Ld_T^(h.α_N)
	end

	demand_G_T = h.spending * (1.-h.ϑ)
	NX = supply_T - aggC_T - demand_G_T

	return NX
end

function update_fiscalrules!(h::Hank)
	unemp = (1.- h.Ld)*100.
	unemp2= unemp.^2
	BoY   = h.bgrid[h.Jgrid[:, 1]] ./ (4 * h.output) * 100
	BoY2  = BoY.^2
	spread= h.spread * 100 # measure in percent
	spr2  = spread.^2
	NX 	  = compute_netexports(h)./(1 * h.output) * 100
	NX2   = NX.^2

	coef_g = [20.6746151785  0.0307028678  0.0021122640  0.0095750707 -0.0002004299  0.0090837466 -0.0001310273]

	coef_B = [ 1.0786829981  0.3342813521  0.0001223178 -0.0102040316  0.0001494438  0.0461321365 -0.0014158229 ]
	g = [ ones(unemp) unemp unemp2 BoY BoY2 NX NX2 ] * coef_g' / 100
	net_iss = [ ones(unemp) unemp unemp2 BoY BoY2 NX NX2 ] * coef_B' / 100

	h.spending = max.(min.(vec(g), 0.35), 0.) .* (1 * h.output)
	h.issuance = min.(0.10,max.(0., vec(net_iss))) .* (4 * h.output) + (1.-h.ρ)*h.bgrid[h.Jgrid[:, 1]]

	def_states = h.ζgrid[h.Jgrid[:, 5]] .!= 1.0
	h.issuance[def_states] = (1.-h.ρ) * h.bgrid[h.Jgrid[def_states, 1]]

	h.issuance = min.(h.issuance, maximum(h.bgrid))

	Void
end

price_index(h::Hank, pN) = (h.ϖ * pN.^(1.0-h.η) + (1.0-h.ϖ)).^(1.0/(1.0-h.η))

function govt_bc(h::Hank, wage_bill)
	"""
	Computes lump-sum taxes from the government's budget constraint.
	`wage_bill` here is w * Lᵈ
	"""
	qᵍ_vec = h.qᵍ
	def_states = h.ζgrid[h.Jgrid[:, 5]] .!= 1.0

	B′ = h.issuance
	B  = h.bgrid[h.Jgrid[:, 1]]

	coupons = (1.0 - def_states) .* h.κ .* B
	g 		= h.spending
	inc_tax = h.τ * wage_bill
	net_iss = qᵍ_vec .* (B′ - (1.0 - h.ρ) .* B)

	T_vec = coupons + g - inc_tax - net_iss
	T_mat = reshape(T_vec, h.Nb, h.Nμ, h.Nσ, h.Nw, h.Nζ, h.Nz)
	return T_mat
end

function _unpackstatefs(h::Hank)

	wL = h.Ld .* h.wage .* (1.0-h.τ)
	jζ = h.Jgrid[:, 5]

	pC = price_index(h, h.pN)

	taxes_mat = govt_bc(h, h.wage .* h.Ld)

	profits_mat = reshape( h.output - h.wage .* h.Ld, h.Nb, h.Nμ, h.Nσ, h.Nw, h.Nζ, h.Nz)

	T_mat = taxes_mat# - profits_mat

	qʰ_mat = reshape(h.qʰ, 	h.Nb, h.Nμ, h.Nσ, h.Nw, h.Nζ, h.Nz)
	qᵍ_mat = reshape(h.qᵍ, 	h.Nb, h.Nμ, h.Nσ, h.Nw, h.Nζ, h.Nz)
	wL_mat = reshape(wL, 	h.Nb, h.Nμ, h.Nσ, h.Nw, h.Nζ, h.Nz)
	pC_mat = reshape(pC, 	h.Nb, h.Nμ, h.Nσ, h.Nw, h.Nζ, h.Nz)

	return qʰ_mat, qᵍ_mat, wL_mat, T_mat, pC_mat, profits_mat
end


function vfi!(h::Hank; tol::Float64=5e-3, verbose::Bool=true, remote::Bool=true, maxiter::Int64=50, bellman_iter::Int64=maxiter)

	print_save("\nSolving household problem: ")
	time_init = time()
	iter = 1
	iter_cycle = 0
	dist, dist_s = 10., 10.

	iterate_qᵍ!(h, verbose = true)
	update_fiscalrules!(h)
	var(h.qʰ) .< 1e-16 || print_save("\nWARNING: qʰ is not constant. $(var(h.qʰ))")
	print_save("\nqᵍ between $(round(minimum(h.qᵍ[h.Jgrid[:,5].==1]),4)) and $(round(maximum(h.qᵍ),4)). risk-free is $(round(mean(h.qʰ),4))")
	print_save("\nspread between $(floor(Int,10000*minimum(h.spread[h.Jgrid[:,5].==1]))) bps and $(floor(Int,10000*maximum(h.spread[h.Jgrid[:,5].==1]))) bps")

	upd_η = 0.5

	dist_statefuncs = Matrix{Float64}(maxiter, 3)
	dist_LoMs = Matrix{Float64}(maxiter, 2)

	μ′_old = copy(h.μ′)
	σ′_old = copy(h.σ′)

	print_save("\nIteration $iter")
	t_old = time()
	while dist > tol && iter < maxiter
		iter_cycle += 1

		qʰ_mat, qᵍ_mat, wL_mat, T_mat, pC_mat, Π_mat = _unpackstatefs(h);
		# println(mean(T_mat))
		v_old = copy(h.vf)
		if iter_cycle <= 5 || iter_cycle % 3 == 0
			warnc0 = bellman_iteration!(h, qʰ_mat, qᵍ_mat, wL_mat, T_mat, pC_mat, Π_mat; resolve=true)
		else
			bellman_iteration!(h, qʰ_mat, qᵍ_mat, wL_mat, T_mat, pC_mat, Π_mat; resolve=false)
		end
		v_new = copy(h.vf)

		dist = sqrt.(sum( (v_new - v_old).^2 )) / sqrt.(sum(v_old.^2))
		norm_v = sqrt.(sum(v_old.^2))
		if verbose && (iter_cycle % 20 == 0 || dist < h.upd_tol)
			t_new = time()
			print_save("\nd(v, v′) = $(@sprintf("%0.3g",dist)) at ‖v‖ = $(@sprintf("%0.3g",norm_v)) after $(time_print(t_new-t_old)) and $iter_cycle iterations ")
			print_save(Dates.format(now(), "HH:MM"))
		end

		if dist < h.upd_tol
			consw = 100*floor(Int,mean(warnc0))
			consw > 0? msg = "\nWARNING: ": msg="\n"
			msg *= "Can't affort consumption $(consw)% of the time"
			print_save(msg)

			t1 = time()
			extend_state_space!(h, qʰ_mat, qᵍ_mat, T_mat)
			print_save(": done in $(time_print(time()-t1))")
			t1 = time()

			# plot_hh_policies(h, remote = remote)
			# plot_hh_policies_b(h, remote = remote)
			# plot_hh_policies_z(h, remote = remote)
			# plot_labor_demand(h, remote = remote)

			print_save("\nUpdating functions of the state")

			exc_dem_prop, exc_sup_prop, mean_excS, max_excS, dists = update_state_functions!(h, upd_η)

			dist_statefuncs[iter, :] = dists

			# plot_state_funcs(h, remote = remote)
			# plot_nontradables(h, remote = remote)
			print_save(": done in $(time_print(time()-t1))")
			t1 = time()

			print_save("\nStates with exc supply, demand = $(round(100*exc_sup_prop,2))%, $(round(100*exc_dem_prop,2))%")
			print_save("\nAverage, max exc supply = $(@sprintf("%0.3g",mean_excS)), $(@sprintf("%0.3g",max_excS))")

			_ = update_grids_pw!(h, exc_dem_prop, exc_sup_prop)
			# update_grids!(h, new_wgrid=new_wgrid)

			print_save("\nNew pN_grid = [$(@sprintf("%0.3g",minimum(h.pngrid))), $(@sprintf("%0.3g",maximum(h.pngrid)))]")
			print_save("\nw_grid = [$(@sprintf("%0.3g",minimum(h.w′))), $(@sprintf("%0.3g",maximum(h.w′)))]")


			print_save("\nDistance in state functions: (dw,dpN,dLd) = ($(@sprintf("%0.3g",mean(dists[1]))),$(@sprintf("%0.3g",mean(dists[2]))),$(@sprintf("%0.3g",mean(dists[3]))))")
			dist_s = maximum(dists)

			dist_exp, new_μgrid, new_σgrid = update_expectations!(h, 0.5 * upd_η)
			# dist_exp = [0. 0.]
			dist_LoMs[iter, :] = dist_exp

			update_grids!(h, new_μgrid = new_μgrid, new_σgrid = new_σgrid)
			print_save("\nDistance in expectations: (dμ,dσ) = ($(@sprintf("%0.3g",mean(dist_exp[1]))),$(@sprintf("%0.3g",mean(dist_exp[2]))))")
			print_save("\nNew μ_grid = [$(@sprintf("%0.3g",minimum(h.μgrid))), $(@sprintf("%0.3g",maximum(h.μgrid)))]")
			print_save("\nNew σ_grid = [$(@sprintf("%0.3g",minimum(h.σgrid))), $(@sprintf("%0.3g",maximum(h.σgrid)))]")

			# plot_LoM(h, remote = remote)

			dist_s = max(dist_s, maximum(dist_exp))
			print_save("\nGrids and expectations updated in $(time_print(time()-t1))")

			# plot_convergence(dist_statefuncs, dist_LoMs, iter, remote = remote)


			iter += 1
			iter_cycle = 0
			print_save("\n\nIteration $iter")
			
			iterate_qᵍ!(h)
			update_fiscalrules!(h)
			var(h.qʰ) .< 1e-16 || print_save("\nWARNING: qʰ is not constant. $(var(h.qʰ))")
			print_save("\nqᵍ between $(round(minimum(h.qᵍ),4)) and $(round(maximum(h.qᵍ),4)). risk-free is $(round(mean(h.qʰ),4))")
			print_save("\nspread between $(floor(Int,10000*minimum(h.spread))) bps and $(floor(Int,10000*maximum(h.spread))) bps")

			h.upd_tol = max(exp(0.85*log(1+h.upd_tol))-1, 1e-6)
			print_save("\nNew update tolerance = $(@sprintf("%0.3g",h.upd_tol))")
			t_old = time()
		end

		dist = max(dist, dist_s)
		if iter % 10 == 0 && !isnan(sum(h.vf)) && !isnan(sum(h.ϕc))
			save(pwd() * "/../../hank.jld", "h", h)
		end

		if isnan.(dist) && iter > 1
			error("NaN encountered")
		end

	end
	# plot_gov_welf(h; remote = remote)
	# plot_aggcons(h; remote = remote)
	# plot_govt_reaction(h; remote = remote)
	# plot_debtprice(h; remote = remote)
	save(pwd() * "/../../hank.jld", "h", h)

	if dist <= tol
		print_save("\nConverged in $iter iterations. ")
	else
		print_save("\nStopping at distance $(@sprintf("%0.3g",dist)), target was $(@sprintf("%0.3g",tol)).")
	end

	print_save("\nTotal time: $(time_print(time()-time_init))\n")

	Void
end
