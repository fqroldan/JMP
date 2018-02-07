using QuantEcon, BasisMatrices, Interpolations, Optim, NLopt, MINPACK, LaTeXStrings, Distributions

include("hh_pb.jl")

function Hank(;	β = (1.0/1.1)^0.25,
				γ = 2.,
				IES = 2.,
				RRA = 2.,
				γw = 0.99,
				τ = 0.3,
				r_star = 1.02^0.25 - 1.0,
				ωmax = 10.,
				curv = .4,
				income_process = "Mendoza-D'Erasmo",
				EpsteinZin = true,
				order = 3,
				Nω_fine = 1000,
				Nω = 5,
				Nϵ = 4,
				Nμ = 2,
				Nσ = 2,
				Nb = 2,
				Nw = 4,
				Nz = 4,
				ρz = 0.9,
				σz = 0.05,
				ℏ = .25,
				Δ = .05,
				θ = .5,
				Nq = 7,
				Np = 4
				)
	ψ = IES
	if EpsteinZin == true
		γ = RRA
	end
	## Prepare discretized processes
	# Aggregate risk
	z_chain = tauchen(Nz, ρz, σz, 0, 1)
	Pz = z_chain.p
	# zgrid = linspace(minimum(z_chain.state_values), maximum(z_chain.state_values), Nz)
	zgrid = exp.(z_chain.state_values)

	# Idiosyncratic risk
	ρϵ, σϵ = 0., 0.
	if income_process == "Floden-Lindé"
		ρϵ = 0.9136			# Floden-Lindé for US
		σϵ = sqrt(0.0426)	# Floden-Lindé for US
	elseif income_process == "Mendoza-D'Erasmo"
		ρϵ = 0.85 			# Mendoza-D'Erasmo for Spain
		σϵ = 0.2498			# Mendoza-D'Erasmo for Spain
	else
		print_save("ERROR: Must specify an income process")
		throw(error("Must specify an income process"))
	end
	ϵ_chain = tauchen(Nϵ, ρϵ, σϵ, 0, 1)
	Pϵ = ϵ_chain.p
	ϵgrid = exp.(ϵ_chain.state_values)

	wgrid = linspace(0.5, 1.2, Nw)
	pngrid = linspace(0.5, 1.4, Np)
	ζgrid = 1:3
	Nζ = 3

	λϵ = stationary_distributions(ϵ_chain)[1]

	χ = 2.0
	Ξ = dot(ϵgrid.^(1.0/χ), λϵ)^χ
	θL = (1.0-τ) * Ξ

	α_T = 0.6
	α_N = 0.6

	η = 0.74 # Taken straight from Anzoategui, from Stockman and Tesar (1995)
	ϖ = 0.80 # Taken from Anzoategui, targets SS output share of nontradables at 88%

	# Grids for endogenous aggregate states
	bgrid = linspace(0.2, 1.0, Nb)
	μgrid = linspace(-0.2, 0.8, Nμ)
	σgrid = linspace(0.5, 3.0, Nσ)

	# Prepare grid for cash in hand.
	ωmin	= -0.5
	ωgrid0	= linspace(0.0, (ωmax-ωmin)^curv, Nω).^(1/curv)
	ωgrid0	= ωgrid0 + ωmin
	ωgrid 	= ωgrid0

	ωgrid_fine	= linspace(0., (ωmax-ωmin)^curv, Nω_fine).^(1/curv)
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

	ϕa = zeros(Nω, Nϵ, Nb, Nμ, Nσ, Nw, Nζ, Nz)
	ϕb = zeros(Nω, Nϵ, Nb, Nμ, Nσ, Nw, Nζ, Nz)
	ϕc = zeros(Nω, Nϵ, Nb, Nμ, Nσ, Nw, Nζ, Nz)

	vf = Array{Float64}(Nω, Nϵ, Nb, Nμ, Nσ, Nw, Nζ, Nz)
	for js in 1:size(Jgrid,1)
		jb = Jgrid[js, 1]
		jμ = Jgrid[js, 2]
		jσ = Jgrid[js, 3]
		jw = Jgrid[js, 4]
		jζ = Jgrid[js, 5]
		jz = Jgrid[js, 6]

		wv = zgrid[jz]
		for (jϵ, ϵv) in enumerate(ϵgrid), (jω, ωv) in enumerate(ωgrid)
	
			Y = ϵv * wv * (1.0-τ) + (ωv-ωmin)
			c = Y * 0.5
			ϕa[jω, jϵ, jb, jμ, jσ, jw, jζ, jz] = Y * 0.5
			ϕb[jω, jϵ, jb, jμ, jσ, jw, jζ, jz] = Y * 0.0
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


	λ = ones(Nω_fine*Nϵ)
	λ = λ/sum(λ)
	
	ϕa_ext = zeros(Nω, Nϵ, Nb, Nμ, Nσ, Nw, Nζ, Nz, Np)
	ϕb_ext = zeros(Nω, Nϵ, Nb, Nμ, Nσ, Nw, Nζ, Nz, Np)
	ϕc_ext = zeros(Nω, Nϵ, Nb, Nμ, Nσ, Nw, Nζ, Nz, Np)

	μ = Array{Float64}(Nb, Nμ, Nσ, Nw, Nζ, Nz)
	for (jμ, μv) in enumerate(μgrid)
		μ[:,jμ,:,:,:,:,:] = μv + ( mean(μgrid) - μv )/2
	end
	μ = reshape(μ, Nb*Nμ*Nσ*Nw*Nζ*Nz)
	σ = Array{Float64}(Nb, Nμ, Nσ, Nw, Nζ, Nz)
	for (jσ, σv) in enumerate(σgrid)
		σ[:,:,jσ,:,:,:,:] = σv + ( mean(σgrid) - σv )/2
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

	pN 		  = Array{Float64}(Nb, Nμ, Nσ, Nw, Nζ, Nz)
	repay 	  = Array{Float64}(Nb, Nμ, Nσ, Nw, Nζ, Nz)
	wage 	  = Array{Float64}(Nb, Nμ, Nσ, Nw, Nζ, Nz)
	spending  =	Array{Float64}(Nb, Nμ, Nσ, Nw, Nζ, Nz)
	issuance  =	Array{Float64}(Nb, Nμ, Nσ, Nw, Nζ, Nz)
	def_thres =	Array{Float64}(Nb, Nμ, Nσ, Nw, Nζ, Nz)
	for (jz, zv) in enumerate(zgrid)
		pN[:,:,:,:,:,jz] = mean(pngrid) - 0.1 * log(zv)
		spending[:,:,:,:,:,jz] = 0.2 - 0.05 * log.(zv)
		for (jb, bv) in enumerate(bgrid)
			issuance[jb,:,:,:,:,jz] = bv + 0.1 * log.(zv) + 0.1 * (mean(bgrid)-bv)
		end
		for (jζ, ζv) in enumerate(ζgrid)
			def = (ζv != 1.0)
			repay[:,:,:,:,jζ,jz] = 1.0 - ℏ * (ζv == 2.0)
			for (jw, wv) in enumerate(wgrid)
				wage[:,:,:,:,jζ,jz] = max(zv * (1.0 - Δ * def), γw*wv)
			end
		end
		def_thres[:,:,:,:,:,jz] = ifelse(jz == Nz, zgrid[1], 0.0)
	end
	pN	 		= reshape(pN, 	 	 Nb*Nμ*Nσ*Nw*Nζ*Nz)
	repay	 	= reshape(repay, 	 Nb*Nμ*Nσ*Nw*Nζ*Nz)
	wage	 	= reshape(wage, 	 Nb*Nμ*Nσ*Nw*Nζ*Nz)
	spending 	= reshape(spending,	 Nb*Nμ*Nσ*Nw*Nζ*Nz)
	issuance 	= min.(max.(reshape(issuance,  Nb*Nμ*Nσ*Nw*Nζ*Nz), minimum(bgrid)), maximum(bgrid))
	def_thres 	= reshape(def_thres, Nb*Nμ*Nσ*Nw*Nζ*Nz)

	pN = ones(Nb*Nμ*Nσ*Nw*Nζ*Nz) - rand(Nb*Nμ*Nσ*Nw*Nζ*Nz) * 0.25

	ξg = zeros(Ns,)
	ξf = zeros(Ns,)
	ξp = zeros(Ns,)

	
	
	return Hank(β, γ, ψ, EpsteinZin, γw, θL, χ, Ξ, ρ, κ, r_star, η, ϖ, α_T, α_N, ϕa, ϕb, ϕc, ϕa_ext, ϕb_ext, ϕc_ext, vf, ρϵ, σϵ, ρz, σz, Nω, Nϵ, Nb, Nμ, Nσ, Nw, Nζ, Nz, Ns, Nω_fine, Pϵ, Pz, λ, λϵ, ℏ, θ, Δ, curv, order, ωmin, ωmax, ωgrid0, ωgrid, ϵgrid, bgrid, μgrid, σgrid, wgrid, ζgrid, zgrid, s, Jgrid, pngrid, basis, bs, Φ, ωgrid_fine, snodes, μ′, σ′, w′, repay, τ, T, issuance, def_thres, spending, wage, Ld, qʰ, qᵍ, pN)
end

function iterate_qᵍ!(h::Hank; verbose::Bool=false)
	""" Uses the government's repayment function to set the price of debt """
	dist, iter = 10.0, 0
	tol, maxiter = 1e-12, 1500

	init_t = time()

	repay_mat = reshape(h.repay, h.Nb, h.Nμ, h.Nσ, h.Nw, h.Nζ, h.Nz)
	qᵍ_mat 	  = reshape(h.qᵍ, 	 h.Nb, h.Nμ, h.Nσ, h.Nw, h.Nζ, h.Nz)
	jζ_mat 	  = reshape(h.Jgrid[:,5], 	h.Nb, h.Nμ, h.Nσ, h.Nw, h.Nζ, h.Nz)

	qᵍ = ones(qᵍ_mat)
	while dist > tol && iter < maxiter
		old_q = copy(qᵍ)
		R = (jζ_mat .== 1) * h.κ + repay_mat .* ((1.0-h.ρ)*qᵍ)
		knots = (h.bgrid, h.μgrid, h.σgrid, h.wgrid, h.ζgrid, h.zgrid)
		itp_R = interpolate(knots, R, Gridded(Linear()))

		for js in 1:size(h.Jgrid,1)
			jb = h.Jgrid[js, 1]
			jμ = h.Jgrid[js, 2]
			jσ = h.Jgrid[js, 3]
			jw = h.Jgrid[js, 4]
			jζ = h.Jgrid[js, 5]
			jz = h.Jgrid[js, 6]

			ζv, zv = h.ζgrid[jζ], h.zgrid[jz]

			jdefault = (ζv != 1.0)

			qᵍv = qᵍ_mat[jb, jμ, jσ, jw, jζ, jz]

			bpv = h.issuance[js]
			wpv = h.w′[js]
			thres = h.def_thres[js]

			E_rep, check = 0.0, 0.0
			for (jzp, zpv) in enumerate(h.zgrid)
				if jdefault == false
					if zpv <= thres
						ζpv = 2.0
						μpv = h.μ′[js, jzp, 1]
						σpv = h.σ′[js, jzp, 1]
						E_rep += h.Pz[jz, jzp] * itp_R[bpv, μpv, σpv, wpv, ζpv, zpv]
						check += h.Pz[jz, jzp]
					else
						ζpv = 1.0
						μpv = h.μ′[js, jzp, 1]
						σpv = h.σ′[js, jzp, 1]
						E_rep += h.Pz[jz, jzp] * itp_R[bpv, μpv, σpv, wpv, ζpv, zpv]
						check += h.Pz[jz, jzp]
					end
				else
					ζ_reent = 1.0
					μpv = h.μ′[js, jzp, 1]
					σpv = h.σ′[js, jzp, 1]
					E_rep += h.Pz[jz, jzp] * itp_R[bpv, μpv, σpv, wpv, ζ_reent, zpv] * h.θ
					check += h.Pz[jz, jzp] * h.θ
					ζ_cont = 3.0
					μpv = h.μ′[js, jzp, 2]
					σpv = h.σ′[js, jzp, 2]
					E_rep += h.Pz[jz, jzp] * itp_R[bpv, μpv, σpv, wpv, ζ_cont, zpv] * (1.0 - h.θ)
					check += h.Pz[jz, jzp] * (1.0 - h.θ)
				end
			end

			isapprox(check, 1.0) || print_save("WARNING: wrong transitions in update_qᵍ!")
			qᵍ[jb, jμ, jσ, jw, jζ, jz] = E_rep / (1.0 + h.r_star)
		end
		iter += 1
		dist = sum( (qᵍ - old_q).^2 ) / sum(old_q.^2)
	end

	h.qᵍ = reshape(qᵍ, h.Nb*h.Nμ*h.Nσ*h.Nw*h.Nζ*h.Nz)

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

price_index(h::Hank, pN) = (h.ϖ * pN.^(1.0-h.η) + (1.0-h.ϖ)).^(1.0/(1.0-h.η))

function govt_bc(h::Hank, wage_bill)
	"""
	Computes lump-sum taxes from the government's budget constraint.
	`w_vec` here is w * Lᵈ
	"""
	qᵍ_vec = h.qᵍ
	def_states = h.ζgrid[h.Jgrid[:, 5]] .!= 1.0
	
	B′ = h.issuance
	B  = h.bgrid[h.Jgrid[:, 1]]

	coupons = (1.0 - def_states) .* h.κ .* B
	g 		= h.spending
	inc_tax = h.τ * wage_bill
	net_iss = qᵍ_vec .* (B′ - (1.0 - h.ρ) * (1.0 - h.ℏ * def_states) .* B)

	T_vec = coupons + g - inc_tax - net_iss
	T_mat = reshape(T_vec, h.Nb, h.Nμ, h.Nσ, h.Nw, h.Nζ, h.Nz)
	return T_mat
end

function _unpackstatefs(h::Hank)

	wL = h.Ld .* h.wage .* (1.0-h.τ)
	R = h.repay .* (h.κ + (1.0-h.ρ)*h.qᵍ)

	pC = price_index(h, h.pN)

	T_mat = govt_bc(h, h.wage .* h.Ld)

	qʰ_mat = reshape(h.qʰ, 	h.Nb, h.Nμ, h.Nσ, h.Nw, h.Nζ, h.Nz)
	qᵍ_mat = reshape(h.qᵍ, 	h.Nb, h.Nμ, h.Nσ, h.Nw, h.Nζ, h.Nz)
	wL_mat = reshape(wL, 	h.Nb, h.Nμ, h.Nσ, h.Nw, h.Nζ, h.Nz)
	R_mat  = reshape(R, 	h.Nb, h.Nμ, h.Nσ, h.Nw, h.Nζ, h.Nz)
	pC_mat = reshape(pC, 	h.Nb, h.Nμ, h.Nσ, h.Nw, h.Nζ, h.Nz)

	return qʰ_mat, qᵍ_mat, wL_mat, T_mat, R_mat, pC_mat
end


function vfi!(h::Hank; tol::Float64=1e-4, verbose::Bool=true, maxiter::Int64=5000, bellman_iter::Int64=maxiter)

	print_save("\n"*Dates.format(now(), "HH:MM")*"\nSolving household problem: ")
	time_init = time()
	t_old = time_init
	iter = 0
	iter_cycle = 0
	dist, dist_s = 10., 10.
	upd_tol = 0.015

	iterate_qᵍ!(h, verbose = true)

	upd_η = 0.25
	dist_statefuncs = [dist]

	while dist > tol && iter < maxiter
		t_old = time()
		iter += 1
		iter_cycle += 1

		qʰ_mat, qᵍ_mat, wL_mat, T_mat, R_mat, pC_mat = _unpackstatefs(h)

		v_old = copy(h.vf)
		if iter <= 5 || iter % 7 == 0 || iter_cycle == 1
			bellman_iteration!(h, qʰ_mat, qᵍ_mat, wL_mat, T_mat, R_mat, pC_mat; resolve=true)
		else
			bellman_iteration!(h, qʰ_mat, qᵍ_mat, wL_mat, T_mat, R_mat, pC_mat; resolve=true) # false
		end
		v_new = copy(h.vf)
		
		dist = sqrt.(sum( (v_new - v_old).^2 )) / sqrt.(sum(v_old.^2))
		norm_v = sqrt.(sum(v_old.^2))
		if verbose
			t_new = time()
			print_save("\nd(cv, cv′) = $(@sprintf("%0.3g",dist)) at ‖v‖ = $(@sprintf("%0.3g",norm_v)) after $(time_print(t_new-t_old)) and $iter iterations ")
			if iter % 5 == 0
				save(pwd() * "/hank.jld", "h", h)
				plot_hh_policies(h, remote = true)
			end
			print_save(Dates.format(now(), "HH:MM"))
		end

		if dist < upd_tol
			t1 = time()
			extend_state_space!(h, qʰ_mat, qᵍ_mat, T_mat, R_mat)
			print_save(": done in $(time_print(time()-t1))")
			t1 = time()

			save(pwd() * "/hank.jld", "h", h)
			plot_labor_demand(h, remote = true)

			print_save("\nUpdating functions of the state")

			err_N, dists = update_state_functions!(h, upd_η)
			print_save(": done in $(time_print(time()-t1)) \nAverage error in nontradables mkt clearing = $(@sprintf("%0.3g",mean(err_N)))")
			print_save("\nDistance in state functions: (dw,dpN,dLd) = ($(@sprintf("%0.3g",mean(dists[1]))),$(@sprintf("%0.3g",mean(dists[2]))),$(@sprintf("%0.3g",mean(dists[3]))))")
			dist_s = maximum(dists)

			update_expectations!(h, 0.1)
			
			iterate_qᵍ!(h)

			var(h.qʰ) .< 1e-16 || print_save("\nWARNING: qʰ is not constant. $(var(h.qʰ))")
			print_save("\n\nqᵍ between $(round(minimum(h.qᵍ),4)) and $(round(maximum(h.qᵍ),4)). risk-free is $(round(mean(h.qʰ),4))")
			iter_cycle = 0
			upd_tol = update_tolerance(upd_tol, dist_s)
			print_save("\nNew update tolerance = $(@sprintf("%0.3g",upd_tol))")
		end

		dist = max(dist, dist_s)

		if iter % 10 == 0
			# plot_hh_policies(h)
		end

		if isnan.(dist) && iter > 1
			error("NaN encountered")
		end

	end

	if dist <= tol
		print_save("\nConverged in $iter iterations. ")
	else
		print_save("\nStopping at distance $(@sprintf("%0.3g",dist)). ")
	end

	if verbose
		# plot_hh_policies(h)
		# plot_state_funcs(h)
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
		if bas == 1
			step_tol, min_tol = 2.5*10.0^(pot-1), 9*10.0^(pot-1)
		end
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

