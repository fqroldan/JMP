include("type_def.jl")

type Path
	data::Matrix{Float64}
	n::Dict{Symbol,Int64}
	y::Function
end
function Path(; T::Int64 = 1)
	n = Dict(
		:B => 1,
		:μ => 2,
		:σ => 3,
		:w => 4,
		:ζ => 5,
		:z => 6,
		)
	data = Matrix{Float64}(T, length(n))

	y(t::Int64, sym::Symbol) = data[t, n[sym]]

	return Path(data, n, y)
end
function fill_path!(p::Path, t::Int64; B::Float64=-Inf, μ::Float64=-Inf, σ::Float64=-Inf, w::Float64=-Inf, ζ::Float64=-Inf, z::Float64=-Inf)
	0 < t <= size(p.data, 1) || throw("t out of bounds")
	B != -Inf? p.data[t, p.n[:B]] = B: Void
	μ != -Inf? p.data[t, p.n[:μ]] = μ: Void
	σ != -Inf? p.data[t, p.n[:σ]] = σ: Void
	w != -Inf? p.data[t, p.n[:w]] = w: Void
	ζ != -Inf? p.data[t, p.n[:ζ]] = ζ: Void
	z != -Inf? p.data[t, p.n[:z]] = z: Void
	Void
end
function trim_path!(p::Path, T_burnin::Int64)
	p.data = p.data[T_burnin+1:end, :]
	Void
end
	
function iter_simul!(h::Hank, p::Path, t, jz_series, itp_ϕa, itp_ϕb, itp_ϕc, itp_B′, itp_G, itp_pN, itp_qᵍ, itp_Zthres, λt, Qϵ)
	# Enter with a state B, μ, σ, w0, ζ, z.
	# h.zgrid[jz] must equal p.y(t, :z)
	# B, ζ, and z are decided at the end of the last period
	jz = jz_series[t]

	Bt = p.y(t, :B)
	μt = p.y(t, :μ)
	σt = p.y(t, :σ)
	w0 = p.y(t, :w)
	ζt = Int(p.y(t, :ζ))
	zt = p.y(t, :z)

	println("$([Bt, μt, σt, w0, ζt, zt])")

	Bprime 	= itp_B′[Bt, μt, σt, w0, ζt, zt]
	G 		= itp_G[Bt, μt, σt, w0, ζt, zt]
	pNg 	= itp_pN[Bt, μt, σt, w0, ζt, zt]
	thres 	= itp_Zthres[Bt, μt, σt, w0, ζt, zt]

	# Find pN at the current state. Deduce w, L, Π, T.
	pNmin, pNmax = minimum(h.pngrid), maximum(h.pngrid)
	jdef = (ζt != 1)
	results, _ = find_prices(h, itp_ϕc, G, Bprime, pNg, pNmin, pNmax, Bt, μt, σt, w0, ζt, jz, jdef)

	wt, pN, Ld, output = results

	# Integrate the household's policy functions to get μ′, σ′

	ϕa = itp_ϕa[h.ωgrid_fine, 1:h.Nϵ, Bt, μt, σt, w0, ζt, zt]
	ϕb = itp_ϕb[h.ωgrid_fine, 1:h.Nϵ, Bt, μt, σt, w0, ζt, zt]

	ϕa = reshape(ϕa, h.Nω_fine*h.Nϵ)
	ϕb = reshape(ϕb, h.Nω_fine*h.Nϵ)

	a  = dot(λt, ϕa)
	a2 = dot(λt, ϕa.^2)
	b  = dot(λt, ϕb)
	b2 = dot(λt, ϕb.^2)
	ab = dot(λt, ϕa.*ϕb)

	var_a  = a2 - a^2
	var_b  = b2 - b^2
	cov_ab = ab - a*b

	μ′, σ′, q′ = compute_stats_logN(h, ζt, a, b, var_a, var_b, cov_ab, itp_qᵍ, Bprime, wt, thres)
	# μ′, σ′, q′ = new_expectations(h, itp_ϕa, itp_ϕb, itp_qᵍ, Bprime, wt, thres, Bt, μt, σt, w0, ζt, zt, jdef) # This would assume that λₜ is lognormal
	println(q′)

	# Draw z and the reentry shock for tomorrow, deduce ζ and correct B, μ, and σ as needed, and update the distribution
	probs = cumsum(h.Pz[jz,:])
	jzp = 1 + findfirst(rand() .> probs)

	zprime = h.zgrid[jzp]

	μprime = μ′[jzp, 1]
	σprime = σ′[jzp, 1]
	qprime = q′[jzp, 1]
	if jdef
		prob_reentry = h.θ
		reentry = (rand() <= prob_reentry)
		if reentry
			ζprime = 1.0
			R = h.κ + (1.0-h.ρ) * qprime
		else
			ζprime = 2.0
			μprime = μ′[jzp, 2]
			σprime = σ′[jzp, 2]
			qprime = q′[jzp, 2]
			R = h.κ + (1.0-h.ρ) * qprime
		end
	else
		if zprime <= thres
			ζprime = 2.0
			Bprime = (1.0 - h.ℏ) * Bprime
			R = (1.0-h.ℏ)*(1.0-h.ρ) * qprime
		else
			ζprime = 1.0
			R = h.κ + (1.0-h.ρ) * qprime
		end
	end

	basis = Basis(LinParams(h.ωgrid_fine, 0))
	Qω = BasisMatrix(basis, Expanded(), ϕa + R*ϕb, 0).vals[1]
	Q = row_kron(Qϵ, Qω)

	λprime = Q' * λt

	jz_series[t+1] = jzp

	# Fill the path for next period
	fill_path!(p,t+1; B = Bprime, μ = μprime, σ = σprime, w = wt, ζ = ζprime, z = zprime)

	return λprime
end

function simul(h::Hank; simul_length::Int64=1, burn_in::Int64=0)
	# Setup
	T = burn_in + simul_length
	p = Path(T = T)
	B0, μ0, σ0, w0, ζ0, z0 = mean(h.bgrid), mean(h.μgrid), mean(h.σgrid), mean(h.wgrid), h.ζgrid[1], mean(h.zgrid)
	fill_path!(p,1; B = B0, μ = μ0, σ = σ0, w = w0, ζ = ζ0, z = z0)

	function itp_all(h::Hank, Y::Array{Float64})
		all_knots = (h.ωgrid, 1:h.Nϵ, h.bgrid, h.μgrid, h.σgrid, h.wgrid, 1:h.Nζ, h.zgrid)
		
		return interpolate(all_knots, Y, (Gridded(Linear()), NoInterp(), Gridded(Linear()), Gridded(Linear()), Gridded(Linear()), Gridded(Linear()), NoInterp(), Gridded(Linear())))
	end

	itp_ϕa = itp_all(h, h.ϕa)
	itp_ϕb = itp_all(h, h.ϕb)
	itp_ϕc = itp_all(h, h.ϕc)

	function itp_agg(h::Hank, Y::Vector{Float64})
		Y_mat = reshape(Y, h.Nb, h.Nμ, h.Nσ, h.Nw, h.Nζ, h.Nz)
		agg_knots = (h.bgrid, h.μgrid, h.σgrid, h.wgrid, 1:h.Nζ, h.zgrid)
		return interpolate(agg_knots, Y_mat, (Gridded(Linear()), Gridded(Linear()), Gridded(Linear()), Gridded(Linear()), NoInterp(), Gridded(Linear())))
	end

	itp_B′		= itp_agg(h, h.issuance)
	itp_G		= itp_agg(h, h.spending)
	itp_pN		= itp_agg(h, h.pN)
	itp_qᵍ 		= itp_agg(h, h.qᵍ)
	itp_Zthres	= itp_agg(h, h.def_thres)

	jz_series = Vector{Int64}(T)
	jz_series[1] = 1

	# Initialize objects for iterating the distribution
	# λ = ones(h.Nω_fine*h.Nϵ) / (h.Nω_fine*h.Nϵ)
	λ = zeros(h.Nω_fine*h.Nϵ)
	js = 0
	for (jϵ, ϵv) in enumerate(h.ϵgrid)
		for jω = 1:length(h.ωgrid_fine)-1
			js += 1
			ωv  = h.ωgrid_fine[jω]
			ω1v = h.ωgrid_fine[jω+1]
			ωmv = 0.5*(ωv+ω1v)

			λ[js] = pdf(LogNormal(μ0, σ0), ωmv-h.ωmin) * h.λϵ[jϵ] * (ω1v - ωv)
		end
	end
	λ = λ / sum(λ)
	Qϵ = kron(h.Pϵ, ones(h.Nω_fine,1))

	# Simulate
	for t in 1:T
		# Advance one step
		λ = iter_simul!(h, p, t, jz_series, itp_ϕa, itp_ϕb, itp_ϕc, itp_B′, itp_G, itp_pN, itp_qᵍ, itp_Zthres, λ, Qϵ)
	end

	# Keep only after the burn_in period
	trim_path!(p, burn_in)
	jz_series = jz_series[burn_in+1:end]

	# Return stuff
	return p, jz_series
end