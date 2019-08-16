using Random, LinearAlgebra
include("type_def.jl")

function iter_simul!(h::Hank, p::Path, t, jz_series, itp_ϕa, itp_ϕb, itp_ϕc, itp_vf, itp_B′, itp_G, itp_pN, itp_qᵍ, itp_repay, itp_W, λt, Qϵ; phase::String="")
	# Enter with a state B, μ, σ, w0, ζ, z.
	# h.zgrid[jz] must equal getfrompath(p, t, :z)
	# B, ζ, and z are decided at the end of the last period
	jz = jz_series[t]

	Bt = getfrompath(p, t, :B)
	μt = getfrompath(p, t, :μ)
	σt = getfrompath(p, t, :σ)
	# w0 = getfrompath(p, t, :w)
	ξt = getfrompath(p, t, :ξ)
	jξ = findfirst(h.ξgrid .== ξt)
	w0 = h.wbar
	ζt = Int(getfrompath(p, t, :ζ))
	zt = getfrompath(p, t, :z)

	zt == h.zgrid[jz] || print("something wrong with the simulator")
	abs(zt - h.zgrid[jz]) < 1e-8 || throw(error("something wrong with the simulator"))

	
	if t % 100 == 0
		print_save("\n$([Bt, μt, σt, ξt, ζt, zt])")
	end

	Bprime 	= itp_B′(Bt, μt, σt, ξt, ζt, jz)
	G 		= itp_G(Bt, μt, σt, ξt, ζt, jz)
	pNg 	= itp_pN(Bt, μt, σt, ξt, ζt, jz)

	if ζt == 2
		Bprime = Bt
	end

	exp_rep = zeros(h.Nξ, h.Nz)
	itp_repay = extrapolate(itp_repay, Interpolations.Flat())
	for jξp in 1:h.Nξ, jzp in 1:h.Nz
		exp_rep[jξp, jzp] = max(0, min(1, itp_repay(Bt, μt, σt, ξt, ζt, jz, jξp, jzp)))
	end

	# Find pN at the current state. Deduce w, L, Π, T.
	pNmin, pNmax = minimum(h.pngrid), maximum(h.pngrid)
	jdef = (ζt != 1)

	results, _ = find_prices(h, itp_ϕc, G, Bprime, pNg, pNmin, pNmax, Bt, μt, σt, ξt, ζt, jz, jdef)

	wt, pN, Ld, output = results
	profits = output - wt*Ld

	pC = price_index(h, pN)
	Ld_N, Ld_T  = labor_demand(h, wt, zt, ζt, pN; get_both=true)
	supply_T = TFP_T(zt, h.Δ, ζt) * Ld_T^(h.α_N)


	# Govt BC
	if ζt == 1
		coupons = h.κ * Bt
		new_debt = Bprime - (1.0 - h.ρ) * Bt
	else
		coupons = 0
		new_debt = 0
	end
	qg = getfrompath(p, t, :qg)

	lumpsumT = coupons + G - h.τ*wt*Ld - qg*new_debt

	
	if t % 100 == 0
		print_save("\npN = $pN, pN^e = $(pNg), u = $(ceil(100*(1-Ld))) at t = $t")
	end

	def_prob = 0.
	if ζt == 1
		for (jξp, ξvp) in enumerate(h.ξgrid), (jzp, zvp) in enumerate(h.zgrid)
			def_prob += h.Pξ[jξ, jξp] * h.Pz[jz, jzp] * (1.0-exp_rep[jξp, jzp])
		end
	end


	# Integrate the household's policy functions to get μ′, σ′ and C
	ϕa = zeros(h.Nω_fine*h.Nϵ)
	ϕb = zeros(h.Nω_fine*h.Nϵ)
	ϕc = zeros(h.Nω_fine*h.Nϵ)
	Crate = zeros(h.Nω_fine*h.Nϵ)

	js = 0
	ωvjϵ = gridmake(h.ωgrid_fine, 1:h.Nϵ)
	λ_mat  = reshape(λt, h.Nω_fine, h.Nϵ)
	mλω = sum(λ_mat, dims=2)
	cdf_ω = cumsum(mλω[:])

	avgω_fromb = zeros(h.Nω_fine*h.Nϵ)

	b25 = zeros(h.Nω_fine*h.Nϵ)
	b90 = zeros(h.Nω_fine*h.Nϵ)
	
	cdf_ω[1] > 0.25  ? q25 = 1 : q25 = findfirst(cdf_ω .<= 0.25)
	cdf_ω[end] < 0.9 ? q90 = length(cdf_ω) : q90 = findfirst(cdf_ω .>= 0.90)

	adjustment = sum(h.λϵ.*exp.(h.ϵgrid))
	for (jϵ, ϵv) in enumerate(h.ϵgrid), (jω, ωv) in enumerate(h.ωgrid_fine)
		js += 1
		ap = itp_ϕa(ωvjϵ[js,1], ωvjϵ[js,2], Bt, μt, σt, ξt, ζt, jz, pN)
		ϕa[js] = max(h.ωmin, ap)
		bp = itp_ϕb(ωvjϵ[js,1], ωvjϵ[js,2], Bt, μt, σt, ξt, ζt, jz, pN)
		ϕb[js] = max(0.0, bp)
		cc = itp_ϕc(ωvjϵ[js,1], ωvjϵ[js,2], Bt, μt, σt, ξt, ζt, jz, pN)
		ϕc[js] = max(0.0, cc)
		yd = (wt*Ld*(1.0-h.τ) + profits) * exp(ϵv)/adjustment + ωv - lumpsumT
		Crate[js] = ϕc[js] / yd

		avgω_fromb[js] = ωv * max(0.0, bp)

		if jω < q25
			b25[js] += max(0.0, bp)
		elseif jω > q90
			b90[js] += max(0.0, bp)
		end
	end

	C = dot(λt, ϕc)
	CoY = C / output
	CoYd = dot(λt, Crate)

	p25 = dot(λt, b25) / dot(λt, ϕb)
	p90 = dot(λt, b90) / dot(λt, ϕb)

	aggC_T = C  * h.ϖ * (1.0./pC)^(-h.η)
	NX = supply_T - aggC_T - G * (1.0.-h.ϑ)

	Eω_fromb = dot(λt, avgω_fromb) / dot(λt, ϕb)

	fill_path!(p,t, Dict(:P => pN, :Pe => pNg, :Y => output, :L => Ld, :π => def_prob, :w => wt, :G => G, :CoY=> CoY, :CoYd => CoYd, :C => C, :T => lumpsumT, :NX => NX, :p25 => p25, :p90 => p90, :avgω => Eω_fromb))

	a  = dot(λt, ϕa)
	a2 = dot(λt, ϕa.^2)
	b  = dot(λt, ϕb)
	b2 = dot(λt, ϕb.^2)
	ab = dot(λt, ϕa.*ϕb)

	var_a  = a2 - a^2
	var_b  = b2 - b^2
	cov_ab = ab - a*b

	prop_domestic = b/Bprime
	Bf = Bprime - b

	# print("\nvar_a, var_b, cov_ab = $([var_a, var_b, cov_ab])")

	μ′, σ′, q′, _ = compute_stats_logN(h, ζt, a, b, var_a, var_b, cov_ab, itp_qᵍ, Bprime, exp_rep)

	lμ = h.μgrid[end] - h.μgrid[1]
	lσ = h.σgrid[end] - h.σgrid[1]
	μ′ = max.(min.(μ′, h.μgrid[end]+0.0*lμ), h.μgrid[1]-0.0*lμ)
	σ′ = max.(min.(σ′, h.σgrid[end]+0.0*lσ), h.σgrid[1]-0.0*lσ)

	# μ′, σ′, q′ = new_expectations(h, itp_ϕa, itp_ϕb, itp_qᵍ, Bprime, wt, thres, Bt, μt, σt, w0, ζt, zt, jdef) # This would assume that λₜ is lognormal
	# print("\n$(q′)")

	# Draw z and the reentry shock for tomorrow, deduce ζ and correct B, μ, and σ as needed, and update the distribution
	probs = cumsum(h.Pz[jz,:])
	jzp = findfirst(probs .> rand())

	probs = cumsum(h.Pξ[jξ,:])
	jξp = findfirst(probs .> rand())

	zprime = h.zgrid[jzp]
	ξprime = h.ξgrid[jξp]

	μprime = μ′[jξp,jzp, 1]
	σprime = σ′[jξp,jzp, 1]
	qprime = q′[jξp,jzp, 1]

	if jdef
		# Compute welfare in case of reentry and remain in default
		Wr = itp_W(Bprime, μ′[jξp,jzp,1], σ′[jξp,jzp,1], ξt, 1, jzp)
		Wd = itp_W(Bprime, μ′[jξp,jzp,2], σ′[jξp,jzp,2], ξt, 2, jzp)
		# Now draw reentry
		prob_reentry = h.θ
		reentry = (rand() <= prob_reentry)
		if reentry
			ζprime = 1.0
			R = h.κ + (1.0-h.ρ) * qprime
		else
			ζprime = 2.0
			μprime = μ′[jξp,jzp, 2]
			σprime = σ′[jξp,jzp, 2]
			qprime = q′[jξp,jzp, 2]
			R = (1.0-h.ρ) * qprime
		end
	else
		# Compute welfare in case repay and default
		Wr = itp_W(Bprime, 			μ′[jξp,jzp,1], σ′[jξp,jzp,1], ξt, 1, jzp)
		Wd = itp_W((1.0.-h.ℏ)*Bprime, μ′[jξp,jzp,2], σ′[jξp,jzp,2], ξt, 2, jzp)
		# Now draw default
		if phase == "no def"
			repay_prime = 1.
		else
			repay_prime = (rand() <= exp_rep[jξp,jzp])
		end
		if repay_prime
			ζprime = 1.0
			R = h.κ + (1.0-h.ρ) * qprime
		else
			ζprime = 2.0
			Bprime = (1.0 - h.ℏ) * Bprime
			qprime = q′[jξp,jzp, 2]
			μprime = μ′[jξp,jzp, 2]
			σprime = σ′[jξp,jzp, 2]
			R = (1.0-h.ℏ)*(1.0-h.ρ) * qprime
		end
	end

	savings = ϕa + R*ϕb
	savings = max.(min.(savings, h.ωmax), h.ωmin)

	basis = Basis(LinParams(h.ωgrid_fine, 0))
	Qω = BasisMatrix(basis, Expanded(), savings, 0).vals[1]
	Q = row_kron(Qϵ, Qω)

	ζt == 1 && ζprime != 1 ? new_def = 1 : new_def = 0

	λprime = Q' * λt

	M, V = unmake_logN(μprime, σprime)

	# print("\n$(sum(λprime))")

	# Fill the path for next period
	if t < length(jz_series)
		jz_series[t+1] = jzp
		fill_path!(p,t+1, Dict(:B => Bprime, :μ => μprime, :σ => σprime, :ζ => ζprime, :ξ => ξprime, :z => zprime, :ψ => prop_domestic, :A => a, :Bh => b, :Bf => Bf, :Wr => Wr, :Wd => Wd, :qg => qprime, :mean => M, :var => V))
	end

	return λprime, new_def
end

function simul(h::Hank; simul_length::Int64=1, burn_in::Int64=0, only_def_end::Bool=false)

	Random.seed!(1)

	# Setup
	T = burn_in + simul_length
	p = Path(T = T)

	jz = 1

	B0, μ0, σ0, ξ0, ζ0, z0 = mean(h.bgrid), mean(h.μgrid), mean(h.σgrid), h.ξgrid[2], h.ζgrid[1], h.zgrid[jz]
	fill_path!(p,1, Dict(:B => B0, :μ => μ0, :σ => σ0, :w=>1.0, :ξ => ξ0, :ζ => ζ0, :z => z0))

	itp_ϕa = make_itp(h, h.ϕa_ext; agg=false)
	itp_ϕa = extrapolate(itp_ϕa, Interpolations.Flat())
	itp_ϕb = make_itp(h, h.ϕb_ext; agg=false)
	itp_ϕb = extrapolate(itp_ϕb, Interpolations.Flat())
	itp_ϕc = make_itp(h, h.ϕc_ext; agg=false)
	itp_ϕc = extrapolate(itp_ϕc, Interpolations.Flat())
	itp_vf = make_itp(h, h.vf; agg=false)
	itp_vf = extrapolate(itp_vf, Interpolations.Flat())


	itp_B′		= make_itp(h, h.issuance; agg=true)
	itp_G		= make_itp(h, h.spending; agg=true)
	itp_pN		= make_itp(h, h.pN; agg=true)
	itp_qᵍ 		= make_itp(h, h.qᵍ; agg=true)
	itp_W 		= make_itp(h, h.welfare; agg=true)

	rep_mat = reshape(h.repay, h.Nb, h.Nμ, h.Nσ, h.Nξ, h.Nζ, h.Nz, h.Nξ, h.Nz)
	knots = (h.bgrid, h.μgrid, h.σgrid, h.ξgrid, 1:h.Nζ, 1:h.Nz, h.ξgrid, 1:h.Nz)
	itp_repay = interpolate(knots, rep_mat, (Gridded(Linear()), Gridded(Linear()), Gridded(Linear()), Gridded(Linear()), NoInterp(), NoInterp(), Gridded(Linear()), NoInterp()))

	jz_series = Vector{Int64}(undef, T)
	jz_series[1] = jz

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
	Ndefs = 0
	for t in 1:T
		# Advance one step
		phase = ""
		if only_def_end
			start_def = floor(Int, burn_in + 3/4 * simul_length)
			if t - burn_in < start_def
				phase = "no_def"
				if t - burn_in > start_def - 4*3.5
					phase = "danger"
				end
			end
		end

		λ, new_def = iter_simul!(h, p, t, jz_series, itp_ϕa, itp_ϕb, itp_ϕc, itp_vf, itp_B′, itp_G, itp_pN, itp_qᵍ, itp_repay, itp_W, λ, Qϵ; phase = phase)

		Ndefs += new_def
	end

	jz_series = jz_series[burn_in+1:end]

	# Return stuff
	return p, jz_series, Ndefs
end

using DataFrames, GLM
function get_AR1(y::Vector)
	y_lag = y[1:end-1]
	y = y[2:end]

	data = DataFrame(yt = y, ylag = y_lag)
	OLS = glm(@formula(yt ~ ylag - 1), data, Normal(), IdentityLink())

	ρ = coef(OLS)[1]

	ϵ = y - ρ * y_lag

	σ = var(ϵ)^0.5

	return ρ, σ
end

get_MV(y::Vector) = mean(y), var(y)^0.5

function simul_stats(path::Path; nodef::Bool=false, ζ_vec::Vector=[])
	T = size(path.data, 1)
	
	if ζ_vec == []
		ζ_vec = series(path,:ζ) .- 1
	end

	conditional = (ζ_vec .> -Inf)
	if nodef
		conditional = (ζ_vec .== 1)
	end

	B_vec = series(path,:B)[conditional]
	μ_vec = series(path,:μ)[conditional]
	σ_vec = series(path,:σ)[conditional]
	w_vec = series(path,:w)[conditional]
	z_vec = exp.(series(path,:z)[conditional])
	Y_vec = series(path,:Y)[conditional]
	C_vec = series(path,:C)[conditional]
	G_vec = series(path,:G)[conditional]
	π_vec = series(path,:π)[conditional]
	ψ_vec = series(path,:ψ)[conditional]
	u_vec = 100.0 * (1.0 .- series(path, :L)[conditional])
	spr_vec = 1.0./series(path, :qg)[conditional] .- 1.0

	print("\nT = $T")

	m_vec, sd_vec = unmake_logN(μ_vec, σ_vec)
	mean_wealth = 100*mean(m_vec./(4*Y_vec))

	ρy, σy = get_AR1(log.(Y_vec.+1e-8))
	ρc, σc = get_AR1(log.(C_vec.+1e-8))
	if var(spr_vec) > 1e-6
		ρs, σs = get_AR1(spr_vec)
	else
		ρs, σs = 1.0, 0.0
	end
	m_unemp, sd_unemp = get_MV(u_vec)
	m_debt, sd_debt = get_MV(100*B_vec./(4*Y_vec))
	m_gspend, sd_gspend = get_MV(100*G_vec./Y_vec)
	ψ_mean = 100 * median(ψ_vec)
	spr_mean = mean(spr_vec)

	v_m = [ρy; σy; ρc; σc; ρs; σs; m_debt; sd_debt; m_unemp; sd_unemp; ψ_mean; mean_wealth]

	return v_m
end


function collect_episodes(path::Path, t_epi, N)

	sample = zeros(size(path.data)[2], 21, N)
	
	for jepi in 1:N	
		jt = t_epi[jepi]
		for jj in 0:20
			sample[:, jj+1, jepi] = vec(getfrompath(path, jt-10+jj))
		end
	end

	return sample
end

function find_times_episodes(path::Path; episode_type::String="default", πthres::Float64=0.975)
	ζ_vec = series(path,:ζ)
	π_vec = series(path,:π)
	π_thres = quantile(π_vec[π_vec.>0], πthres)
	ψ_vec = series(path, :ψ)
	ψ_thres = quantile(ψ_vec, 0.025)

	N = 0
	t_epi = []
	for jt in 20:length(ζ_vec) - 10
		if episode_type=="default"
			if ζ_vec[jt-1] == 1 && ζ_vec[jt] != 1
				N += 1
				push!(t_epi, jt)
			end
		elseif episode_type=="highspread"
			if minimum(π_vec[jt-2+1:jt+2]) >= π_thres
				N += 1
				push!(t_epi, jt)
			end
		elseif episode_type=="onlyspread"
			if minimum(π_vec[jt-2+1:jt+2]) >= π_thres && maximum(ζ_vec[jt-3+1:jt+3]) == 1# && ψ_vec[jt-3+1] >= ψ_thres
				N += 1
				push!(t_epi, jt)
			end
		else
			throw(error("Please select episode_type = 'default' or 'highspread' or 'onlyspread'"))
		end
	end

	if N == 0
		print_save("WARNING: No episodes of $(episode_type) found")
		# return sample
	else
		println("$N episodes found.")
	end
	# println(t_epi)
	return t_epi, N
end

function find_episodes(path::Path; episode_type::String="default", πthres::Float64=0.975)

	t_epi, N = find_times_episodes(path; episode_type=episode_type, πthres=πthres)

	sample = collect_episodes(path, t_epi, N)

	return sample, N
end

function IRF(p::Path, h::Hank; shock::String="z")
	if shock == "z"
		zvec = series(path, :z)
		Ez = zvec * h.ρz
	end
	t0 = 11
	T = length(zvec)

	eps_z = [zvec[jt] - Ez[jt-1] for jt in t0:T]

end