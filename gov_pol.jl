function integrate_itp(h::Hank, bv, μv, σv, ξv, jζ, jz, itp_obj)

	ωmin_int, ωmax_int = quantile.(LogNormal(μv, σv), [.0005; .9995]) + h.ωmin
	ωmax_int = min(ωmax_int, h.ωmax)
	W, sum_prob = 0., 0.
	for (jϵ, ϵv) in enumerate(h.ϵgrid)
		f_pdf(ω) = pdf(LogNormal(μv, σv), ω-h.ωmin)
		(val_pdf, err) = hquadrature(f_pdf, ωmin_int, ωmax_int, reltol=1e-10, abstol=1e-12, maxevals=0)
		sum_prob += val_pdf * h.λϵ[jϵ]

	    f(ω) = f_pdf(ω) * itp_obj[ω, jϵ, bv, μv, σv, ξv, jζ, jz]
	    (val, err) = hquadrature(f, ωmin_int, ωmax_int, reltol=1e-12, abstol=0, maxevals=0)
	    W += val * h.λϵ[jϵ]
	end
	W = W / sum_prob

	# val, sum_prob = 0., 0.
	# for (jϵ, ϵv) in enumerate(h.ϵgrid)
	# 	for jω = 1:length(h.ωgrid_fine)-1
	# 		ωv  = h.ωgrid_fine[jω]
	# 		ω1v = h.ωgrid_fine[jω+1]
	# 		ωmv = 0.5*(ωv+ω1v)

	# 		prob = pdf(LogNormal(μv, σv), ωmv-h.ωmin) * h.λϵ[jϵ] * (ω1v - ωv)

	# 		Y = itp_obj[ωmv, jϵ, bv, μv, σv, wv, jζ, jz]

	# 		val  += prob * Y
	# 		sum_prob += prob
	# 	end
	# end
	# W = val / sum_prob
	return W
end

function update_govpol(h::Hank; η_rep::Float64=0.5)
	itp_vf = make_itp(h, h.vf; agg=false)

	B′_mat = reshape(h.issuance, h.Nb, h.Nμ, h.Nσ, h.Nξ, h.Nζ, h.Nz)
	μ′_mat = reshape(h.μ′, h.Nb, h.Nμ, h.Nσ, h.Nξ, h.Nζ, h.Nz, h.Nξ, h.Nz, 2)
	σ′_mat = reshape(h.σ′, h.Nb, h.Nμ, h.Nσ, h.Nξ, h.Nζ, h.Nz, h.Nξ, h.Nz, 2)
	
	itp_W = make_itp(h, h.welfare; agg=true)

	μ_gov = 0.01 #* 0.0
	σ_gov = 0.0008

	repay = reshape(h.repay, h.Nb, h.Nμ, h.Nσ, h.Nξ, h.Nζ, h.Nz, h.Nξ, h.Nz)
	diff_W = Array{Float64}(h.Nb, h.Nμ, h.Nσ, h.Nξ, h.Nζ, h.Nz, h.Nξ, h.Nz)
	diff_R = Array{Float64}(h.Nb, h.Nμ, h.Nσ, h.Nξ, h.Nζ, h.Nz, h.Nξ, h.Nz)
	rep_prob = zeros(h.Nb, h.Nμ, h.Nσ, h.Nξ, h.Nζ, h.Nz, h.Nξ, h.Nz)
	for js in 1:size(h.Jgrid, 1)
		jb = h.Jgrid[js, 1]
		jμ = h.Jgrid[js, 2]
		jσ = h.Jgrid[js, 3]
		jξ = h.Jgrid[js, 4]
		jζ = h.Jgrid[js, 5]
		jz = h.Jgrid[js, 6]

		bvp = B′_mat[jb, jμ, jσ, jξ, jζ, jz]
		for (jξp, ξpv) in enumerate(h.ξgrid), jzp in 1:h.Nz
			if jζ == 1
				μvp = μ′_mat[jb, jμ, jσ, jξ, jζ, jz, jξp, jzp, 1]
				σvp = σ′_mat[jb, jμ, jσ, jξ, jζ, jz, jξp, jzp, 1]
				# Wr = integrate_itp(h, bvp, μvp, σvp, ξpv, 1, jzp, itp_vf)
				Wr = itp_W[bvp, μvp, σvp, ξpv, 1, jzp]

				μvp = μ′_mat[jb, jμ, jσ, jξ, jζ, jz, jξp, jzp, 2]
				σvp = σ′_mat[jb, jμ, jσ, jξ, jζ, jz, jξp, jzp, 2]
				# Wd = integrate_itp(h, (1.-h.ℏ)*bvp, μvp, σvp, ξpv, 2, jzp, itp_vf)
				Wd = itp_W[(1.-h.ℏ)*bvp, μvp, σvp, ξpv, 2, jzp]

				diff_W[jb, jμ, jσ, jξ, jζ, jz, jξp, jzp] = Wr - Wd
				if Wr > Wd && repay[jb, jμ, jσ, jξ, jζ, jz, jξp, jzp] < 0.5
					diff_R[jb, jμ, jσ, jξ, jζ, jz, jξp, jzp] = 1.
				elseif Wr < Wd && repay[jb, jμ, jσ, jξ, jζ, jz, jξp, jzp] > 0.5
					diff_R[jb, jμ, jσ, jξ, jζ, jz, jξp, jzp] = -1.
				end
				rep_prob[jb, jμ, jσ, jξ, jζ, jz, jξp, jzp] = 1.0 - cdf(Normal(μ_gov, σ_gov), Wd-Wr)
			else
				diff_W[jb, jμ, jσ, jξ, jζ, jz, jξp, jzp] = 0.
				diff_R[jb, jμ, jσ, jξ, jζ, jz, jξp, jzp] = 0.
				rep_prob[jb, jμ, jσ, jξ, jζ, jz, jξp, jzp] = 0.
			end
		end
	end

	if maximum(diff_R) == 1.
		threshold = quantile(diff_W[diff_R .== 1.], 1.-η_rep)
		ind_change = (diff_W .> threshold) .& (diff_R .== 1.)
		
		repay[ind_change] = 1.
	end

	if minimum(diff_R) == -1.
		threshold = quantile(diff_W[diff_R .== -1.], η_rep)
		ind_change = (diff_W .< threshold) .& (diff_R .== -1.)
		
		repay[ind_change] = 0.
	end
	
	# rep_new = reshape(repay, length(repay))
	rep_new = reshape(rep_prob, length(rep_prob))

	return rep_new
end

function update_W(h::Hank)
	itp_vf = make_itp(h, h.vf; agg=false)

	W = zeros(size(h.Jgrid, 1))
	for js in 1:size(h.Jgrid, 1)
		bv = h.bgrid[h.Jgrid[js, 1]]
		μv = h.μgrid[h.Jgrid[js, 2]]
		σv = h.σgrid[h.Jgrid[js, 3]]
		ξv = h.ξgrid[h.Jgrid[js, 4]]
		jζ = h.Jgrid[js, 5]
		jz = h.Jgrid[js, 6]
		
		W[js] = integrate_itp(h, bv, μv, σv, ξv, jζ, jz, itp_vf)
	end

	return W
end


function mpe_iter!(h::Hank; remote::Bool=false, maxiter::Int64=150, tol::Float64=1e-4, nodef::Bool=false, rep_agent::Bool=false, run_number::Int64=1)
	print_save("\nIterating on the government's policy: ")
	time_init = time()
	t_old = time_init
	out_iter = 1
	dist = 10.

	upd_η = 1.
	tol_vfi = 5e-2
	h.upd_tol = max(h.upd_tol, 1e-3)

	while dist > tol && out_iter < maxiter
		print_save("\n\nOuter Iteration $out_iter\n")
		vfi!(h, verbose = true, remote = remote, tol = tol_vfi, maxiter = 15)
		h.upd_tol = max(min(h.upd_tol*10, tol_vfi/10), 1e-5)
		
		W_new = update_W(h)

		h.welfare = upd_η * W_new + (1.-upd_η) * h.welfare
		upd_η = 0.5

		if isnan.(sum(h.welfare))
			print_save("\nWARNING: ||welf|| = NaN")
		end

		old_rep = copy(h.repay)

		if nodef
			h.repay = ones(h.repay)
			dist = 0.
		else
			new_rep = update_govpol(h; η_rep = 0.25)
			old_norm = sqrt.(sum(old_rep.^2))
			print_save("\n||rep₀|| = $(old_norm)")
			if isapprox(old_norm, 0.0)
				old_norm = 1.0
			end
			new_norm = sqrt.(sum(new_rep.^2))
			print_save("\n||repₜ|| = $(new_norm)")
			dist = sqrt.(sum( (new_rep - old_rep).^2 )) / old_norm
			h.repay = 0.4*upd_η * new_rep + (1.-0.4*upd_η) * old_rep
		end


		tol_vfi = max(exp(0.8*log(1+tol_vfi))-1, 1e-6)
		t_new = time()
		print_save("\n$(Dates.format(now(), "HH:MM")) Distance = $(@sprintf("%0.3g",dist)) after $(time_print(t_new-t_old)) and $out_iter iterations. New tol = $(@sprintf("%0.3g",tol_vfi))")

		dist = max(dist, tol_vfi)
		
		push!(h.outer_dists, dist)
		# plot_outerdists(h; remote = remote)


		if out_iter % 5 == 0 && out_iter > 5
			t_sim = time()
			print_save("\nSimulating")
			make_simulated_path(h, run_number)
			print_save(": done in $(time_print(time()-t_sim))")
		end

		out_iter += 1
		time_old = time()
	end
	if dist <= tol
		print_save("\nConverged in $out_iter iterations. ")
	else
		print_save("\nStopping at distance $(@sprintf("%0.3g",dist)). ")
	end

	print_save("\nTotal time: $(time_print(time()-time_init))\n")
end
