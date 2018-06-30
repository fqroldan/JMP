function integrate_itp(h::Hank, bv, μv, σv, wv, jζ, jz, itp_obj)

    # ωmin_int, ωmax_int = quantile.(LogNormal(μv, σv), [.0001; .9999]) + h.ωmin
    # W = 0.
    # for (jϵ, ϵv) in enumerate(h.ϵgrid)
    #     f(ω) = pdf(LogNormal(μv, σv), ω-h.ωmin) * h.λϵ[jϵ] * itp_obj[ω, jϵ, bv, μv, σv, wv, jζ, jz]
    #     (val, err) = hquadrature(f, ωmin_int, ωmax_int, reltol=1e-10, abstol=0, maxevals=0)
    #     W += val
    # end
    val, sum_prob = 0., 0.
    for (jϵ, ϵv) in enumerate(h.ϵgrid)
		for jω = 1:length(h.ωgrid_fine)-1
			ωv  = h.ωgrid_fine[jω]
			ω1v = h.ωgrid_fine[jω+1]
			ωmv = 0.5*(ωv+ω1v)

			prob = pdf(LogNormal(μv, σv), ωmv-h.ωmin) * h.λϵ[jϵ] * (ω1v - ωv)

			Y = itp_obj[ωmv, jϵ, bv, μv, σv, wv, jζ, jz]

			val  += prob * ϕa
			sum_prob += prob
		end
	end
    W = val / sum_prob
    return W
end

function update_govpol(h::Hank)
    itp_vf = make_itp(h, h.vf; agg=false)

    B′_mat = reshape(h.issuance, h.Nb, h.Nμ, h.Nσ, h.Nw, h.Nζ, h.Nz)
    μ′_mat = reshape(h.μ′, h.Nb, h.Nμ, h.Nσ, h.Nw, h.Nζ, h.Nz, h.Nz, 2)
    σ′_mat = reshape(h.σ′, h.Nb, h.Nμ, h.Nσ, h.Nw, h.Nζ, h.Nz, h.Nz, 2)
    w′_mat = reshape(h.wage, h.Nb, h.Nμ, h.Nσ, h.Nw, h.Nζ, h.Nz)

    rep_new = zeros(h.repay)

    jsz = 0
    for js in 1:size(h.Jgrid, 1)
        jb = h.Jgrid[js, 1]
        jμ = h.Jgrid[js, 2]
        jσ = h.Jgrid[js, 3]
        jw = h.Jgrid[js, 4]
        jζ = h.Jgrid[js, 5]
        jz = h.Jgrid[js, 6]

        bvp = B′_mat[jb, jμ, jσ, jw, jζ, jz]
        wvp = w′_mat[jb, jμ, jσ, jw, jζ, jz]
        for jzp in 1:h.Nz
            jsz += 1
            if jζ == 1
                μvp = μ′_mat[jb, jμ, jσ, jw, jζ, jz, jzp, 1]
                σvp = σ′_mat[jb, jμ, jσ, jw, jζ, jz, jzp, 1]
                Wr = integrate_itp(h, bvp, μvp, σvp, wvp, 1, jzp, itp_vf)
                μvp = μ′_mat[jb, jμ, jσ, jw, jζ, jz, jzp, 2]
                σvp = σ′_mat[jb, jμ, jσ, jw, jζ, jz, jzp, 2]
                Wd = integrate_itp(h, (1.-h.ℏ)*bvp, μvp, σvp, wvp, 2, jzp, itp_vf)
                if Wr > Wd
                    rep_new[jsz] = 1.
                else
                    rep_new[jsz] = 0.
                end
            else
                rep_new[jsz] = 1.
            end
        end
    end
    return rep_new
end

function mpe_iter!(h::Hank; remote::Bool=false, maxiter::Int64=100, tol::Float64=1e-4)
    print_save("\nIterating on the government's policy: ")
    time_init = time()
    t_old = time_init
    out_iter = 1
    dist = 10.

    upd_η = 0.33
    tol_vfi = 2e-2

	while dist > tol && out_iter < maxiter
        print_save("\n\nOuter Iteration $out_iter\n")
        vfi!(h, verbose = true, remote = remote, tol = tol_vfi, maxiter = 40)
        h.upd_tol = 1e-3

        old_rep = copy(h.repay)
        new_rep = update_govpol(h)

        dist = sqrt.(sum( (new_rep - old_rep).^2 )) / sqrt.(sum(old_rep.^2))
        h.repay = upd_η * new_rep + (1.-upd_η) * old_rep

        tol_vfi = max(exp(0.9*log(1+tol_vfi))-1, 1e-6)
        t_new = time()
        print_save("\n$(Dates.format(now(), "HH:MM")) Distance = $(@sprintf("%0.3g",dist)) after $(time_print(t_new-t_old)) and $out_iter iterations. New tol = $(@sprintf("%0.3g",tol_vfi))")

        push!(h.outer_dists, dist)
        plot_outerdists(h; remote = remote)

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
