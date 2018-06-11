function welfare(h::Hank, bv, μv, σv, wv, jζ, jz, itp_vf)

    W = 0.
    for (jϵ, ϵv) in enumerate(h.ϵgrid)

        f(ω) = pdf(LogNormal(μv, σv), ω-h.ωmin) * h.λϵ[jϵ] * itp_vf[ω, jϵ, bv, μv, σv, wv, jζ, jz]

        (val, err) = hquadrature(f, h.ωmin, h.ωmax,
                                    reltol=1e-8, abstol=0, maxevals=0)

        W += val
    end
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
                Wr = welfare(h, bvp, μvp, σvp, wvp, 1, jzp, itp_vf)
                μvp = μ′_mat[jb, jμ, jσ, jw, jζ, jz, jzp, 2]
                σvp = σ′_mat[jb, jμ, jσ, jw, jζ, jz, jzp, 2]
                Wd = welfare(h, bvp, μvp, σvp, wvp, 2, jzp, itp_vf)
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
    iter = 1
    dist = 10.

    upd_η = 0.5

    print_save("\nOuter Iteration $iter")
	while dist > tol && iter < maxiter
        vfi!(h, verbose = true, remote = remote)

        old_rep = copy(h.repay)
        new_rep = update_govpol(h)

        dist = sqrt.(sum( (new_rep - old_rep).^2 )) / sqrt.(sum(old_rep.^2))
        h.repay = upd_η * new_rep + (1.-upd_η) * old_rep

        t_new = time()
        print_save("\nDistance = $(@sprintf("%0.3g",dist)) after $(time_print(t_new-t_old)) and $iter iterations ")
        print_save(Dates.format(now(), "HH:MM"))

        iter += 1
        time_old = time()
    end
    if dist <= tol
		print_save("\nConverged in $iter iterations. ")
	else
		print_save("\nStopping at distance $(@sprintf("%0.3g",dist)). ")
	end

	print_save("\nTotal time: $(time_print(time()-time_init))\n")
end
