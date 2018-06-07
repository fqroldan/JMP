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
