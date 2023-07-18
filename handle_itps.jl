function reshape_long(sd::SOEdef, y::Array)
    if length(y) == N(sd, :b) * N(sd, :μ) * N(sd, :σ) * N(sd, :ξ) * N(sd, :ζ) * N(sd, :z)
        y_mat = reshape(y, N(sd, :b), N(sd, :μ), N(sd, :σ), N(sd, :ξ), N(sd, :ζ), N(sd, :z))
    end
    return y_mat
end
function reshape_long_shocks(sd::SOEdef, y::Array)
    if length(y) == N(sd, :b) * N(sd, :μ) * N(sd, :σ) * N(sd, :ξ) * N(sd, :ζ) * N(sd, :z) * N(sd, :ξ) * N(sd, :z)
        y_mat = reshape(y, N(sd, :b), N(sd, :μ), N(sd, :σ), N(sd, :ξ), N(sd, :ζ), N(sd, :z), N(sd, :ξ), N(sd, :z))
    end
    return y_mat
end

flatten(y::Array) = reshape(y, length(y))

function make_itp(sd::SOEdef, y::Array; agg::Bool=false)
    if agg
        if length(size(y)) == 1
            length(y) == N(sd, :b) * N(sd, :μ) * N(sd, :σ) * N(sd, :ξ) * N(sd, :ζ) * N(sd, :z) || throw(error("wrong dimensions of interpolated aggregate vector"))
            y = reshape(y, N(sd, :b), N(sd, :μ), N(sd, :σ), N(sd, :ξ), N(sd, :ζ), N(sd, :z))
        end
        knots = (sd.gr[:b], sd.gr[:μ], sd.gr[:σ], sd.gr[:ξ], sd.gr[:ζ], sd.gr[:z])
    else
        knots = (sd.gr[:ω], sd.gr[:ϵ], sd.gr[:b], sd.gr[:μ], sd.gr[:σ], sd.gr[:ξ], sd.gr[:ζ], sd.gr[:z])
    end
    itp_obj = interpolate(knots, y, Gridded(Linear()))
    itp_obj = extrapolate(itp_obj, Interpolations.Line())
    return itp_obj

end

# N = dimensions of array; K = dimension of interpolation object
const Arr_itp_VF{N,K} = Array{Interpolations.GriddedInterpolation{Float64,K,Float64,Gridded{Linear},NTuple{K,Array{Float64,1}}},N}
const Arr_ext_vf{N,K} = Array{Interpolations.Extrapolation{Float64,K,Interpolations.GriddedInterpolation{Float64,K,Array{Float64,K},Gridded{Linear{Throw{OnGrid}}},NTuple{K,Array{Float64,1}}},Gridded{Linear{Throw{OnGrid}}},Line{Nothing}},N}

function integrate_itp(sd::SOEdef, bv, μv, σv, ξv, ζv, zv, itp_obj)
    pars, gr = sd.pars, sd.gr

    ωmin_int, ωmax_int = quantile.(LogNormal(μv, σv), [1e-6; 1 - 1e-6]) .+ pars[:ωmin]
    ωmax_int = min(ωmax_int, maximum(gr[:ω]))
    ωmin_int = pars[:ωmin]
    W, sum_prob = 0.0, 0.0
    # itp_obj = extrapolate(itp_obj, Interpolations.Line())

    λϵ = ergodic_ϵ(sd)

    for (jϵ, ϵv) in enumerate(gr[:ϵ])
        f_pdf(ω) = pdf(LogNormal(μv, σv), ω .- pars[:ωmin])
        (val_pdf, err) = hquadrature(f_pdf, ωmin_int, ωmax_int, rtol=1e-10, atol=1e-12, maxevals=0)
        sum_prob += val_pdf * λϵ[jϵ]

        f(ω) = f_pdf(ω) * itp_obj(ω, ϵv, bv, μv, σv, ξv, ζv, zv)
        (val, err) = hquadrature(f, ωmin_int, ωmax_int, rtol=1e-10, atol=1e-12, maxevals=0)
        W += val * λϵ[jϵ]
    end
    W = W / sum_prob

    return W
end