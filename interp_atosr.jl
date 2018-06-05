const Arr_itp_VF =  Array{Interpolations.ScaledInterpolation{Float64,2,Interpolations.BSplineInterpolation{Float64,2,Array{Float64,2},Tuple{Interpolations.BSpline{Interpolations.Quadratic{Interpolations.Line}},Interpolations.NoInterp},Interpolations.OnGrid,(1, 0)},Tuple{Interpolations.BSpline{Interpolations.Quadratic{Interpolations.Line}},Interpolations.NoInterp},Interpolations.OnGrid,Tuple{StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}},UnitRange{Int64}}}, 2}

# const Arr_itp_VF = Array{Interpolations.GriddedInterpolation{Float64,2,Float64,Tuple{Interpolations.Gridded{Interpolations.Linear},Interpolations.NoInterp},Tuple{Array{Float64,1},Array{Int64,1}},0}}

function make_itp(h::Hank, Y; agg::Bool=true)
	ωrange = linspace(h.ωgrid[1], h.ωgrid[end], h.Nω)
	brange = linspace(h.bgrid[1], h.bgrid[end], h.Nb)
	μrange = linspace(h.μgrid[1], h.μgrid[end], h.Nμ)
	σrange = linspace(h.σgrid[1], h.σgrid[end], h.Nσ)
	wrange = linspace(h.wgrid[1], h.wgrid[end], h.Nw)
	# zrange = linspace(h.zgrid[1], h.zgrid[end], h.Nz)

	ext = false
	if length(size(Y)) == 9
		agg = false
		ext = true
	end

	if agg
		if length(size(Y)) == 1
			length(Y) == h.Nb * h.Nμ * h.Nσ * h.Nw * h.Nζ * h.Nz || throw(error("wrong dimensions of interpolated aggregate vector"))
			Y = reshape(Y, h.Nb, h.Nμ, h.Nσ, h.Nw, h.Nζ, h.Nz)
		end

		# knots = (h.bgrid, h.μgrid, h.σgrid, h.wgrid, 1:h.Nζ, 1:h.Nz)
		# itp_obj = interpolate(knots, Y, (Gridded(Linear()), Gridded(Linear()), Gridded(Linear()), Gridded(Linear()), NoInterp(), NoInterp()))
		unscaled_itp  = interpolate(Y, (BSpline(Linear()), BSpline(Linear()), BSpline(Linear()), BSpline(Linear()), NoInterp(), NoInterp()), OnGrid())
		itp_obj = Interpolations.scale(unscaled_itp, brange, μrange, σrange, wrange, 1:h.Nζ, 1:h.Nz)
	elseif ext
		# knots = (h.bgrid, h.μgrid, h.σgrid, h.wgrid, 1:h.Nζ, 1:h.Nz, h.pngrid)
		# itp_obj = interpolate(knots, Y, (Gridded(Linear()), Gridded(Linear()), Gridded(Linear()), Gridded(Linear()), NoInterp(), NoInterp(), Gridded(Linear())))
		pnrange = linspace(h.pngrid[1], h.pngrid[end], length(h.pngrid))
		unscaled_itp = interpolate(Y, (BSpline(Quadratic(Line())), NoInterp(), BSpline(Linear()), BSpline(Linear()), BSpline(Linear()), BSpline(Linear()), NoInterp(), NoInterp(), BSpline(Linear())), OnGrid())
		itp_obj = Interpolations.scale(unscaled_itp, ωrange, 1:h.Nϵ, brange, μrange, σrange, wrange, 1:h.Nζ, 1:h.Nz, pnrange)
	else
		# knots = (h.ωgrid, 1:h.Nϵ, h.bgrid, h.μgrid, h.σgrid, h.wgrid, 1:h.Nζ, 1:h.Nz)
		# itp_obj = interpolate(knots, Y, (Gridded(Linear()), NoInterp(), Gridded(Linear()), Gridded(Linear()), Gridded(Linear()), Gridded(Linear()), NoInterp(), NoInterp()))
		unscaled_itp = interpolate(Y, (BSpline(Quadratic(Line())), NoInterp(), BSpline(Linear()), BSpline(Linear()), BSpline(Linear()), BSpline(Linear()), NoInterp(), NoInterp()), OnGrid())
		itp_obj = Interpolations.scale(unscaled_itp, ωrange, 1:h.Nϵ, brange, μrange, σrange, wrange, 1:h.Nζ, 1:h.Nz)
	end

	return itp_obj
end
