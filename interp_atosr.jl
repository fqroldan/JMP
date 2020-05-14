const Arr_itp_VF = Array{Interpolations.GriddedInterpolation{Float64,2,Float64,Tuple{Gridded{Linear},NoInterp},Tuple{Array{Float64,1},UnitRange{Int64}}},3}

function make_itp(h::Hank, Y; agg::Bool=true)
	# ωrange = range(h.ωgrid[1], h.ωgrid[end], length=h.Nω)
	# brange = range(h.bgrid[1], h.bgrid[end], length=h.Nb)
	# μrange = range(h.μgrid[1], h.μgrid[end], length=h.Nμ)
	# σrange = range(h.σgrid[1], h.σgrid[end], length=h.Nσ)
	# ξrange = range(h.ξgrid[1], h.ξgrid[end], length=h.Nξ)
	# zrange = range(h.zgrid[1], h.zgrid[end], length=h.Nz)

	interp_gridded = true

	ext = false
	if length(size(Y)) == 9
		agg = false
		ext = true
	end

	if agg
		if length(size(Y)) == 1
			length(Y) == h.Nb * h.Nμ * h.Nσ * h.Nξ * h.Nζ * h.Nz || throw(error("wrong dimensions of interpolated aggregate vector"))
			Y = reshape(Y, h.Nb, h.Nμ, h.Nσ, h.Nξ, h.Nζ, h.Nz)
		end

		if interp_gridded
			knots = (h.bgrid, h.μgrid, h.σgrid, h.ξgrid, 1:h.Nζ, 1:h.Nz)
			itp_obj = interpolate(knots, Y, (Gridded(Linear()), Gridded(Linear()), Gridded(Linear()), Gridded(Linear()), NoInterp(), NoInterp()))
		# else
		# 	unscaled_itp  = interpolate(Y, (BSpline(Linear()), BSpline(Linear()), BSpline(Linear()), BSpline(Linear()), NoInterp(), NoInterp()), OnGrid())
		# 	itp_obj = Interpolations.scale(unscaled_itp, brange, μrange, σrange, ξrange, 1:h.Nζ, 1:h.Nz)
		end
	# elseif ext
	# 	if interp_gridded
	# 		knots = (h.ωgrid, 1:h.Nϵ, h.bgrid, h.μgrid, h.σgrid, h.ξgrid, 1:h.Nζ, 1:h.Nz, h.pngrid)
	# 		itp_obj = interpolate(knots, Y, (Gridded(Linear()), NoInterp(), Gridded(Linear()), Gridded(Linear()), Gridded(Linear()), Gridded(Linear()), NoInterp(), NoInterp(), Gridded(Linear())))
	# 	else
	# 		pnrange = linspace(h.pngrid[1], h.pngrid[end], length(h.pngrid))
	# 		unscaled_itp = interpolate(Y, (BSpline(Quadratic(Line())), NoInterp(), BSpline(Linear()), BSpline(Linear()), BSpline(Linear()), BSpline(Linear()), NoInterp(), NoInterp(), BSpline(Linear())), OnGrid())
	# 		itp_obj = Interpolations.scale(unscaled_itp, ωrange, 1:h.Nϵ, brange, μrange, σrange, ξrange, 1:h.Nζ, 1:h.Nz, pnrange)
	# 	end
	else
		if interp_gridded
			knots = (h.ωgrid, 1:h.Nϵ, h.bgrid, h.μgrid, h.σgrid, h.ξgrid, 1:h.Nζ, 1:h.Nz)
			itp_obj = interpolate(knots, Y, (Gridded(Linear()), NoInterp(), Gridded(Linear()), Gridded(Linear()), Gridded(Linear()), Gridded(Linear()), NoInterp(), NoInterp()))
		# else
		# 	unscaled_itp = interpolate(Y, (BSpline(Quadratic(Line())), NoInterp(), BSpline(Linear()), BSpline(Linear()), BSpline(Linear()), BSpline(Linear()), NoInterp(), NoInterp()), OnGrid())
		# 	itp_obj = Interpolations.scale(unscaled_itp, ωrange, 1:h.Nϵ, brange, μrange, σrange, ξrange, 1:h.Nζ, 1:h.Nz)
		end
	end

	return itp_obj
end
