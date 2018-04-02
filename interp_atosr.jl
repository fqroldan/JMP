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

		unscaled_itp  = interpolate(Y, (BSpline(Linear()), BSpline(Linear()), BSpline(Linear()), BSpline(Linear()), NoInterp(), NoInterp()), OnGrid())
		itp_obj = Interpolations.scale(unscaled_itp, brange, μrange, σrange, wrange, 1:h.Nζ, 1:h.Nz)
	elseif ext
		pnrange = linspace(h.pngrid[1], h.pngrid[end], length(h.pngrid))
		unscaled_itp = interpolate(Y, (BSpline(Quadratic(Line())), NoInterp(), BSpline(Linear()), BSpline(Linear()), BSpline(Linear()), BSpline(Linear()), NoInterp(), NoInterp(), BSpline(Linear())), OnGrid())
		itp_obj = Interpolations.scale(unscaled_itp, ωrange, 1:h.Nϵ, brange, μrange, σrange, wrange, 1:h.Nζ, 1:h.Nz, pnrange)
	else
		unscaled_itp = interpolate(Y, (BSpline(Quadratic(Line())), NoInterp(), BSpline(Linear()), BSpline(Linear()), BSpline(Linear()), BSpline(Linear()), NoInterp(), NoInterp()), OnGrid())
		itp_obj = Interpolations.scale(unscaled_itp, ωrange, 1:h.Nϵ, brange, μrange, σrange, wrange, 1:h.Nζ, 1:h.Nz)
	end

	return itp_obj
end