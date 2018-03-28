function make_itps(h::Hank, Y; agg::Bool=true)
	ωrange = linspace(h.ωgrid[1], h.ωgrid[end], h.Nω)
	brange = linspace(h.bgrid[1], h.bgrid[end], h.Nb)
	μrange = linspace(h.μgrid[1], h.μgrid[end], h.Nμ)
	σrange = linspace(h.σgrid[1], h.σgrid[end], h.Nσ)
	wrange = linspace(h.wgrid[1], h.wgrid[end], h.Nw)
	zrange = linspace(h.zgrid[1], h.zgrid[end], h.Nz)

	ext = false
	if length(size(Y)) == 9
		agg = false
		ext = true
	end

	if agg
		Y_mat = reshape(Y, h.Nb, h.Nμ, h.Nσ, h.Nw, h.Nζ, h.Nz)

		unscaled_itp  = interpolate(Y_mat, (BSpline(Linear()), BSpline(Linear()), BSpline(Linear()), BSpline(Linear()), NoInterp(), BSpline(Linear())), OnGrid())
		itp_obj = Interpolations.scale(unscaled_itp, brange, μrange, σrange, wrange, 1:h.Nζ, zrange)
	elseif ext
		pnrange = linspace(h.pngrid[1], h.pngrid[end], length(h.pngrid))
		unscaled_itp = interpolate(Y, (BSpline(Quadratic(Line())), NoInterp(), BSpline(Linear()), BSpline(Linear()), BSpline(Linear()), BSpline(Linear()), NoInterp(), BSpline(Linear()), BSpline(Linear())), OnGrid())
		itp_obj = Interpolations.scale(unscaled_itp, ωrange, 1:h.Nϵ, brange, μrange, σrange, wrange, 1:h.Nζ, zrange, pnrange)
	else
		unscaled_itp = interpolate(Y, (BSpline(Quadratic(Line())), NoInterp(), BSpline(Linear()), BSpline(Linear()), BSpline(Linear()), BSpline(Linear()), NoInterp(), BSpline(Linear())), OnGrid())
		itp_obj = Interpolations.scale(unscaled_itp, ωrange, 1:h.Nϵ, brange, μrange, σrange, wrange, 1:h.Nζ, zrange)
	end

	return itp_obj
end