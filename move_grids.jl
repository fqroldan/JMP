function update_grids(h::Hank, new_μgrid, new_σgrid)

	all_knots = (h.ωgrid, h.ϵgrid, h.bgrid, h.μgrid, h.σgrid, h.wgrid, h.ζgrid, h.zgrid)
	agg_knots = (h.bgrid, h.μgrid, h.σgrid, h.wgrid, h.ζgrid, h.zgrid)

	itp_ϕa			= interpolate(all_knots, ϕa, Gridded(Linear()))
	itp_ϕb			= interpolate(all_knots, ϕb, Gridded(Linear()))
	itp_ϕc			= interpolate(all_knots, ϕc, Gridded(Linear()))
	itp_vf			= interpolate(all_knots, vf, Gridded(Linear()))

	itp_Ld			= interpolate(agg_knots, Ld, Gridded(Linear()))
	itp_wage		= interpolate(agg_knots, wage, Gridded(Linear()))

	itp_repay		= interpolate(agg_knots, repay, Gridded(Linear()))
	itp_issuance	= interpolate(agg_knots, issuance, Gridded(Linear()))
	itp_spending	= interpolate(agg_knots, spending, Gridded(Linear()))

	itp_pN			= interpolate(agg_knots, pN, Gridded(Linear()))

	itp_μ′			= interpolate(agg_knots, μ′, Gridded(Linear()))
	itp_σ′			= interpolate(agg_knots, σ′, Gridded(Linear()))
	itp_w′			= interpolate(agg_knots, w′, Gridded(Linear()))

	h.ϕa			= itp_ϕa[h.ωgrid, h.ϵgrid, h.bgrid, new_μgrid, new_σgrid, h.wgrid, h.ζgrid, h.zgrid]
	h.ϕb			= itp_ϕb[h.ωgrid, h.ϵgrid, h.bgrid, new_μgrid, new_σgrid, h.wgrid, h.ζgrid, h.zgrid]
	h.ϕc			= itp_ϕc[h.ωgrid, h.ϵgrid, h.bgrid, new_μgrid, new_σgrid, h.wgrid, h.ζgrid, h.zgrid]
	h.vf			= itp_vf[h.ωgrid, h.ϵgrid, h.bgrid, new_μgrid, new_σgrid, h.wgrid, h.ζgrid, h.zgrid]
	h.Ld			= itp_Ld[h.bgrid, new_μgrid, new_σgrid, h.wgrid, h.ζgrid, h.zgrid]
	h.wage			= itp_wage[h.bgrid, new_μgrid, new_σgrid, h.wgrid, h.ζgrid, h.zgrid]
	h.repay			= itp_repay[h.bgrid, new_μgrid, new_σgrid, h.wgrid, h.ζgrid, h.zgrid]
	h.issuance		= itp_issuance[h.bgrid, new_μgrid, new_σgrid, h.wgrid, h.ζgrid, h.zgrid]
	h.spending		= itp_spending[h.bgrid, new_μgrid, new_σgrid, h.wgrid, h.ζgrid, h.zgrid]
	h.pN			= itp_pN[h.bgrid, new_μgrid, new_σgrid, h.wgrid, h.ζgrid, h.zgrid]
	h.μ′			= itp_μ′[h.bgrid, new_μgrid, new_σgrid, h.wgrid, h.ζgrid, h.zgrid]
	h.σ′			= itp_σ′[h.bgrid, new_μgrid, new_σgrid, h.wgrid, h.ζgrid, h.zgrid]
	h.w′			= itp_w′[h.bgrid, new_μgrid, new_σgrid, h.wgrid, h.ζgrid, h.zgrid]

	Void
end