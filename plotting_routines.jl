using Interpolations, PlotlyJS

col = [	"#1f77b4",  # muted blue
		"#ff7f0e",  # safety orange
		"#2ca02c",  # cooked asparagus green
		"#d62728",  # brick red
		"#9467bd",  # muted purple
		"#8c564b",  # chestnut brown
		"#e377c2",  # raspberry yogurt pink
		"#7f7f7f",  # middle gray
		"#bcbd22",  # curry yellow-green
		"#17becf"   # blue-teal
		]

function plot_hh_policies(h::Hank; remote::Bool=false)
	# leg = Array{LaTeXStrings.LaTeXString}(1, h.Nϵ)
	leg = Array{String}(1, h.Nϵ)
	for jϵ in 1:h.Nϵ
		# leg[jϵ] = latexstring("\\epsilon = $(round(h.ϵgrid[jϵ],2))")
		leg[jϵ] = "ϵ = $(round(h.ϵgrid[jϵ],2))"
	end

	show_b, show_μ, show_σ, show_w, show_ζ, show_z = mean(h.bgrid), mean(h.μgrid), mean(h.σgrid), mean(h.wgrid), h.ζgrid[1], h.zgrid[end]

	function hh_pol(h::Hank, show_b, show_μ, show_σ, show_w, show_ζ, show_z)
		knots = (h.ωgrid, h.ϵgrid, h.bgrid, h.μgrid, h.σgrid, h.wgrid, h.ζgrid, h.zgrid)
		itp_ϕa  = interpolate(knots, h.ϕa, Gridded(Linear()))
		itp_ϕb  = interpolate(knots, h.ϕb, Gridded(Linear()))
		itp_ϕc  = interpolate(knots, h.ϕc, Gridded(Linear()))
		itp_vf  = interpolate(knots, h.vf, Gridded(Linear()))

		qᵍ_mat  = reshape(h.qᵍ, h.Nb, h.Nμ, h.Nσ, h.Nw, h.Nζ, h.Nz)
		pN_mat  = reshape(h.pN, h.Nb, h.Nμ, h.Nσ, h.Nw, h.Nζ, h.Nz)
		agg_knots = (h.bgrid, h.μgrid, h.σgrid, h.wgrid, h.ζgrid, h.zgrid)
		itp_qᵍ  = interpolate(agg_knots, qᵍ_mat, Gridded(Linear()))
		itp_pN  = interpolate(agg_knots, pN_mat, Gridded(Linear()))


		ϕc_mat = itp_ϕc[h.ωgrid, h.ϵgrid, show_b, show_μ, show_σ, show_w, show_ζ, show_z]
		ϕa_mat = itp_ϕa[h.ωgrid, h.ϵgrid, show_b, show_μ, show_σ, show_w, show_ζ, show_z]
		ϕb_mat = itp_ϕb[h.ωgrid, h.ϵgrid, show_b, show_μ, show_σ, show_w, show_ζ, show_z]
		vf_mat = itp_vf[h.ωgrid, h.ϵgrid, show_b, show_μ, show_σ, show_w, show_ζ, show_z]
		qᵍ_mat = itp_qᵍ[show_b, show_μ, show_σ, show_w, show_ζ, show_z]

		qᵍ_all = zeros(vf_mat)
		for jω in 1:h.Nω
			for jϵ in 1:h.Nϵ
				qᵍ_all[jω, jϵ, :,:,:,:,:,:] = qᵍ_mat
			end
		end

		ωg_mat = 1.0/(1.0+h.r_star) * ϕa_mat + qᵍ_all .* ϕb_mat
		θg_mat = 1.0/(1.0+h.r_star) * (ϕa_mat - h.ωmin) ./ (ωg_mat - 1.0/(1.0+h.r_star)*h.ωmin)
		θg_mat[isapprox.(ωg_mat, h.ωmin)] = 1.0

		l = Array{PlotlyBase.GenericTrace{Dict{Symbol,Any}}}(h.Nϵ, 4)
		for (jϵ, ϵv) in enumerate(h.ϵgrid)
			l_new = scatter(;x=h.ωgrid, y=ϕc_mat[:,jϵ,1,1,1,1,1,1], line_shape="spline", name = "ϵ = $(round(exp(ϵv),4))", showlegend=false, marker_color=col[jϵ])
			l[jϵ,1] = l_new
			l_new = scatter(;x=h.ωgrid, y=vf_mat[:,jϵ,1,1,1,1,1,1], line_shape="spline", name = "ϵ = $(round(exp(ϵv),4))", showlegend=false, marker_color=col[jϵ])
			l[jϵ,2] = l_new
			l_new = scatter(;x=h.ωgrid, y=ωg_mat[:,jϵ,1,1,1,1,1,1], showlegend=false, name = "ϵ = $(round(exp(ϵv),4))", marker_color=col[jϵ])
			l[jϵ,3] = l_new
			# l_new = scatter(;x=h.ωgrid, y=ϕb_mat[:,jϵ,1,1,1,1,1,1], showlegend=false, marker_color=col[jϵ])
			l_new = scatter(;x=h.ωgrid, y=θg_mat[:,jϵ,1,1,1,1,1,1], showlegend=false, name = "ϵ = $(round(exp(ϵv),4))", marker_color=col[jϵ])
			l[jϵ,4] = l_new
		end

		ωmax_show = min(h.ωmax, quantile(LogNormal(show_μ, show_σ), 0.999)+h.ωmin)

		pc = plot([l[jϵ, 1] for jϵ in 1:h.Nϵ], Layout(; xaxis=attr(title="ω", zeroline=true, range=[h.ωmin, ωmax_show]), font_size=16, title="Consumption"))
		pv = plot([l[jϵ, 2] for jϵ in 1:h.Nϵ], Layout(; xaxis=attr(title="ω", zeroline=true, range=[h.ωmin, ωmax_show]), font_size=16, title="Value function"))
		pb = plot([l[jϵ, 3] for jϵ in 1:h.Nϵ], Layout(; xaxis=attr(title="ω", zeroline=true, range=[h.ωmin, ωmax_show]), font_size=16, title="Savings"))
		pθ = plot([l[jϵ, 4] for jϵ in 1:h.Nϵ], Layout(; xaxis=attr(title="ω", zeroline=true, range=[h.ωmin, ωmax_show]), font_size=16, title="Proportion risk-free debt"))

		p = [pc pv; pb pθ]
		p.plot.layout["xlabel"] = "ω"
		p.plot.layout["width"] = 800
		p.plot.layout["height"] = 600
		p.plot.layout["font_family"] = "Fira Sans Light"

		return p
	end

	p = hh_pol(h, show_b, show_μ, show_σ, show_w, show_ζ, show_z)

	if remote
		path = pwd() * "/../../Graphs/"
		save(path * "p_hh.jld", "p", p)
	else
		savefig(p, pwd() * "/../Graphs/hh.pdf")
	end

	show_b, show_μ, show_σ, show_w, show_ζ, show_z = mean(h.bgrid), mean(h.μgrid), mean(h.σgrid), mean(h.wgrid), h.ζgrid[2], h.zgrid[1]

	p = hh_pol(h, show_b, show_μ, show_σ, show_w, show_ζ, show_z)

	if remote
		path = pwd() * "/../../Graphs/"
		save(path * "p_hh_def.jld", "p", p)
	else
		path = pwd() * "/../Graphs/"
		savefig(p, path * "hh_def.pdf")
	end

	return Void
end

function plot_hh_policies_z(h::Hank; remote::Bool=false)
	show_ϵ, show_b, show_μ, show_σ, show_w, show_ζ = mean(h.ϵgrid), mean(h.bgrid), mean(h.μgrid), mean(h.σgrid), mean(h.wgrid), h.ζgrid[1]

	knots = (h.ωgrid, h.ϵgrid, h.bgrid, h.μgrid, h.σgrid, h.wgrid, h.ζgrid, h.zgrid)
	itp_ϕc  = interpolate(knots, h.ϕc, Gridded(Linear()))
	itp_vf  = interpolate(knots, h.vf, Gridded(Linear()))
	knots = (h.ωgrid, h.ϵgrid, h.bgrid, h.μgrid, h.σgrid, h.wgrid, h.ζgrid, h.zgrid, h.pngrid)
	itp_ϕc_ext  = interpolate(knots, h.ϕc_ext, Gridded(Linear()))

	itp_pN = make_itp(h, h.pN, agg=true)


	l = Array{PlotlyBase.GenericTrace{Dict{Symbol,Any}}}(h.Nz, 4)
	Cz = Vector{Float64}(h.Nz)
	Cz_fix = Vector{Float64}(h.Nz)
	for (jz, zv) in enumerate(h.zgrid)
		ϕc_vec = zeros(h.Nω)
		ϕce_vec = zeros(h.Nω)
		ϕce_vecfix = zeros(h.Nω)
		vf_vec = zeros(h.Nω)
		show_pN = itp_pN[show_b, show_μ, show_σ, show_w, 1., jz]
		for (jω, ωv) in enumerate(h.ωgrid)
			ϕc_vec[jω] = itp_ϕc[ωv, show_ϵ, show_b, show_μ, show_σ, show_w, show_ζ, zv]
			ϕce_vec[jω] = itp_ϕc_ext[ωv, show_ϵ, show_b, show_μ, show_σ, show_w, show_ζ, zv, show_pN]
			ϕce_vecfix[jω] = itp_ϕc_ext[ωv, show_ϵ, show_b, show_μ, show_σ, show_w, show_ζ, zv, mean(h.pngrid)]
			vf_vec[jω] = itp_vf[ωv, show_ϵ, show_b, show_μ, show_σ, show_w, show_ζ, zv]
		end
		l_new = scatter(;x=h.ωgrid, y=ϕc_vec, line_shape="spline", name="z = $(round(exp(zv),2))", marker_color=col[ceil(Int,10*jz/h.Nz)])
		l[jz,1] = l_new
		l_new = scatter(;x=h.ωgrid, y=ϕce_vec, line_shape="spline", name="z = $(round(exp(zv),2))", showlegend=false, marker_color=col[ceil(Int,10*jz/h.Nz)])
		l[jz,2] = l_new
		l_new = scatter(;x=h.ωgrid, y=ϕce_vecfix, line_shape="spline", name="z = $(round(exp(zv),2))", showlegend=false, marker_color=col[ceil(Int,10*jz/h.Nz)])
		l[jz,3] = l_new
		l_new = scatter(;x=h.ωgrid, y=vf_vec, line_shape="spline", name="z = $(round(exp(zv),2))", showlegend=false, marker_color=col[ceil(Int,10*jz/h.Nz)])
		l[jz,4] = l_new

		ωmin_int, ωmax_int = quantile.(LogNormal(show_μ, show_σ), [.005; .995]) + h.ωmin
		val_int_C, val_int_Cfix = 0., 0.
		for (jϵ, ϵv) in enumerate(h.ϵgrid)
			f(ω) = pdf(LogNormal(show_μ, show_σ), ω-h.ωmin) * itp_ϕc_ext[ω, ϵv, show_b, show_μ, show_σ, show_w, show_ζ, zv, show_pN]
			ffix(ω) = pdf(LogNormal(show_μ, show_σ), ω-h.ωmin) * itp_ϕc_ext[ω, ϵv, show_b, show_μ, show_σ, show_w, show_ζ, zv, mean(h.pngrid)]
			(val, err) = hquadrature(f, ωmin_int, ωmax_int, reltol=1e-12, abstol=0, maxevals=0)
			(valfix, err) = hquadrature(ffix, ωmin_int, ωmax_int, reltol=1e-12, abstol=0, maxevals=0)
			val_int_C += val * h.λϵ[jϵ] 
			val_int_Cfix += valfix * h.λϵ[jϵ] 
		end

		Cz[jz], Cz_fix[jz] = val_int_C, val_int_Cfix
	end

	ωmax_show = min(h.ωmax, quantile(LogNormal(show_μ, show_σ), 0.999)+h.ωmin)

	pc = plot([l[jz, 1] for jz in 1:h.Nz], Layout(; xaxis=attr(title="ω", zeroline=true, range=[h.ωmin, ωmax_show]), font_size=16, title="Consumption"))
	pce = plot([l[jz, 2] for jz in 1:h.Nz], Layout(; xaxis=attr(title="ω", zeroline=true), font_size=16, title="Cons from ext ϕ"))
	pcef = plot([l[jz, 3] for jz in 1:h.Nz], Layout(; xaxis=attr(title="ω", zeroline=true, range=[h.ωmin, ωmax_show]), font_size=16, title="Cons from ext ϕ, fixed pN"))
	pv = plot([l[jz, 4] for jz in 1:h.Nz], Layout(; xaxis=attr(title="ω", zeroline=true, range=[h.ωmin, ωmax_show]), font_size=16, title="Value function"))

	pC = plot(scatter(;x=h.zgrid, y=Cz, showlegend=false), Layout(;xaxis_title="Z", font_size=16, title="Agg Consumption"))
	pCf = plot(scatter(;x=h.zgrid, y=Cz_fix, showlegend=false), Layout(;xaxis_title="Z", font_size=16, title="Agg Consumption with fixed pN"))

	p = [pc pv; pce pcef; pC pCf]
	p.plot.layout["xlabel"] = "ω"
	p.plot.layout["width"] = 800
	p.plot.layout["height"] = 1200
	p.plot.layout["font_family"] = "Fira Sans Light"
	if remote
		path = pwd() * "/../../Graphs/"
		save(path * "p_hh_z.jld", "p", p)
	else
		path = pwd() * "/../Graphs/"
		savefig(p, path * "hh_z.pdf")
	end
	Void
end

function plot_hh_policies_b(h::Hank; remote::Bool=false)
	show_ϵ, show_μ, show_σ, show_w, show_ζ, show_z = mean(h.ϵgrid), mean(h.μgrid), mean(h.σgrid), mean(h.wgrid), h.ζgrid[1], mean(h.zgrid)

	knots = (h.ωgrid, h.ϵgrid, h.bgrid, h.μgrid, h.σgrid, h.wgrid, h.ζgrid, h.zgrid)
	itp_ϕc  = interpolate(knots, h.ϕc, Gridded(Linear()))
	itp_vf  = interpolate(knots, h.vf, Gridded(Linear()))
	knots = (h.ωgrid, h.ϵgrid, h.bgrid, h.μgrid, h.σgrid, h.wgrid, h.ζgrid, h.zgrid, h.pngrid)
	itp_ϕc_ext  = interpolate(knots, h.ϕc_ext, Gridded(Linear()))

	itp_pN = make_itp(h, h.pN, agg=true)

	l = Array{PlotlyBase.GenericTrace{Dict{Symbol,Any}}}(h.Nb, 4)
	Cb = Vector{Float64}(h.Nb)
	Cb_fix = Vector{Float64}(h.Nb)
	for (jb, bv) in enumerate(h.bgrid)
		ϕc_vec = zeros(h.Nω)
		ϕce_vec = zeros(h.Nω)
		ϕce_vecfix = zeros(h.Nω)
		vf_vec = zeros(h.Nω)
		show_pN = itp_pN[bv, show_μ, show_σ, show_w, 1., ceil(Int, h.Nz/2)]
		for (jω, ωv) in enumerate(h.ωgrid)
			ϕc_vec[jω] = itp_ϕc[ωv, show_ϵ, bv, show_μ, show_σ, show_w, show_ζ, show_z]
			ϕce_vec[jω] = itp_ϕc_ext[ωv, show_ϵ, bv, show_μ, show_σ, show_w, show_ζ, show_z, show_pN]
			ϕce_vecfix[jω] = itp_ϕc_ext[ωv, show_ϵ, bv, show_μ, show_σ, show_w, show_ζ, show_z, mean(h.pngrid)]
			vf_vec[jω] = itp_vf[ωv, show_ϵ, bv, show_μ, show_σ, show_w, show_ζ, show_z]
		end
		l_new = scatter(;x=h.ωgrid, y=ϕc_vec, line_shape="spline", name="b = $(round(bv,2))", marker_color=col[ceil(Int,10*jb/h.Nb)])
		l[jb,1] = l_new
		l_new = scatter(;x=h.ωgrid, y=ϕce_vec, line_shape="spline", name="b = $(round(bv,2))", showlegend=false, marker_color=col[ceil(Int,10*jb/h.Nb)])
		l[jb,2] = l_new
		l_new = scatter(;x=h.ωgrid, y=ϕce_vecfix, line_shape="spline", name="b = $(round(bv,2))", showlegend=false, marker_color=col[ceil(Int,10*jb/h.Nb)])
		l[jb,3] = l_new
		l_new = scatter(;x=h.ωgrid, y=vf_vec, line_shape="spline", name="b = $(round(bv,2))", showlegend=false, marker_color=col[ceil(Int,10*jb/h.Nb)])
		l[jb,4] = l_new

		ωmin_int, ωmax_int = quantile.(LogNormal(show_μ, show_σ), [.005; .995]) + h.ωmin
		val_int_C, val_int_Cfix = 0., 0.
		for (jϵ, ϵv) in enumerate(h.ϵgrid)
			f(ω) = pdf(LogNormal(show_μ, show_σ), ω-h.ωmin) * itp_ϕc_ext[ω, ϵv, bv, show_μ, show_σ, show_w, show_ζ, show_z, show_pN]
			ffix(ω) = pdf(LogNormal(show_μ, show_σ), ω-h.ωmin) * itp_ϕc_ext[ω, ϵv, bv, show_μ, show_σ, show_w, show_ζ, show_z, mean(h.pngrid)]
			(val, err) = hquadrature(f, ωmin_int, ωmax_int, reltol=1e-12, abstol=0, maxevals=0)
			(valfix, err) = hquadrature(ffix, ωmin_int, ωmax_int, reltol=1e-12, abstol=0, maxevals=0)
			val_int_C += val * h.λϵ[jϵ] 
			val_int_Cfix += valfix * h.λϵ[jϵ] 
		end

		Cb[jb], Cb_fix[jb] = val_int_C, val_int_Cfix
	end
	ωmax_show = min(h.ωmax, quantile(LogNormal(show_μ, show_σ), 0.999)+h.ωmin)

	pc = plot([l[jb, 1] for jb in 1:h.Nb], Layout(; xaxis=attr(title="ω", zeroline=true, range=[h.ωmin, ωmax_show]), font_size=16, title="Consumption"))
	pce = plot([l[jb, 2] for jb in 1:h.Nb], Layout(; xaxis=attr(title="ω", zeroline=true), font_size=16, title="Cons from ext ϕ"))
	pcef = plot([l[jb, 3] for jb in 1:h.Nb], Layout(; xaxis=attr(title="ω", zeroline=true, range=[h.ωmin, ωmax_show]), font_size=16, title="Cons from ext ϕ, fixed pN"))
	pv = plot([l[jb, 4] for jb in 1:h.Nb], Layout(; xaxis=attr(title="ω", zeroline=true, range=[h.ωmin, ωmax_show]), font_size=16, title="Value function"))

	pC = plot(scatter(;x=h.bgrid, y=Cb, showlegend=false), Layout(;xaxis_title="B", font_size=16, title="Agg Consumption"))
	pCf = plot(scatter(;x=h.bgrid, y=Cb_fix, showlegend=false), Layout(;xaxis_title="B", font_size=16, title="Agg Consumption with fixed pN"))

	p = [pc pv; pce pcef; pC pCf]
	p.plot.layout["xlabel"] = "ω"
	p.plot.layout["width"] = 800
	p.plot.layout["height"] = 1200
	p.plot.layout["font_family"] = "Fira Sans Light"

	if remote
		path = pwd() * "/../../Graphs/"
		save(path * "p_hh_b.jld", "p", p)
	else
		path = pwd() * "/../Graphs/"
		savefig(p, path * "hh_b.pdf")
	end
	Void
end

function lines(h::Hank, y, x_dim, name=""; custom_w::Int=0)
	jshow_b, jshow_μ, jshow_σ, jshow_w, jshow_ζ, jshow_z = ceil(Int, h.Nb/2), ceil(Int, h.Nμ/2), ceil(Int, h.Nσ/2), floor(Int, h.Nw/2), 1, ceil(Int, h.Nz/2)

	if custom_w != 0
		jshow_w = custom_w
	end

	x = h.bgrid
	xlabel = "B"
	if x_dim == 1
		y = y[:, jshow_μ, jshow_σ, jshow_w, jshow_ζ, jshow_z]
	elseif x_dim == 2
		x, xlabel = h.μgrid, "μ"
		y = y[jshow_b, :, jshow_σ, jshow_w, jshow_ζ, jshow_z]
	elseif x_dim == 3
		x, xlabel = h.σgrid, "σ"
		y = y[jshow_b, jshow_μ, :, jshow_w, jshow_ζ, jshow_z]
	elseif x_dim == 4
		x, xlabel = h.wgrid, "w"
		y = y[jshow_b, jshow_μ, jshow_σ, :, jshow_ζ, jshow_z]
	elseif x_dim == 6
		x, xlabel = h.zgrid, "z"
		y = y[jshow_b, jshow_μ, jshow_σ, jshow_w, jshow_ζ, :]
	else
		print_save("x_dim wrong")
	end


	layout = Layout(;	xaxis=attr(title=xlabel, zeroline=true),
						yaxis=attr(zeroline=true),
						font_size=16, font_family="Fira Sans Light")

	l = scatter(;x=x, y=y, showlegend=false)
	p = plot(l, layout)
	if name == ""
	else
		p.plot.layout["title"]=name
	end
	return p
end

function plot_gov_welf(h::Hank; remote::Bool=false)
	itp_vf = make_itp(h, h.vf; agg=false)

	B′_mat = reshape(h.issuance, h.Nb, h.Nμ, h.Nσ, h.Nw, h.Nζ, h.Nz)
	μ′_mat = reshape(h.μ′, h.Nb, h.Nμ, h.Nσ, h.Nw, h.Nζ, h.Nz, h.Nz, 2)
	σ′_mat = reshape(h.σ′, h.Nb, h.Nμ, h.Nσ, h.Nw, h.Nζ, h.Nz, h.Nz, 2)
	w′_mat = reshape(h.wage, h.Nb, h.Nμ, h.Nσ, h.Nw, h.Nζ, h.Nz)

	Wr_vec = zeros(size(h.Jgrid, 1))
	Wd_vec = zeros(size(h.Jgrid, 1))
	for js in 1:length(Wr_vec)
		jb = h.Jgrid[js, 1]
		jμ = h.Jgrid[js, 2]
		jσ = h.Jgrid[js, 3]
		jw = h.Jgrid[js, 4]
		jζ = 1
		jz = h.Jgrid[js, 6]

		EWr, EWd = 0., 0.

		bvp = B′_mat[jb, jμ, jσ, jw, jζ, jz]
		wvp = w′_mat[jb, jμ, jσ, jw, jζ, jz]
		for jzp in 1:h.Nz
			μvp = μ′_mat[jb, jμ, jσ, jw, jζ, jz, jzp, 1]
			σvp = σ′_mat[jb, jμ, jσ, jw, jζ, jz, jzp, 1]
			EWr += h.Pz[jz, jzp] * integrate_itp(h, bvp, μvp, σvp, wvp, 1, jzp, itp_vf)
			μvp = μ′_mat[jb, jμ, jσ, jw, jζ, jz, jzp, 2]
			σvp = σ′_mat[jb, jμ, jσ, jw, jζ, jz, jzp, 2]
			EWd += h.Pz[jz, jzp] * integrate_itp(h, (1.-h.ℏ)*bvp, μvp, σvp, wvp, 2, jzp, itp_vf)
		end

		Wr_vec[js] = EWr
		Wd_vec[js] = EWd
	end

	Wr_mat = reshape(Wr_vec, h.Nb, h.Nμ, h.Nσ, h.Nw, h.Nζ, h.Nz)
	Wd_mat = reshape(Wd_vec, h.Nb, h.Nμ, h.Nσ, h.Nw, h.Nζ, h.Nz)
	tempt  = Wd_mat - Wr_mat

	pWr1 = lines(h, Wr_mat, 1, "Expected welfare in repayment")
	pWr2 = lines(h, Wr_mat, 2)
	pWr3 = lines(h, Wr_mat, 3)
	pWr4 = lines(h, Wr_mat, 4)
	pWr6 = lines(h, Wr_mat, 6)
	pWd1 = lines(h, Wd_mat, 1, "Expected welfare in default")
	pWd2 = lines(h, Wd_mat, 2)
	pWd3 = lines(h, Wd_mat, 3)
	pWd4 = lines(h, Wd_mat, 4)
	pWd6 = lines(h, Wd_mat, 6)
	ptp1 = lines(h, tempt,  1, "Expected default incentives")
	ptp2 = lines(h, tempt,  2)
	ptp3 = lines(h, tempt,  3)
	ptp4 = lines(h, tempt,  4)
	ptp6 = lines(h, tempt,  6)

	p = [pWr1 pWr2 pWr3 pWr4 pWr6; pWd1 pWd2 pWd3 pWd4 pWd6; ptp1 ptp2 ptp3 ptp4 ptp6]
	p.plot.layout["width"] = 800
	p.plot.layout["height"] = 800/1.5
	p.plot.layout["font_family"] = "Fira Sans Light"
	if remote
		path = pwd() * "/../../Graphs/"
		save(path * "p_objfunc.jld", "p", p)
	else
		savefig(p, pwd() * "/../Graphs/objfunc.pdf")
	end
	Void
end

function plot_govt_reaction(h::Hank; remote::Bool=false)
	jμ, jσ, jw = ceil(Int, h.Nμ/2), ceil(Int, h.Nσ/2), ceil(Int, h.Nw/2)
	μv, σv, wv = h.μgrid[jμ], h.σgrid[jσ], h.wgrid[jw]
	jζ = 1

	itp_vf = make_itp(h, h.vf; agg=false)

	B′_mat = reshape(h.issuance, h.Nb, h.Nμ, h.Nσ, h.Nw, h.Nζ, h.Nz)
	μ′_mat = reshape(h.μ′, h.Nb, h.Nμ, h.Nσ, h.Nw, h.Nζ, h.Nz, h.Nz, 2)
	σ′_mat = reshape(h.σ′, h.Nb, h.Nμ, h.Nσ, h.Nw, h.Nζ, h.Nz, h.Nz, 2)
	w′_mat = reshape(h.wage, h.Nb, h.Nμ, h.Nσ, h.Nw, h.Nζ, h.Nz)

	states = gridmake([1; h.Nb], [1; h.Nz])
	# p_vec = Array{PlotlyJS.SyncPlot{PlotlyJS.ElectronDisplay}}(size(states,1))
	p_vec = Array{PlotlyJS.SyncPlot}(size(states,1))
	for js in 1:size(states,1)
		Wr = zeros(h.Nz)
		Wd = zeros(h.Nz)
		jb, jz = states[js, :]
		bvp = B′_mat[jb, jμ, jσ, jw, jζ, jz]
		wvp = w′_mat[jb, jμ, jσ, jw, jζ, jz]
		for jzp in 1:h.Nz
			μvp = μ′_mat[jb, jμ, jσ, jw, jζ, jz, jzp, 1]
			σvp = σ′_mat[jb, jμ, jσ, jw, jζ, jz, jzp, 1]
			Wr[jzp] = integrate_itp(h, bvp, μvp, σvp, wvp, 1, jzp, itp_vf)
			μvp = μ′_mat[jb, jμ, jσ, jw, jζ, jz, jzp, 2]
			σvp = σ′_mat[jb, jμ, jσ, jw, jζ, jz, jzp, 2]
			Wd[jzp] = integrate_itp(h, (1.-h.ℏ)*bvp, μvp, σvp, wvp, 2, jzp, itp_vf)
		end
		p_vec[js] = plot(  [scatter(;x=h.zgrid, y=Wr, marker_color=col[1], showlegend=false),
						scatter(;x=h.zgrid, y=Wd, marker_color=col[4], showlegend=false, line_dash="dashdot")],
						Layout(;title="B=$(h.bgrid[jb]), z=$(exp(h.zgrid[jz]))"))
	end

	p = [p_vec[1] p_vec[2]; p_vec[3] p_vec[4]]
	p.plot.layout["width"] = 800
	p.plot.layout["height"] = 800
	p.plot.layout["font_family"] = "Fira Sans Light"
	if remote
		path = pwd() * "/../../Graphs/"
		save(path * "p_reactions.jld", "p", p)
	else
		savefig(p, pwd() * "/../Graphs/reactions.pdf")
	end
	Void
end

function plot_debtprice(h::Hank; remote::Bool=false)

	qʰ_mat, qᵍ_mat, wL_mat, T_mat, pC_mat, Π_mat = _unpackstatefs(h)
	T_vec  = reshape(T_mat, length(T_mat))
	Π_vec  = reshape(Π_mat, length(Π_mat))
	wL_vec = reshape(wL_mat, length(wL_mat))

	ϕc_mat = h.ϕc
	yd_mat = zeros(h.ϕc)
	pC_big = zeros(h.ϕc)

	adj = sum(h.λϵ.*exp.(h.ϵgrid))
	agg_income = wL_vec + Π_vec / adj

	def_prob = zeros(h.Nb, h.Nμ, h.Nσ, h.Nw, h.Nζ, h.Nz)
	rep_mat = reshape(h.repay, h.Nb, h.Nμ, h.Nσ, h.Nw, h.Nζ, h.Nz, h.Nz)

	jpn = ceil(Int, length(h.pngrid)/2)
	pnv = h.pngrid[jpn]
	N = size(h.Jgrid, 1)
	wage_pn, labor_pn, profits_pn = Array{Float64, 1}(N), Array{Float64, 1}(N), Array{Float64, 1}(N)
	for js in 1:N
		jw = h.Jgrid[js, 4]
		jζ = h.Jgrid[js, 5]
		jz = h.Jgrid[js, 6]

		wv = h.wgrid[jw]
		ζv = h.ζgrid[jζ]
		zv = h.zgrid[jz]

		labor_pn[js], wage_pn[js], profits_pn[js], _ = labor_market(h, ζv, zv, wv, pnv)
	end

	pC = price_index(h, pnv)
	pC_fix = ones(h.Nb, h.Nμ, h.Nσ, h.Nw, h.Nζ, h.Nz) * pC

	T_fix = govt_bc(h, wage_pn.*labor_pn)# - reshape(profits_pn - h.profits, h.Nb, h.Nμ, h.Nσ, h.Nw, h.Nζ, h.Nz)
	Π_fix = reshape(profits_pn, h.Nb, h.Nμ, h.Nσ, h.Nw, h.Nζ, h.Nz)

	wL_fix  = reshape(wage_pn.*labor_pn, h.Nb, h.Nμ, h.Nσ, h.Nw, h.Nζ, h.Nz) * (1.0 - h.τ)

	yd_fix = zeros(h.ϕc)
    pC_bigfix = zeros(h.ϕc)
	for js in 1:size(h.Jgrid, 1)
		jb = h.Jgrid[js, 1]
		jμ = h.Jgrid[js, 2]
		jσ = h.Jgrid[js, 3]
		jw = h.Jgrid[js, 4]
		jζ = h.Jgrid[js, 5]
		jz = h.Jgrid[js, 6]
		for (jϵ, ϵv) in enumerate(h.ϵgrid), (jω, ωv) in enumerate(h.ωgrid)
			yd_mat[jω, jϵ, jb, jμ, jσ, jw, jζ, jz] = ωv + agg_income[js] * exp(ϵv) - T_vec[js]
			pC_big[jω, jϵ, jb, jμ, jσ, jw, jζ, jz] = pC_mat[js]

			yd_fix[jω, jϵ, jb, jμ, jσ, jw, jζ, jz] = ωv + (wL_fix[jb, jμ, jσ, jw, jζ, jz] + Π_fix[jb, jμ, jσ, jw, jζ, jz]/adj) * exp(ϵv) - T_fix[jb, jμ, jσ, jw, jζ, jz]
			pC_bigfix[jω, jϵ, jb, jμ, jσ, jw, jζ, jz] = pC_fix[jb, jμ, jσ, jw, jζ, jz]
		end
		for jzp in 1:h.Nz
			def_prob[jb, jμ, jσ, jw, jζ, jz] += h.Pz[jz, jzp] * (1.-rep_mat[jb, jμ, jσ, jw, jζ, jz, jzp])
		end
	end
	ϕc_ext_mat = h.ϕc_ext[:,:,:,:,:,:,:,:,jpn]
	
	Srate = 1. - pC_big .* ϕc_mat ./ yd_mat
	Sratef= 1. - pC_bigfix .* ϕc_ext_mat ./ yd_mat

	pq1 = lines(h, qᵍ_mat,  1, "Price of government debt")
	pq2 = lines(h, qᵍ_mat,  2)
	pq3 = lines(h, qᵍ_mat,  3)
	pq4 = lines(h, qᵍ_mat,  4)
	pq6 = lines(h, qᵍ_mat,  6)
	
	pd1 = lines(h, def_prob,  1, "One-period def prob")
	pd2 = lines(h, def_prob,  2)
	pd3 = lines(h, def_prob,  3)
	pd4 = lines(h, def_prob,  4)
	pd6 = lines(h, def_prob,  6)

	jω1, jω2 = 1, ceil(Int, h.Nω / 2)
	jϵ_show = ceil(Int, h.Nϵ/2)
	pc1p = lines(h, Sratef[jω1, jϵ_show,:,:,:,:,:,:],  1, "Saving rate at ω = $(round(h.ωgrid[jω1],2))")
	pc2p = lines(h, Sratef[jω1, jϵ_show,:,:,:,:,:,:],  2)
	pc3p = lines(h, Sratef[jω1, jϵ_show,:,:,:,:,:,:],  3)
	pc4p = lines(h, Sratef[jω1, jϵ_show,:,:,:,:,:,:],  4)
	pc6p = lines(h, Sratef[jω1, jϵ_show,:,:,:,:,:,:],  6)

	pc1r = lines(h, Sratef[jω2, jϵ_show,:,:,:,:,:,:],  1, "Saving rate at ω = $(round(h.ωgrid[jω2],2))")
	pc2r = lines(h, Sratef[jω2, jϵ_show,:,:,:,:,:,:],  2)
	pc3r = lines(h, Sratef[jω2, jϵ_show,:,:,:,:,:,:],  3)
	pc4r = lines(h, Sratef[jω2, jϵ_show,:,:,:,:,:,:],  4)
	pc6r = lines(h, Sratef[jω2, jϵ_show,:,:,:,:,:,:],  6)
	
#	pc1pf = lines(h, Sratef[jω1, jϵ_show,:,:,:,:,:,:],  1, "S/Y at ω = $(round(h.ωgrid[jω1],2)), fixed pN")
#	pc2pf = lines(h, Sratef[jω1, jϵ_show,:,:,:,:,:,:],  2)
#	pc3pf = lines(h, Sratef[jω1, jϵ_show,:,:,:,:,:,:],  3)
#	pc4pf = lines(h, Sratef[jω1, jϵ_show,:,:,:,:,:,:],  4)
#	pc6pf = lines(h, Sratef[jω1, jϵ_show,:,:,:,:,:,:],  6)

#	pc1rf = lines(h, Sratef[jω2, jϵ_show,:,:,:,:,:,:],  1, "S/Y at ω = $(round(h.ωgrid[jω2],2)), fixed pN")
#	pc2rf = lines(h, Sratef[jω2, jϵ_show,:,:,:,:,:,:],  2)
#	pc3rf = lines(h, Sratef[jω2, jϵ_show,:,:,:,:,:,:],  3)
#	pc4rf = lines(h, Sratef[jω2, jϵ_show,:,:,:,:,:,:],  4)
#	pc6rf = lines(h, Sratef[jω2, jϵ_show,:,:,:,:,:,:],  6)


	p = [pq1 pq2 pq3 pq4 pq6; pd1 pd2 pd3 pd4 pd6; pc1p pc2p pc3p pc4p pc6p; pc1r pc2r pc3r pc4r pc6r]#; pc1pf pc2pf pc3pf pc4pf pc6pf; pc1rf pc2rf pc3rf pc4rf pc6rf]
	p.plot.layout["width"] = 800
	p.plot.layout["height"] = 800/1.15
	p.plot.layout["font_family"] = "Fira Sans Light"
	if remote
		path = pwd() * "/../../Graphs/"
		save(path * "p_debtprice.jld", "p", p)
	else
		savefig(p, pwd() * "/../Graphs/debtprice.pdf")
	end
	Void
end

function plot_eulereq(h::Hank; remote::Bool=false)
	ExpRealRet = zeros(h.Ns, h.Nz, 2)
	EZ = zeros(h.Nω, h.Nϵ, 1, h.Nϵ, h.Nz, 2)
	EIS = zeros(h.Nω, h.Nϵ, 1, h.Nϵ, h.Nz, 2)

	pC_vec = price_index(h, h.pN)

	itp_qᵍ = make_itp(h, h.qᵍ, agg=true)
	itp_pC = make_itp(h, pC_vec, agg=true)
	itp_ϕc = make_itp(h, h.ϕc, agg=false)
	itp_vf = make_itp(h, h.vf, agg=false)

	for (js, js_show) in enumerate(572)
		jb = h.Jgrid[js_show, 1]
		jμ = h.Jgrid[js_show, 2]
		jσ = h.Jgrid[js_show, 3]
		jw = h.Jgrid[js_show, 4]
		jζ = h.Jgrid[js_show, 5]
		jz = h.Jgrid[js_show, 6]

		bp = h.issuance[js_show]
		μp = h.μ′[js_show,:,:]
		σp = h.σ′[js_show,:,:]
		wpv = h.wage[js_show]

		pCv = price_index(h, h.pN[js_show])
		for jzp in 1:h.Nz
			# In repayment
			bpv = bp
			μpv = μp[jzp, 1]
			σpv = σp[jzp, 1]
			Rb = h.κ + (1.-h.ρ) * itp_qᵍ[bpv, μpv, σpv, wpv, 1, jzp]
			ExpRealRet[js, jzp, 1] = Rb * itp_pC[bpv, μpv, σpv, wpv, 1, jzp] / pCv

			# In default
			haircut = (1.-h.ℏ*(jζ==1))
			bpv = haircut * bp
			μpv = μp[jzp, 2]
			σpv = σp[jzp, 2]
			Rb = h.κ + (1.-h.ρ) * haircut * itp_qᵍ[bpv, μpv, σpv, wpv, 2, jzp]
			ExpRealRet[js, jzp, 1] = Rb * itp_pC[bpv, μpv, σpv, wpv, 2, jzp] / pCv
		end

		rep_mat = reshape(h.repay, h.Nb, h.Nμ, h.Nσ, h.Nw, h.Nζ, h.Nz, h.Nz)
		for (jω, ωv) in enumerate(h.ωgrid)
			for jϵ in 1:h.Nϵ
				Tvf = 0.
				V = zeros(h.Nϵ, h.Nz, 2)
				Cv = h.ϕc[jω, jϵ, jb, jμ, jσ, jw, jζ, jz]
				Vf = h.vf[jω, jϵ, jb, jμ, jσ, jw, jζ, jz]

				A = h.ϕa[jω, jϵ, jb, jμ, jσ, jw, jζ, jz]
				B = h.ϕb[jω, jϵ, jb, jμ, jσ, jw, jζ, jz]
				for jzp in 1:h.Nz
					rep_prob = rep_mat[jb, jμ, jσ, jw, jζ, jz, jzp] * (jζ == 1) + h.θ * (jζ == 2)

					# First in repayment
					bpv = bp
					μpv = μp[jzp, 1]
					σpv = σp[jzp, 1]
					R = h.κ + (1.-h.ρ) * itp_qᵍ[bpv, μpv, σpv, wpv, 1, jzp]
					
					ωpv = A + R * B
					for jϵp in 1:h.Nϵ
						V_t = itp_vf[ωpv, jϵp, bpv, μpv, σpv, wpv, 1, jzp]
						V[jϵp, jzp, 1] = V_t
						EIS[jω, jϵ, js, jϵp, jzp, 1] = (itp_ϕc[ωpv, jϵp, bpv, μpv, σpv, wpv, 1, jzp] / Cv)^(-1./h.ψ)
						Tvf += V_t^(1.-h.γ) * h.Pϵ[jϵ, jϵp] * h.Pz[jz, jzp] * rep_prob
					end
					
					# Then in default
					haircut = (1.-h.ℏ*(jζ==1))
					bpv = haircut * bp
					μpv = μp[jzp, 2]
					σpv = σp[jzp, 2]
					R = h.κ + (1.-h.ρ) * haircut * itp_qᵍ[bpv, μpv, σpv, wpv, 2, jzp]
					
					ωpv = A + R * B
					for jϵp in 1:h.Nϵ
						V_t = itp_vf[ωpv, jϵp, bpv, μpv, σpv, wpv, 2, jzp]
						V[jϵp, jzp, 2] = V_t
						EIS[jω, jϵ, js, jϵp, jzp, 2] = (itp_ϕc[ωpv, jϵp, bpv, μpv, σpv, wpv, 2, jzp] / Cv)^(-1./h.ψ)
						Tvf += V_t^(1.-h.γ) * h.Pϵ[jϵ, jϵp] * h.Pz[jz, jzp] * (1.-rep_prob)
						
						EZ[jω, jϵ, js, jϵp, jzp, 1] = (V[jϵp, jzp, 1] ./ Tvf).^(1./h.ψ - h.γ)
						EZ[jω, jϵ, js, jϵp, jzp, 2] = (V[jϵp, jzp, 2] ./ Tvf).^(1./h.ψ - h.γ)
					end
				end
			end
		end
	end
	SDF = EZ .* EIS

	Void
end

function plot_aggcons(h::Hank; remote::Bool=false)
	jμ, jσ, jw = ceil(Int, h.Nμ/2), ceil(Int, h.Nσ/2), ceil(Int, h.Nw/2)
	μv, σv, wv = h.μgrid[jμ], h.σgrid[jσ], h.wgrid[jw]
	jζ = 1

	itp_ϕc = make_itp(h, h.ϕc; agg=false)
	itp_ϕc2 = make_itp(h, h.ϕc.^2; agg=false)

	B′_mat = reshape(h.issuance, h.Nb, h.Nμ, h.Nσ, h.Nw, h.Nζ, h.Nz)
	μ′_mat = reshape(h.μ′, h.Nb, h.Nμ, h.Nσ, h.Nw, h.Nζ, h.Nz, h.Nz, 2)
	σ′_mat = reshape(h.σ′, h.Nb, h.Nμ, h.Nσ, h.Nw, h.Nζ, h.Nz, h.Nz, 2)
	w′_mat = reshape(h.wage, h.Nb, h.Nμ, h.Nσ, h.Nw, h.Nζ, h.Nz)

	states = gridmake([1; h.Nb], [1; h.Nz])
	# p_vec = Array{PlotlyJS.SyncPlot{PlotlyJS.ElectronDisplay}}(size(states,1))
	p_vec = Array{PlotlyJS.SyncPlot}(size(states,1))
	# p2_vec = Array{PlotlyJS.SyncPlot{PlotlyJS.ElectronDisplay}}(size(states,1))
	p2_vec = Array{PlotlyJS.SyncPlot}(size(states,1))
	for js in 1:size(states,1)
		C_r = zeros(h.Nz)
		VarCr = zeros(h.Nz)
		C_d = zeros(h.Nz)
		VarCd = zeros(h.Nz)
		jb, jz = states[js, :]
		bvp = B′_mat[jb, jμ, jσ, jw, jζ, jz]
		wvp = w′_mat[jb, jμ, jσ, jw, jζ, jz]
		for jzp in 1:h.Nz
			μvp = μ′_mat[jb, jμ, jσ, jw, jζ, jz, jzp, 1]
			σvp = σ′_mat[jb, jμ, jσ, jw, jζ, jz, jzp, 1]
			C_r[jzp] = integrate_itp(h, bvp, μvp, σvp, wvp, 1, jzp, itp_ϕc)
			VarCr[jzp] = integrate_itp(h, bvp, μvp, σvp, wvp, 1, jzp, itp_ϕc2) - C_r[jzp]^2
			μvp = μ′_mat[jb, jμ, jσ, jw, jζ, jz, jzp, 2]
			σvp = σ′_mat[jb, jμ, jσ, jw, jζ, jz, jzp, 2]
			C_d[jzp] = integrate_itp(h, (1.-h.ℏ)*bvp, μvp, σvp, wvp, 2, jzp, itp_ϕc)
			VarCd[jzp] = integrate_itp(h, (1.-h.ℏ)*bvp, μvp, σvp, wvp, 2, jzp, itp_ϕc2) - C_d[jzp]^2
		end
		p_vec[js] = plot(  [scatter(;x=h.zgrid, y=C_r, marker_color=col[1], showlegend=false),
						scatter(;x=h.zgrid, y=C_d, marker_color=col[4], showlegend=false, line_dash="dashdot")],
						Layout(;title="B=$(h.bgrid[jb]), z=$(exp(h.zgrid[jz]))"))
		p2_vec[js] = plot(  [scatter(;x=h.zgrid, y=VarCr, marker_color=col[1], showlegend=false),
						scatter(;x=h.zgrid, y=VarCd, marker_color=col[4], showlegend=false, line_dash="dashdot")],
						Layout(;title="B=$(h.bgrid[jb]), z=$(exp(h.zgrid[jz]))"))
	end

	p = [p_vec[1] p_vec[2]; p_vec[3] p_vec[4]]
	p.plot.layout["width"] = 800
	p.plot.layout["height"] = 800
	p.plot.layout["font_family"] = "Fira Sans Light"
	if remote
		path = pwd() * "/../../Graphs/"
		save(path * "p_aggcons.jld", "p", p)
	else
		savefig(p, pwd() * "/../Graphs/aggcons.pdf")
	end
	p2 = [p2_vec[1] p2_vec[2]; p2_vec[3] p2_vec[4]]
	p2.plot.layout["width"] = 800
	p2.plot.layout["height"] = 800
	p2.plot.layout["font_family"] = "Fira Sans Light"
	if remote
		path = pwd() * "/../../Graphs/"
		save(path * "p_varcons.jld", "p2", p2)
	else
		savefig(p2, pwd() * "/../Graphs/varcons.pdf")
	end
	Void
end

function plot_state_funcs(h::Hank; remote::Bool=false)

	pN_mat = reshape(h.pN,     h.Nb, h.Nμ, h.Nσ, h.Nw, h.Nζ, h.Nz)
	w_mat  = reshape(h.wage,   h.Nb, h.Nμ, h.Nσ, h.Nw, h.Nζ, h.Nz)
	Ld_mat = reshape(h.Ld,     h.Nb, h.Nμ, h.Nσ, h.Nw, h.Nζ, h.Nz)
	Y_mat  = reshape(h.output, h.Nb, h.Nμ, h.Nσ, h.Nw, h.Nζ, h.Nz)
	Π_mat  = reshape(h.profits,h.Nb, h.Nμ, h.Nσ, h.Nw, h.Nζ, h.Nz)

	T_mat  = govt_bc(h, h.wage.*h.Ld)


	for (jp, jw) in enumerate([1; h.Nw])
		ppN1 = lines(h, pN_mat, 1, "Price of nontradables"; custom_w = jw)
		pw1  = lines(h, w_mat, 1, "Wage"; custom_w = jw)
		pLd1 = lines(h, Ld_mat, 1, "Labor supply"; custom_w = jw)
		pY1  = lines(h, Y_mat, 1, "Output"; custom_w = jw)
		pΠ1  = lines(h, Π_mat, 1, "Profits"; custom_w = jw)
		pT1  = lines(h, T_mat, 1, "Taxes"; custom_w = jw)

		ppN2 = lines(h, pN_mat, 2; custom_w = jw)
		pw2  = lines(h, w_mat, 2; custom_w = jw)
		pLd2 = lines(h, Ld_mat, 2; custom_w = jw)
		pY2  = lines(h, Y_mat, 2; custom_w = jw)
		pΠ2  = lines(h, Π_mat, 2; custom_w = jw)
		pT2  = lines(h, T_mat, 2; custom_w = jw)

		ppN3 = lines(h, pN_mat, 3; custom_w = jw)
		pw3  = lines(h, w_mat, 3; custom_w = jw)
		pLd3 = lines(h, Ld_mat, 3; custom_w = jw)
		pY3  = lines(h, Y_mat, 3; custom_w = jw)
		pΠ3  = lines(h, Π_mat, 3; custom_w = jw)
		pT3  = lines(h, T_mat, 3; custom_w = jw)

		ppN4 = lines(h, pN_mat, 4; custom_w = jw)
		pw4  = lines(h, w_mat, 4; custom_w = jw)
		pLd4 = lines(h, Ld_mat, 4; custom_w = jw)
		pY4  = lines(h, Y_mat, 4; custom_w = jw)
		pΠ4  = lines(h, Π_mat, 4; custom_w = jw)
		pT4  = lines(h, T_mat, 4; custom_w = jw)

		ppN6 = lines(h, pN_mat, 6; custom_w = jw)
		pw6  = lines(h, w_mat, 6; custom_w = jw)
		pLd6 = lines(h, Ld_mat, 6; custom_w = jw)
		pY6  = lines(h, Y_mat, 6; custom_w = jw)
		pΠ6  = lines(h, Π_mat, 6; custom_w = jw)
		pT6  = lines(h, T_mat, 6; custom_w = jw)

		# p = [ppN1 pw1 pLd1; ppN2 pw2 pLd2; ppN3 pw3 pLd3; ppN4 pw4 pLd4; ppN6 pw6 pLd6]
		p = [ppN1 ppN2 ppN3 ppN4 ppN6; pw1 pw2 pw3 pw4 pw6; pLd1 pLd2 pLd3 pLd4 pLd6; pY1 pY2 pY3 pY4 pY6; pΠ1 pΠ2 pΠ3 pΠ4 pΠ6; pT1 pT2 pT3 pT4 pT6]
		p.plot.layout["width"] = 800
		p.plot.layout["height"] = 640/4*6
		p.plot.layout["font_family"] = "Fira Sans Light"
		if remote
			path = pwd() * "/../../Graphs/"
			save(path * "p_statefuncs$(jp).jld", "p", p)
		else
			savefig(p, pwd() * "/../Graphs/statefuncs$(jp).pdf")
		end
	end
	Void
end

function plot_LoM(h::Hank; remote::Bool=false)
	jz = ceil(Int, h.Nz/2)

	μ′_mat = zeros(h.Nb*h.Nμ*h.Nσ*h.Nw*h.Nζ*h.Nz)
	σ′_mat = zeros(h.Nb*h.Nμ*h.Nσ*h.Nw*h.Nζ*h.Nz)

	for js in 1:size(h.Jgrid, 1)
		jz = h.Jgrid[js, 6]

		μ′_mat[js] = h.μ′[js,jz,1]
		σ′_mat[js] = h.σ′[js,jz,1]
	end

	μ′_mat = reshape(μ′_mat, h.Nb, h.Nμ, h.Nσ, h.Nw, h.Nζ, h.Nz)
	σ′_mat = reshape(σ′_mat, h.Nb, h.Nμ, h.Nσ, h.Nw, h.Nζ, h.Nz)

	pμ1 = lines(h, μ′_mat, 1, "Next period μ")
	pσ1 = lines(h, σ′_mat, 1, "Next period σ")

	pμ2 = lines(h, μ′_mat, 2)
	pσ2 = lines(h, σ′_mat, 2)

	pμ3 = lines(h, μ′_mat, 3)
	pσ3 = lines(h, σ′_mat, 3)

	pμ4 = lines(h, μ′_mat, 4)
	pσ4 = lines(h, σ′_mat, 4)

	pμ6 = lines(h, μ′_mat, 6)
	pσ6 = lines(h, σ′_mat, 6)

	p = [pμ1 pμ2 pμ3 pμ4 pμ6; pσ1 pσ2 pσ3 pσ4 pσ6]
	p.plot.layout["width"] = 800
	p.plot.layout["height"] = 500
	p.plot.layout["font_family"] = "Fira Sans Light"

	if remote
		path = pwd() * "/../../Graphs/"
		save(path * "p_LoMs.jld", "p", p)
	else
		savefig(p, pwd() * "/../Graphs/LoMs.pdf")
	end
	Void
end

function labor_demand(h::Hank, w, tfp, pN; get_both::Bool = false)

	Ld_nontradables = (h.α_N * pN  ./ w).^(1.0/(1.0-h.α_N))
	Ld_tradables    = (h.α_T * tfp ./ w).^(1.0/(1.0-h.α_T))

	if get_both
		return Ld_nontradables, Ld_tradables
	else
		return Ld_nontradables + Ld_tradables
	end
end

function plot_labor_demand(h::Hank; remote::Bool=false)
	jz = ceil(Int, h.Nz/2)
	z_show = h.zgrid[jz]

	vl = 1e8
	l = scatter(;y=h.wgrid, x=ones(h.wgrid), line_dash="dashdot", marker_color="black", showlegend=false, mode="lines", title="Labor market")
	for (jpN, pNv) in enumerate(h.pngrid)
		Ld = labor_demand(h, h.wgrid, exp(z_show), pNv)
		label = "pₙ = $(round(pNv,2))"
		l = hcat(l, scatter(;y=h.wgrid, x=Ld, name=label, marker_color=col[jpN], line_shape="spline"))
		if minimum(Ld) < vl
			vl = minimum(Ld)
		end
	end
	shapes = [hline(minimum(h.wgrid), line_width=1)]
	layout = Layout(;	xaxis=attr(title="L", zeroline=true, range=[0., 3.]),
						yaxis=attr(title="w", zeroline=true),
						title="Labor Market",
						annotations=[attr(x=1, y=maximum(h.wgrid),text="Lˢ", xanchor="center", yanchor="bottom", showarrow=false, font_size=18)],
						shapes=shapes,
						font_size=16, font_family="Fira Sans Light")

	p = plot([l[jj] for jj in 1:length(l)], layout)
	p.plot.layout["width"] = 800
	p.plot.layout["height"] = 500

	if remote
		path = pwd() * "/../../Graphs/"
		save(path * "p_labordemand.jld", "p", p)
	else
		path = pwd() * "/../Graphs/"
		savefig(p, path * "labordemand.pdf")
	end
	Void
end

function plot_nontradables(h::Hank; remote::Bool=false)
	jb = ceil(Int, h.Nb/2)
	jμ = ceil(Int, h.Nμ/2)
	jσ = ceil(Int, h.Nσ/2)
	jw = ceil(Int, h.Nw/2)
	jζ = ceil(Int, h.Nζ/2)
	jz = ceil(Int, h.Nz/2)

	bv, μv, σv, wv, ζv, zv = h.bgrid[jb], h.μgrid[jμ], h.σgrid[jσ], h.wgrid[jw], h.ζgrid[jζ], h.zgrid[jz]

	G_mat = reshape(h.spending, h.Nb, h.Nμ, h.Nσ, h.Nw, h.Nζ, h.Nz)
	B_mat = reshape(h.issuance, h.Nb, h.Nμ, h.Nσ, h.Nw, h.Nζ, h.Nz)

	itp_ϕc = make_itp(h, h.ϕc_ext; agg = false)

	pNmin, pNmax = minimum(h.pngrid), maximum(h.pngrid)

	l = Array{PlotlyBase.GenericTrace{Dict{Symbol,Any}}}(2*h.Nb)
	maxq = 0.
	minq = 10.
	for (jb, bv) in enumerate(h.bgrid)
		sup = zeros(h.pngrid)
		dem = zeros(h.pngrid)
		G   = G_mat[jb, jμ, jσ, jw, jζ, jz]
		Bpv = B_mat[jb, jμ, jσ, jw, jζ, jz]
		for (jpn, pnv) in enumerate(h.pngrid)
			sup[jpn], dem[jpn] = mkt_clearing(h, itp_ϕc, G, Bpv, pnv, pNmin, pNmax, bv, μv, σv, wv, jζ, jz, (jζ!=1); get_both=true)
		end
		l[jb] = scatter(; y=h.pngrid, x=sup, marker_color=col[jb], name="B = $(round(bv, 2))")
		l[h.Nb+jb] = scatter(; y=h.pngrid, x=dem, marker_color=col[jb], name="B = $(round(bv, 2))", showlegend=false)
		maxq = max(max(maximum(dem), maximum(sup)), maxq)
		minq = min(min(minimum(dem), minimum(sup)), minq)
	end
	maxq = min(maxq * 1.10, 3.)
	minq = minq * 0.9

	p = plot([l[jb] for jb in 1:2*h.Nb], Layout(; yaxis_title="pₙ", xaxis_title="Q", xaxis_range=[0., maxq]))
	if remote
		path = pwd() * "/../../Graphs/"
		save(path * "p_nontradables_B.jld", "p", p)
	else
		path = pwd() * "/../Graphs/"
		savefig(p, path * "nontradables_B.pdf")
	end

	jb = ceil(Int, h.Nb/2)
	bv, μv, σv, wv, ζv, zv = h.bgrid[jb], h.μgrid[jμ], h.σgrid[jσ], h.wgrid[jw], h.ζgrid[jζ], h.zgrid[jz]
	l = Array{PlotlyBase.GenericTrace{Dict{Symbol,Any}}}(2*h.Nz,2)
	maxq = 0.
	minq = 10.
	for (jz, zv) in enumerate(h.zgrid)
		sup = zeros(h.pngrid)
		dem = zeros(h.pngrid)
		supN = zeros(h.pngrid)
		G   = G_mat[jb, jμ, jσ, jw, jζ, jz]
		Bpv = B_mat[jb, jμ, jσ, jw, jζ, jz]
		for (jpn, pnv) in enumerate(h.pngrid)
			sup[jpn], dem[jpn] = mkt_clearing(h, itp_ϕc, G, Bpv, pnv, pNmin, pNmax, bv, μv, σv, wv, jζ, jz, (jζ!=1); get_both=true)

			zv = h.zgrid[jz]
			Ld, w_new, profits, output = labor_market(h, jζ, zv, wv, pnv)
			Ld_N, _  = labor_demand(h, w_new, zv, jζ, pnv; get_both=true)
			supN[jpn] = TFP_N(zv, h.Δ, jζ) * Ld_N^(h.α_N)
		end
		l[jz,1] = scatter(; y=h.pngrid, x=sup, marker_color=col[ceil(Int,10*jz/h.Nz)], name="z = $(round(exp(zv), 2))")
		l[h.Nz+jz,1] = scatter(; y=h.pngrid, x=dem, marker_color=col[ceil(Int,10*jz/h.Nz)], name="z = $(round(exp(zv), 2))", showlegend=false)
		l[jz,2] = scatter(; x=supN, y=h.pngrid, marker_color=col[ceil(Int,10*jz/h.Nz)], name="z = $(round(exp(zv), 2))")
		maxq = max(max(maximum(dem), maximum(sup)), maxq)
		minq = min(min(minimum(dem), minimum(sup)), minq)
	end
	maxq = min(maxq * 1.10, 3.)
	minq = minq * 0.9

	p = plot([l[jz,1] for jz in 1:2*h.Nz], Layout(; yaxis_title="pₙ", xaxis_title="Q", xaxis_range=[0., maxq]))

	if remote
		path = pwd() * "/../../Graphs/"
		save(path * "p_nontradables_z.jld", "p", p)
		p = plot([l[jz,2] for jz in 1:h.Nz], Layout(;xaxis_title="Q", yaxis_title="pₙ", xaxis_range=[0., maxq]))
		save(path * "p_nontradables_z2.jld", "p", p)
	else
		path = pwd() * "/../Graphs/"
		savefig(p, path * "nontradables_z.pdf")
	end
	Void
end


function plot_convergence(dist_statefuncs, dist_LoMs, T::Int64; remote::Bool=false)

	l_funcs1 = scatter(; x=1:T, y = (dist_statefuncs[1:T, 1]), name = "w")
	l_funcs2 = scatter(; x=1:T, y = (dist_statefuncs[1:T, 2]), name = "pN")
	l_funcs3 = scatter(; x=1:T, y = (dist_statefuncs[1:T, 3]), name = "Ld")
	l_LoM1   = scatter(; x=1:T, y = (dist_LoMs[1:T,1]), name = "μ′")
	l_LoM2   = scatter(; x=1:T, y = (dist_LoMs[1:T,2]), name = "σ′")

	layout = Layout(;	xaxis=attr(title="𝑡", zeroline=true),
						yaxis_type="log",
						font_size=16, font_family="Fira Sans Light")


	p = plot([l_funcs1, l_funcs2, l_funcs3, l_LoM1, l_LoM2], layout)
	p.plot.layout["width"] = 800
	p.plot.layout["height"] = 500

	if remote
		path = pwd() * "/../../Graphs/"
		save(path * "p_conv.jld", "p", p)
	else
		path = pwd() * "/../Graphs/"
		savefig(p, path * "conv.pdf")
	end
	Void
end

function plot_outerdists(h; remote::Bool=false)
	T = length(h.outer_dists)

	p = plot(scatter(; x=1:T, y = h.outer_dists, showlegend=false),
		Layout(; xaxis_title="𝑡", yaxis_type="log", font_size=16, font_family="Fira Sans Light", width=800, height=500))

	if remote
		path = pwd() * "/../../Graphs/"
		save(path * "p_outconv.jld", "p", p)
	else
		path = pwd() * "/../Graphs/"
		savefig(p, path * "outconv.pdf")
	end
	Void
end

function plot_simul(path::Path; remote::Bool=false, trim::Int=0)
	name = ""
	if trim > 0
		trim_path!(path, trim)
	else
		name = "_full"
	end

	T = size(path.data, 1)

	B_vec = series(path,:B)
	μ_vec = series(path,:μ)
	σ_vec = series(path,:σ)
	w_vec = series(path,:w)
	ζ_vec = series(path,:ζ)-1
	z_vec = exp.(series(path,:z))
	Y_vec = series(path,:Y)
	L_vec = series(path,:L)
	π_vec = series(path,:π)
	P_vec = series(path,:P)
	Pe_vec= series(path,:Pe)
	ψ_vec = series(path,:ψ)
	A_vec = series(path,:A)
	Bf_vec= series(path,:Bf)
	Wr_vec= series(path,:Wr)
	Wd_vec= series(path,:Wd)

	shiftζ = [0; ζ_vec[1:end-1]]

	defaults = find((ζ_vec.==1) .* (shiftζ.==0))./4
	exits    = find((ζ_vec.==0) .* (shiftζ.==1))./4

	times = (1:T)./4

	default_shades = rect(defaults, exits, 0, 1; fillcolor="#d3d3d3", opacity=0.5, line_width=0, xref="x", yref="paper")

	pB = plot([	scatter(; x=times, y=B_vec, marker_color=col[1], showlegend=false),
				# scatter(; x=times, y=ones(times)*minimum(h.bgrid), showlegend=false, line_dash="dashdot", marker_color="black", line_width=0.5),
				scatter(; x=times, y=ones(times)*maximum(h.bgrid), showlegend=false, line_dash="dashdot", marker_color="black", line_width=0.5)],
						Layout(; title="Bonds", xaxis=attr(title="𝑡")));
	pμ = plot([ scatter(; x=times, y=μ_vec, marker_color=col[1], showlegend=false),
				# scatter(; x=times, y=ones(times)*minimum(h.μgrid), showlegend=false, line_dash="dashdot", marker_color="black", line_width=0.5),
				scatter(; x=times, y=ones(times)*maximum(h.μgrid), showlegend=false, line_dash="dashdot", marker_color="black", line_width=0.5)],
						Layout(; title="μ", xaxis=attr(title="𝑡")));
	pσ = plot([ scatter(; x=times, y=σ_vec, marker_color=col[1], showlegend=false),
				# scatter(; x=times, y=ones(times)*maximum(h.σgrid), showlegend=false, line_dash="dashdot", marker_color="black", line_width=0.5),
				scatter(; x=times, y=ones(times)*minimum(h.σgrid), showlegend=false, line_dash="dashdot", marker_color="black", line_width=0.5)],
						Layout(; title="σ", xaxis=attr(title="𝑡")));
	pw = plot([ scatter(; x=times, y=w_vec, marker_color=col[1], showlegend=false),
				scatter(; x=times, y=ones(times)*minimum(h.wgrid), showlegend=false, line_dash="dashdot", marker_color="black", line_width=0.5),
				scatter(; x=times, y=ones(times)*maximum(h.wgrid), showlegend=false, line_dash="dashdot", marker_color="black", line_width=0.5)],
						Layout(; title="Wage", xaxis=attr(title="𝑡")));
	pz = plot(scatter(; x=times, y=z_vec, marker_color=col[1], showlegend=false), Layout(; title="TFP", xaxis=attr(title="𝑡")));
	pY = plot([ scatter(; x=times, y=Y_vec, marker_color=col[1], showlegend=false),
				scatter(; x=times, y=L_vec, marker_color=col[2], showlegend=false, line_dash="dashdot")],
			Layout(; title="Output", xaxis=attr(title="𝑡")));
	pπ = plot([scatter(; x=times, y=ζ_vec, marker_color=col[1], showlegend=false),
				scatter(; x=times, y=π_vec, marker_color=col[2], showlegend=false, line_dash="dashdot")],
			Layout(; title="Default prob", xaxis=attr(title="𝑡")));
	pP = plot([ scatter(; x=times, y=P_vec, marker_color=col[1], showlegend=false),
				# scatter(; x=times, y=ones(times)*maximum(h.pngrid), showlegend=false, line_dash="dashdot", marker_color="black", line_width=0.5),
				# scatter(; x=times, y=ones(times)*minimum(h.pngrid), showlegend=false, line_dash="dashdot", marker_color="black", line_width=0.5),
				scatter(; x=times, y=Pe_vec,marker_color=col[4], showlegend=false, line_dash="dashdot")],
			Layout(; title="Price of nontradables", xaxis=attr(title="𝑡")));
	pψ = plot(scatter(; x=times, y=ψ_vec, marker_color=col[1],  showlegend=false), Layout(; title="Fraction domestic", xaxis=attr(title="𝑡")));
	pA = plot(scatter(; x=times, y=A_vec, marker_color=col[1],  showlegend=false), Layout(; title="Domestic risk-free debt", xaxis_title="𝑡"));
	pBf= plot(scatter(; x=times, y=Bf_vec, marker_color=col[1], showlegend=false), Layout(; title="Foreign debt", xaxis_title="𝑡"));
	pW = plot([ scatter(;x=times, y=Wr_vec, marker_color=col[1], showlegend=false),
				scatter(;x=times, y=Wd_vec, marker_color=col[2], showlegend=false, line_dash="dashdot")], Layout(;title="Welfare", xaxis_title="𝑡"));

	p = [pB pw pz; pY pμ pσ; pA pBf pψ; pπ pW pP]
	# p.plot.layout["shapes"] = default_shades
	p.plot.layout["width"] = 850
	p.plot.layout["height"] = 850
	p.plot.layout["font_family"] = "Fira Sans Light"

	name = "simul"*name
	if remote
		path = pwd() * "/../../Graphs/"
		save(path * "p_"*name*".jld", "p", p)
	else
		path = pwd() * "/../Graphs/"
		savefig(p, path*name*".pdf")
		savefig(p, path*name*".png")
	end

	Void
end
