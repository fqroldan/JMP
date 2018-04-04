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

	show_b, show_μ, show_σ, show_w, show_ζ, show_z = mean(h.bgrid), mean(h.μgrid), mean(h.σgrid), mean(h.wgrid), h.ζgrid[1], h.zgrid[end]

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
	for jϵ in 1:h.Nϵ
		l_new = scatter(;x=h.ωgrid, y=ϕc_mat[:,jϵ,1,1,1,1,1,1], line_shape="spline", showlegend=false, marker_color=col[jϵ])
		l[jϵ,1] = l_new
		l_new = scatter(;x=h.ωgrid, y=vf_mat[:,jϵ,1,1,1,1,1,1], line_shape="spline", showlegend=false, marker_color=col[jϵ])
		l[jϵ,2] = l_new
		l_new = scatter(;x=h.ωgrid, y=ϕb_mat[:,jϵ,1,1,1,1,1,1], showlegend=false, marker_color=col[jϵ])
		l[jϵ,3] = l_new
		# l_new = scatter(;x=h.ωgrid, y=ϕb_mat[:,jϵ,1,1,1,1,1,1], showlegend=false, marker_color=col[jϵ])
		l_new = scatter(;x=h.ωgrid, y=θg_mat[:,jϵ,1,1,1,1,1,1], showlegend=false, marker_color=col[jϵ])
		l[jϵ,4] = l_new
	end
	pc = plot([l[jϵ, 1] for jϵ in 1:h.Nϵ], Layout(; font_size=16, title="Consumption"))
	pv = plot([l[jϵ, 2] for jϵ in 1:h.Nϵ], Layout(; font_size=16, title="Value function"))
	pb = plot([l[jϵ, 3] for jϵ in 1:h.Nϵ], Layout(; font_size=16, title="Debt purchases"))
	pθ = plot([l[jϵ, 4] for jϵ in 1:h.Nϵ], Layout(; font_size=16, title="Proportion risk-free debt"))

	p = [pc pv; pb pθ]
	p.plot.layout["xlabel"] = "ω"
	p.plot.layout["width"] = 800
	p.plot.layout["height"] = 600
	p.plot.layout["font_family"] = "Fira Sans Light"

	if remote
		path = pwd() * "/../../Graphs/"
		save(path * "p_hh.jld", "p", p)
	else
		savefig(p, pwd() * "/../Graphs/hh.pdf")
	end

	show_b, show_μ, show_σ, show_w, show_ζ, show_z = mean(h.bgrid), mean(h.μgrid), mean(h.σgrid), mean(h.wgrid), h.ζgrid[2], h.zgrid[1]

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
	for jϵ in 1:h.Nϵ
		l_new = scatter(;x=h.ωgrid, y=ϕc_mat[:,jϵ,1,1,1,1,1,1], line_shape="spline", showlegend=false, marker_color=col[jϵ])
		l[jϵ,1] = l_new
		l_new = scatter(;x=h.ωgrid, y=vf_mat[:,jϵ,1,1,1,1,1,1], line_shape="spline", showlegend=false, marker_color=col[jϵ])
		l[jϵ,2] = l_new
		l_new = scatter(;x=h.ωgrid, y=ϕb_mat[:,jϵ,1,1,1,1,1,1], showlegend=false, marker_color=col[jϵ])
		l[jϵ,3] = l_new
		# l_new = scatter(;x=h.ωgrid, y=ϕb_mat[:,jϵ,1,1,1,1,1,1], showlegend=false, marker_color=col[jϵ])
		l_new = scatter(;x=h.ωgrid, y=θg_mat[:,jϵ,1,1,1,1,1,1], showlegend=false, marker_color=col[jϵ])
		l[jϵ,4] = l_new
	end
	pc = plot([l[jϵ, 1] for jϵ in 1:h.Nϵ], Layout(; font_size=16, title="Consumption"))
	pv = plot([l[jϵ, 2] for jϵ in 1:h.Nϵ], Layout(; font_size=16, title="Value function"))
	pb = plot([l[jϵ, 3] for jϵ in 1:h.Nϵ], Layout(; font_size=16, title="Debt purchases"))
	pθ = plot([l[jϵ, 4] for jϵ in 1:h.Nϵ], Layout(; font_size=16, title="Proportion risk-free debt"))

	p = [pc pv; pb pθ]
	p.plot.layout["xlabel"] = "ω"
	p.plot.layout["width"] = 800
	p.plot.layout["height"] = 600
	p.plot.layout["font_family"] = "Fira Sans Light"

	if remote
		path = pwd() * "/../../Graphs/"
		save(path * "p_hh_def.jld", "p", p)
	else
		path = pwd() * "/../Graphs/"
		savefig(p, path * "hh_def.pdf")
	end


	# leg = Array{String}(1, h.Nω)
	# for jω in 1:h.Nω
	# 	leg[jω] = "ω = $(round(h.ωgrid[jω],2))"
	# end

	# leg = "ω = $(round(mean(h.ωgrid),2))"

	# ϕc_mat = itp_ϕc[h.ωgrid, mean(h.ϵgrid), show_b, show_μ, show_σ, show_w, show_ζ, h.zgrid]
	# ϕa_mat = itp_ϕa[h.ωgrid, mean(h.ϵgrid), show_b, show_μ, show_σ, show_w, show_ζ, h.zgrid]
	# ϕb_mat = itp_ϕb[h.ωgrid, mean(h.ϵgrid), show_b, show_μ, show_σ, show_w, show_ζ, h.zgrid]
	# vf_mat = itp_vf[h.ωgrid, mean(h.ϵgrid), show_b, show_μ, show_σ, show_w, show_ζ, h.zgrid]
    #
	# pc = plot(h.zgrid, ϕc_mat[:,1,1,1,1,1,1,:]', title = "Consumption", label = "")
	# pa = plot(h.zgrid, ϕa_mat[:,1,1,1,1,1,1,:]', title = "Private Savings", label = "")
	# pb = plot(h.zgrid, ϕb_mat[:,1,1,1,1,1,1,:]', title = "Debt Purchases", label = leg)
	# pv = plot(h.zgrid, vf_mat[:,1,1,1,1,1,1,:]', title = "Value Function", label = "")

	# if remote
	path = pwd() * "/../../Graphs/"
	# 	path = path * "Graphs/"
	# 	savefig(path * "hh_z.pdf")
	# else
	# 	path = pwd() * "/../Graphs/"
	# end



	return Void
end

function lines(h::Hank, y, x_dim, name="")
	jshow_b, jshow_μ, jshow_σ, jshow_w, jshow_ζ, jshow_z = ceil(Int, h.Nb/2), ceil(Int, h.Nμ/2), ceil(Int, h.Nσ/2), floor(Int, h.Nw/2), 1, ceil(Int, h.Nz/2)

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

function plot_state_funcs(h::Hank; remote::Bool=false)

	pN_mat = reshape(h.pN,     h.Nb, h.Nμ, h.Nσ, h.Nw, h.Nζ, h.Nz)
	w_mat  = reshape(h.wage,   h.Nb, h.Nμ, h.Nσ, h.Nw, h.Nζ, h.Nz)
	Ld_mat = reshape(h.Ld,     h.Nb, h.Nμ, h.Nσ, h.Nw, h.Nζ, h.Nz)
	Y_mat  = reshape(h.output, h.Nb, h.Nμ, h.Nσ, h.Nw, h.Nζ, h.Nz)

	ppN1 = lines(h, pN_mat, 1, "Price of nontradables")
	pw1  = lines(h, w_mat, 1, "Wage")
	pLd1 = lines(h, Ld_mat, 1, "Labor supply")
	pY1  = lines(h, Y_mat, 1, "Output")

	ppN2 = lines(h, pN_mat, 2)
	pw2  = lines(h, w_mat, 2)
	pLd2 = lines(h, Ld_mat, 2)
	pY2  = lines(h, Y_mat, 2)

	ppN3 = lines(h, pN_mat, 3)
	pw3  = lines(h, w_mat, 3)
	pLd3 = lines(h, Ld_mat, 3)
	pY3  = lines(h, Y_mat, 3)

	ppN4 = lines(h, pN_mat, 4)
	pw4  = lines(h, w_mat, 4)
	pLd4 = lines(h, Ld_mat, 4)
	pY4  = lines(h, Y_mat, 4)

	ppN6 = lines(h, pN_mat, 6)
	pw6  = lines(h, w_mat, 6)
	pLd6 = lines(h, Ld_mat, 6)
	pY6  = lines(h, Y_mat, 6)

	# p = [ppN1 pw1 pLd1; ppN2 pw2 pLd2; ppN3 pw3 pLd3; ppN4 pw4 pLd4; ppN6 pw6 pLd6]
	p = [ppN1 ppN2 ppN3 ppN4 ppN6; pw1 pw2 pw3 pw4 pw6; pLd1 pLd2 pLd3 pLd4 pLd6; pY1 pY2 pY3 pY4 pY6]
	p.plot.layout["width"] = 800
	p.plot.layout["height"] = 640
	p.plot.layout["font_family"] = "Fira Sans Light"
	if remote
		path = pwd() * "/../../Graphs/"
		save(path * "p_statefuncs.jld", "p", p)
	else
		savefig(p, pwd() * "/../Graphs/statefuncs.pdf")
	end
	Void
end

function plot_LoM(h::Hank; remote::Bool=false)
	jz = ceil(Int, h.Nz/2)
	μ′_mat = reshape(h.μ′[:,jz,1], h.Nb, h.Nμ, h.Nσ, h.Nw, h.Nζ, h.Nz)
	σ′_mat = reshape(h.σ′[:,jz,1], h.Nb, h.Nμ, h.Nσ, h.Nw, h.Nζ, h.Nz)

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
	layout = Layout(;	xaxis=attr(title="L", zeroline=true),
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

function plot_nontradables_demand(h::Hank; remote::Bool=false)
	jb = ceil(Int, h.Nb/2)
	jμ = ceil(Int, h.Nμ/2)
	jσ = ceil(Int, h.Nσ/2)
	jw = ceil(Int, h.Nw/2)
	jζ = ceil(Int, h.Nζ/2)
	jz = ceil(Int, h.Nz/2)

	bv, μv, σv, wv, ζv, zv = h.bgrid[jb], h.μgrid[jμ], h.σgrid[jσ], h.wgrid[jw], h.ζgrid[jζ], h.zgrid[jz]

	if remote
		path = pwd() * "/../../Graphs/"
	else
		path = pwd() * "/../Graphs/"
	end
	savefig(path * "N_demand.pdf")
	Void
end


function plot_convergence(dist_statefuncs, dist_LoMs, T::Int64; remote::Bool=false)

	l_funcs1 = scatter(; x=1:T, y = (dist_statefuncs[1:T, 1]), name = "w")
	l_funcs2 = scatter(; x=1:T, y = (dist_statefuncs[1:T, 2]), name = "pN")
	l_funcs3 = scatter(; x=1:T, y = (dist_statefuncs[1:T, 3]), name = "Ld")
	l_LoM1   = scatter(; x=1:T, y = (dist_LoMs[1:T,1]), name = "μ′")
	l_LoM2   = scatter(; x=1:T, y = (dist_LoMs[1:T,2]), name = "σ′")

	layout = Layout(;	xaxis=attr(title="t", zeroline=true),
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


function plot_simul(p::Path; remote::Bool=false)
	T = size(p.data, 1)

	B_vec = series(p,:B)
	μ_vec = series(p,:μ)
	σ_vec = series(p,:σ)
	w_vec = series(p,:w)
	ζ_vec = series(p,:ζ)-1
	z_vec = exp.(series(p,:z))
	Y_vec = series(p,:Y)
	L_vec = series(p,:L)
	π_vec = series(p,:π)
	P_vec = series(p,:P)
	ψ_vec = series(p,:ψ)

	shiftζ = [0; ζ_vec[1:end-1]]

	defaults = find((ζ_vec.==1) .* (shiftζ.==0))./4
	exits    = find((ζ_vec.==0) .* (shiftζ.==1))./4

	times = (1:T)./4

	default_shades = rect(defaults, exits, 0, 1; fillcolor="#d3d3d3", opacity=0.5, line_width=0, xref="x", yref="paper")

	pB = plot(scatter(; x=times, y=B_vec, showlegend=false), Layout(; shapes=default_shades, title="Bonds", xaxis=attr(title="t")));
	pμ = plot(scatter(; x=times, y=μ_vec, showlegend=false), Layout(; shapes=default_shades, title="μ", xaxis=attr(title="t")));
	pσ = plot(scatter(; x=times, y=σ_vec, showlegend=false), Layout(; shapes=default_shades, title="σ", xaxis=attr(title="t")));
	pw = plot(scatter(; x=times, y=w_vec, showlegend=false), Layout(; shapes=default_shades, title="Wage", xaxis=attr(title="t")));
	pζ = plot(scatter(; x=times, y=ζ_vec, showlegend=false), Layout(; shapes=default_shades, title="Default", xaxis=attr(title="t")));
	pz = plot(scatter(; x=times, y=z_vec, showlegend=false), Layout(; shapes=default_shades, title="TFP", xaxis=attr(title="t")));
	pY = plot([ scatter(; x=times, y=Y_vec, showlegend=false);
				scatter(; x=times, y=L_vec, showlegend=false)],
			Layout(; shapes=default_shades, title="Output", xaxis=attr(title="t")));
	pπ = plot(scatter(; x=times, y=π_vec, showlegend=false), Layout(; shapes=default_shades, title="Default prob", xaxis=attr(title="t")));
	pP = plot(scatter(; x=times, y=P_vec, showlegend=false), Layout(; shapes=default_shades, title="Price of nontradables", xaxis=attr(title="t")));
	pψ = plot(scatter(; x=times, y=ψ_vec, showlegend=false), Layout(; shapes=default_shades, title="Fraction domestic", xaxis=attr(title="t")));


	p = [pB pw; pμ pσ; pζ pz]
	# p.plot.layout["shapes"] = default_shades
	p.plot.layout["width"] = 800
	p.plot.layout["height"] = 600
	p.plot.layout["font_family"] = "Fira Sans Light"

	name = "simul"
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
