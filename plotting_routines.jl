abstract type Hank
end

using PlotlyJS, ColorSchemes, ORCA

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

""" Define styles """
def_style = let
	axis = attr(showgrid = true, gridcolor="#e2e2e2", gridwidth=0.5, zeroline=false)
	layout = Layout(xaxis = axis, yaxis=axis)
	Style(layout=layout)
end

slides_def = let
	layout = Layout(plot_bgcolor="#fafafa", paper_bgcolor="#fafafa",
		width=1920*0.45, height=1080*0.45, font_size=16, font_family="Lato",
		legend = attr(orientation = "h", x=0.05))
	Style(def_style, layout=layout)
end

dark_bg = let
	axis = attr(gridcolor="#1b1b1b")
	layout = Layout(plot_bgcolor="#020202", paper_bgcolor="#020202", font_color="white", xaxis=axis,yaxis=axis)
	Style(layout=layout)
end
slides_dark = Style(slides_def, dark_bg)

paper = let
	layout = Layout(width = 1920 * 0.5, height = 1080 * 0.35, font_size=16, font_family = "Linux Libertine",
		legend = attr(orientation = "h", x=0.05))
	Style(def_style, layout=layout)
end

contsty(;div::Bool=false) = let
	colpal = ifelse(div, ColorSchemes.oleron, ColorSchemes.lajolla)
	colscale = [[vv, get(colpal, vv)] for vv in range(0,1,length=100)]
	c_att = attr(colorscale=colscale, autocontour=false)
	Style(trace=Dict(:contour=>c_att))
end

default_eval_points(sd::SOEdef) = floor(Int, N(sd,:b)*0.9), floor(Int, N(sd,:μ)*0.95), max(1,floor(Int, N(sd,:σ)*0.8)), 1, 2, floor(Int, N(sd,:z)*0.15)

function plot_hh_policies(sd::SOEdef; style::Style=slides_def)
	jb, jμ, jσ, jξ, jζ, jz = default_eval_points(sd)
	jb, jμ, jσ, jξ, jζ, jz = 2, 3, 3, 2, 2, 7
	μv, σv = sd.gr[:μ][jμ], sd.gr[:σ][jσ]

	ϕ = Dict(key => [sd.ϕ[key][jω, jϵ, jb, jμ, jσ, jξ, jζ, jz] for jω in 1:N(sd,:ω), jϵ in 1:N(sd,:ϵ)] for key in keys(sd.ϕ))
	ϕ[:v] = [sd.v[:v][jω, jϵ, jb, jμ, jσ, jξ, jζ, jz] for jω in 1:N(sd,:ω), jϵ in 1:N(sd,:ϵ)]

	colpal = ColorSchemes.lajolla

	l = Array{PlotlyBase.GenericTrace}(undef, N(sd,:ϵ), 4)
	for (jϵ, ϵv) in enumerate(sd.gr[:ϵ])
		colv = get(colpal, 0.8*jϵ/N(sd,:ϵ))
		l_new = scatter(;x=sd.gr[:ω], y=ϕ[:c][:,jϵ], xaxis="x1", yaxis="y1", line_shape="spline", name="ϵ = $(round(exp(ϵv),digits=4))", showlegend=false, marker_color=colv)
		l[jϵ,1] = l_new
		l_new = scatter(;x=sd.gr[:ω], y=ϕ[:v][:,jϵ], xaxis="x2", yaxis="y2", line_shape="spline", name="ϵ = $(round(exp(ϵv),digits=4))", showlegend=false, marker_color=colv)
		l[jϵ,2] = l_new
		l_new = scatter(;x=sd.gr[:ω], y=ϕ[:s][:,jϵ], xaxis="x3", yaxis="y3", showlegend=false, name="ϵ = $(round(exp(ϵv),digits=4))", marker_color=colv)
		l[jϵ,3] = l_new
		l_new = scatter(;x=sd.gr[:ω], y=ϕ[:θ][:,jϵ], xaxis="x4", yaxis="y4", showlegend=false, name="ϵ = $(round(exp(ϵv),digits=4))", marker_color=colv)
		l[jϵ,4] = l_new
	end

	ωmax_show = min(sd.pars[:ωmax], quantile(LogNormal(μv, σv), 0.999)+sd.pars[:ωmin])

	# pc = plot([l[jϵ, 1] for jϵ in 1:N(sd,:ϵ)], style=style, Layout(title="<i>Consumption"))
	# pv = plot([l[jϵ, 2] for jϵ in 1:N(sd,:ϵ)], style=style, Layout(title="<i>Value function"))
	# ps = plot([l[jϵ, 3] for jϵ in 1:N(sd,:ϵ)], style=style, Layout(title="<i>Savings"))
	# pθ = plot([l[jϵ, 4] for jϵ in 1:N(sd,:ϵ)], style=style, Layout(title="<i>Proportion risk-free"))
	# p1 = [pc pv; ps pθ]
	# return p1

	shapes = []
	annotations = [
		attr(text = "<i>Consumption", x = 0.45/2, xanchor="center", xref="paper", y = 1, yref="paper", showarrow=false, font_size=18)
		attr(text = "<i>Value function", x = 0.55 + 0.45/2, xanchor="center", xref="paper", y = 1, yref="paper", showarrow=false, font_size=18)
		attr(text = "<i>Savings", x = 0.45/2, xanchor="center", xref="paper", y = 0.435, yref="paper", showarrow=false, font_size=18)
		attr(text = "<i>Proportion risk-free", x = 0.55 + 0.45/2, xanchor="center", xref="paper", y = 0.435, yref="paper", showarrow=false, font_size=18)
		]

	layout = Layout(shapes = shapes, annotations=annotations,
		height = 1080*0.55,
		xaxis1 = attr(domain=[0, 0.425], anchor="y1", zeroline=false),
		xaxis2 = attr(domain=[0.525, 0.975], anchor="y2", zeroline=false),
		xaxis3 = attr(domain=[0, 0.425], anchor="y3", zeroline=false),
		xaxis4 = attr(domain=[0.525, 0.975], anchor="y4", zeroline=false),
		yaxis1 = attr(domain=[0.525, 0.975], anchor="x1", zeroline=false),
		yaxis2 = attr(domain=[0.525, 0.975], anchor="x2", zeroline=false),
		yaxis3 = attr(domain=[0, 0.425], anchor="x3", zeroline=false),
		yaxis4 = attr(domain=[0, 0.425], anchor="x4", zeroline=false),
		)

	plot(l[:], layout, style=style)
end	

function makecontour(sd::SOEdef, y::Matrix, dim1::Symbol, dim2::Symbol, min_z, max_z; f1::Function=identity, f2::Function=identity, divergent::Bool=false, reversescale::Bool=false, suffix="")

	if divergent
		colpal = ColorSchemes.broc
	else
		colpal = ColorSchemes.lajolla
	end

	colscale = [[vv, get(colpal, vv)] for vv in range(0,1,length=100)]

	contour(;x = f1.(sd.gr[dim1]), y = f2.(sd.gr[dim2]), z=y, colorscale = colscale, reversescale=reversescale, autocontour=false, colorbar_ticksuffix=suffix, colorbar_showticksuffix="all", contours=Dict(:start=>min_z,:end=>max_z), xaxis="x1", yaxis="y1")
end

function makecontour_μσ(sd::SOEdef, y::Matrix, min_z, max_z; divergent::Bool=false, reversescale::Bool=false, fz::Function=identity, suffix="")
	knots = (sd.gr[:μ], sd.gr[:σ])
	itp_y = extrapolate(interpolate(knots, y, Gridded(Linear())), Interpolations.Line())

	y_mat, xgrid, ygrid = reeval_mat_MV(sd, itp_y, lb = 0)
	xax, yax = "Mean", "Variance"

	if divergent
		colpal = ColorSchemes.broc
	else
		colpal = ColorSchemes.lajolla
	end
	colscale = [[vv, get(colpal, vv)] for vv in range(0,1,length=100)]

	ctμσ = contour(;
	x = xgrid, y = ygrid,
	z = fz.(y_mat), xaxis="x2", yaxis="y2",
	colorscale = colscale, reversescale=reversescale, colorbar_ticksuffix=suffix, colorbar_showticksuffix="all", autocontour=false, contours=Dict(:start=>min_z,:end=>max_z))
end

function reeval_mat_MV(sd::SOEdef, itp_obj; lb=-Inf, ub=Inf)
	lb < ub || throw(error("Must specify upper bound greater than lower bound"))
	m_min, v_min = unmake_logN(sd.gr[:μ][1], sd.gr[:σ][1])
	m_max, v_max = unmake_logN(sd.gr[:μ][3], sd.gr[:σ][4])

	# itp_obj = extrapolate(itp_obj, Interpolations.Flat())
	
	Nm = max(4*N(sd,:μ), 4*N(sd,:σ))

	mgrid = range(m_min, m_max, length=Nm)
	vgrid = range(v_min, v_max, length=Nm)
	
	mat = zeros(Nm, Nm)
	for (jm, m) in enumerate(mgrid)
		for (jv, v) in enumerate(vgrid)
			μv, σv = make_logN(m, v)
			Y = itp_obj(μv, σv)
			mat[jm, jv] = max(lb, min(ub, Y))
		end
	end
	
	return mat, mgrid, vgrid
end


function make_unemp(sd::SOEdef; style::Style=slides_def, leg=true)
	jb, jμ, jσ, jξ, jζ, jz = [10, 3,2,1,2,5]

	unemp = (1 .- reshape_long(sd, sd.eq[:Ld])) * 100

	U_matbz = [unemp[jb1, jμ, jσ, jξ, jζ, jz1] for (jb1, bv) in enumerate(sd.gr[:b]), (jz1,zv) in enumerate(sd.gr[:z])]

	U_matμσ = [unemp[jb, jμ, jσ, jξ, jζ, jz] for (jμ, μv) in enumerate(sd.gr[:μ]), (jσ,σv) in enumerate(sd.gr[:σ])]

	min_z1, max_z1 = extrema(U_matbz)
	min_z2, max_z2 = extrema(U_matμσ)
	min_z = min(min_z1, min_z2)
	max_z = max(max_z1, max_z2)

	data1 = makecontour(sd, U_matbz, :b, :z, min_z, max_z, f2=x->100x, suffix="%")

	data2 = makecontour_μσ(sd, U_matμσ, min_z, max_z, suffix="%")

	data = [data1, data2]
	# data = data1

	layout = Layout(title=ifelse(leg,"<i>Unemployment", ""),
		xaxis1 = attr(domain=[0, 0.45], anchor="y1", title="<i>B"),
		xaxis2 = attr(domain=[0.55, 1], anchor="y2", title="<i>Mean"),
		yaxis1 = attr(anchor="x1", title="<i>z"),
		yaxis2 = attr(anchor="x2", title="<i>Variance"),
		)
	plot(data, layout, style=style)
	# data1
end

function make_debtprice(sd::SOEdef; style::Style=slides_def, leg=true)
	jb, jμ, jσ, jξ, jζ, jz = [7, 4,3,1,2,5]

	qg_matbz = [reshape_long(sd, sd.eq[:qᵍ])[jb, jμ, jσ, jξ, jζ, jz] for (jb, bv) in enumerate(sd.gr[:b]), (jz,zv) in enumerate(sd.gr[:z])]
	qg_matμσ = [reshape_long(sd, sd.eq[:qᵍ])[jb, jμ, jσ, jξ, jζ, jz] for (jμ, μv) in enumerate(sd.gr[:μ]), (jσ,σv) in enumerate(sd.gr[:σ])]
	

	data1 = makecontour(sd, qg_matbz, :b, :z, 0.5, 1, f2=x->100x, reversescale=true)
	data2 = makecontour_μσ(sd, qg_matμσ, 0.5, 1, reversescale=true)

	data = [data1, data2]
	# data = data1

	layout = Layout(title=ifelse(leg,"<i>Price of Debt", ""),
		xaxis1 = attr(domain=[0, 0.45], anchor="y1", title="<i>B"),
		xaxis2 = attr(domain=[0.55, 1], anchor="y2", title="<i>Mean"),
		yaxis1 = attr(anchor="x1", title="<i>z"),
		yaxis2 = attr(anchor="x2", title="<i>Variance"),
		)

	plot(data, layout, style=style)
end

function make_def_incentive(sd::SOEdef; style::Style=slides_def)
	jb, jμ, jσ, jξ, jζ, jz = default_eval_points(sd)
	μv = sd.gr[:μ][jμ]
	σv = sd.gr[:σ][jσ]
	ξv = sd.gr[:ξ][jξ]

	vR = [reshape_long(sd, sd.eq[:welfare])[jb, jμ, jσ, jξ, 2, jz] for (jb, bv) in enumerate(sd.gr[:b]), (jz,zv) in enumerate(sd.gr[:z])]
	itp_vD = make_itp(sd, sd.eq[:welfare], agg=true)
	vD = [itp_vD(1-sd.pars[:ℏ]*bv, μv, σv, ξv, 0.0, zv) for (jb, bv) in enumerate(sd.gr[:b]), (jz,zv) in enumerate(sd.gr[:z])]

	p1 = makecontour(sd, vR-vD, :b, :z, f2=x->100x, style=style, title="<i>V<sup>R</sup> - V<sup>D</sup>", reversescale=false)
end

scats_crises(pv::Vector{T}, tvv::Vector{Vector{Int64}}, key::Symbol, X::Float64) where T <: AbstractPath = scats_crises(pv, tvv, key, x->x/X)
function scats_crises(pv::Vector{T}, tvv::Vector{Vector{Int64}}, key::Symbol, f::Function=identity; axis::Int64=1, indiv=false, k=8) where T <: AbstractPath

	ymat, y_up, y_me, y_lo, y_av = series_crises(pv, tvv, key, k)

	line_up = scatter(x=(-2k:k)/4, y=f(y_up),hoverinfo="skip",showlegend=false,legendgroup=3,mode="lines",line=attr(color="rgb(31,119,180)", width=0.001), xaxis="x$axis", yaxis="y$axis")
	line_lo = scatter(x=(-2k:k)/4, y=f(y_lo),hoverinfo="skip",showlegend=(axis==1),legendgroup=3, name="IQR",mode="lines",line=attr(color="rgb(31,119,180)",width=0.001), fill="tonexty", fillcolor="rgba(31,119,180,0.25)", xaxis="x$axis", yaxis="y$axis")
	line_me = scatter(x=(-2k:k)/4, y=f(y_me), name="Median",legendgroup=1,showlegend=(axis==1), mode="lines", line=attr(color=col[4]), xaxis="x$axis", yaxis="y$axis")
	line_av = scatter(x=(-2k:k)/4, y=f(y_av), name="Mean",legendgroup=2,showlegend=(axis==1), mode="lines", line=attr(color=col[3], dash="dash"), xaxis="x$axis", yaxis="y$axis")


	s1 = [line_up, line_lo, line_me, line_av]
	if indiv
		s1 =[
		[scatter(x=(-2k:k)/4, y=f.(ymat)[:,jc], name="Individual", legendgroup=3, mode="lines", line_color="rgb(31,119,180)", opacity=0.75, line_width=0.5, showlegend=((axis==1)&(jc==1)), xaxis="x$axis", yaxis="y$axis", hoverinfo="skip") for jc in 1:size(ymat,2)]
		line_up
		line_lo
		line_me
		line_av]
	end
	s1
end

function plot_crises(pv::Vector{T}, πthres::Float64, key::Symbol, f::Function=identity; style::Style=slides_def) where T <: AbstractPath
	Nc, tvv = get_crises(pv, πthres, 8)
	shapes = [vline(0, line_width=1, marker_color="black")]
	scats_crises(pv, tvv, key, f)
	plot(scats_crises(pv, tvv, key, f), style=style, Layout(shapes=shapes))
end

scats_comp(pv_bench::Vector{T}, pv_nodef::Vector{T}, tvv::Vector{Vector{Int64}}, key::Symbol, x1::Float64, x2::Float64=x1; CI::Bool=false, avg::Bool=false) where T <: AbstractPath = scats_comp(pv_bench, pv_nodef, tvv, key, x->x/x1, x->x/x2, CI=CI, avg=avg)
function scats_comp(pvb::Vector{T}, pvn::Vector{T}, tvv::Vector{Vector{Int64}}, key::Symbol, f1::Function=identity, f2::Function=f1; CI::Bool=false, avg::Bool=false, axis::Int64=1, k=8) where T <: AbstractPath
	ybench, bench_up, bench_me, bench_lo, bench_av = series_crises(pvb, tvv, key, k)
	ynodef, nodef_up, nodef_me, nodef_lo, nodef_av = series_crises(pvn, tvv, key, k)

	colbench = get(ColorSchemes.vik, 0.3)
	colnodef = get(ColorSchemes.vik, 0.7)

	line_bench = scatter(x=(-2k:k)/4, y=f1.(bench_me), name="Benchmark", mode="lines", line_color=colbench, fill="tonexty", showlegend=(axis==1), legendgroup = 1, xaxis="x$axis", yaxis="y$axis")
	# lb_avg = scatter(x=(-2k:k)/4, y=f1.(bench_av), name="Benchmark", mode="lines", line_color=col[fill="tonexty", 1], showlegend=(axis==1), legendgroup = 1, xaxis="x$axis", yaxis="y$axis")
	lb_up = scatter(x=(-2k:k)/4, y=f1.(bench_up), hoverinfo="skip", mode="lines", line_width=0.001, line_color=colbench, showlegend=false, legendgroup = 1, xaxis="x$axis", yaxis="y$axis")
	lb_lo = scatter(x=(-2k:k)/4, y=f1.(bench_lo), hoverinfo="skip", mode="lines", line_width=0.001, line_color=colbench, fill="tonexty", showlegend=false, legendgroup = 1, xaxis="x$axis", yaxis="y$axis")
	line_nodef = scatter(x=(-2k:k)/4, y=f2.(nodef_me), name="No default", line_dash="dashdot", mode="lines", line_color=colnodef, fill="tonexty", showlegend=(axis==1), legendgroup = 2, xaxis="x$axis", yaxis="y$axis")
	# ln_avg = scatter(x=(-2k:k)/4, y=f2.(nodef_av), name="No default", mode="lines", line_color=colnodef, fill="tonexty", showlegend=(axis==1), legendgroup = 2, xaxis="x$axis", yaxis="y$axis")
	ln_up = scatter(x=(-2k:k)/4, y=f2.(nodef_up), hoverinfo="skip", mode="lines", line_width=0.001, line_color=colnodef, showlegend=false, legendgroup = 2, xaxis="x$axis", yaxis="y$axis")
	ln_lo = scatter(x=(-2k:k)/4, y=f2.(nodef_lo), hoverinfo="skip", mode="lines", line_width=0.001, line_color=colnodef, fill="tonexty", showlegend=false, legendgroup = 2, xaxis="x$axis", yaxis="y$axis")

	s1 = [line_bench, line_nodef]
	if avg
		s1 = [lb_avg, ln_avg]
	elseif CI
		s1 = [lb_up, line_bench, lb_lo, ln_up, line_nodef, ln_lo]
	end
	s1
end


function scats_MIT(pvb::Vector{T}, pvh::Vector{T}, key::Symbol, f::Function, g::Function; axis=1) where T<:AbstractPath
	xvec = (1:periods(first(pvb))) / 4
	ybmat = [g(series(pvb[jk], key)[jt]) for jt in 1:length(xvec), jk in 1:length(pvb)]
	yhmat = [g(series(pvh[jk], key)[jt]) for jt in 1:length(xvec), jk in 1:length(pvh)]

	ymat = [f(ybmat[jt,jk], yhmat[jt,jk]) for jt in 1:length(xvec), jk in 1:length(pvh)]

	[
	scatter(x=xvec, y=[quantile(ybmat[jt,:], 0.75) for jt in 1:length(xvec)], name="Benchmark", legendgroup=1, line_color=col[1], showlegend=false, xaxis="x$axis", yaxis="y$axis", hoverinfo="skip", line_dash="dash")
	scatter(x=xvec, y=[quantile(ybmat[jt,:], 0.25) for jt in 1:length(xvec)], name="Benchmark", legendgroup=1, line_color=col[1], showlegend=false, xaxis="x$axis", yaxis="y$axis", hoverinfo="skip", line_dash="dash")
	scatter(x=xvec, y=[mean(ybmat[jt,:]) for jt in 1:length(xvec)], name="Benchmark", legendgroup=1, showlegend=(axis==1), xaxis="x$axis", yaxis="y$axis", line_color=col[1])

	scatter(x=xvec, y=[quantile(yhmat[jt,:], 0.75) for jt in 1:length(xvec)], name="Shock", legendgroup=2, line_color=col[2], showlegend=false, xaxis="x$axis", yaxis="y$axis", hoverinfo="skip", line_dash="dash")
	scatter(x=xvec, y=[quantile(yhmat[jt,:], 0.25) for jt in 1:length(xvec)], name="Shock", legendgroup=2, line_color=col[2], showlegend=false, xaxis="x$axis", yaxis="y$axis", hoverinfo="skip", line_dash="dash")
	scatter(x=xvec, y=[mean(yhmat[jt,:]) for jt in 1:length(xvec)], name="Shock", legendgroup=2, showlegend=(axis==1), xaxis="x$axis", yaxis="y$axis", line_color=col[2])
	]
end

function make_MIT_shock(sd::SOEdef, B0 = mean(sd.gr[:b]), ϵb = 0.05; K=100, T=4*5, burn_in = 4*100, yh = 0.65, verbose=false, nodef=true, style::Style=slides_def)

	pvb, pvh = MIT_shock(sd, B0, ϵb; K=K, T=T, burn_in=burn_in, verbose=verbose)

	if nodef
		indices = [jj for jj in 1:length(pvb) if minimum(series(pvb[jj],:ζ))*minimum(series(pvh[jj],:ζ)) == 1]

		pvb = pvb[indices]
		pvh = pvh[indices]
	end

	keyvec = [:z, :Y, :C, :B, :G, :T, :L, :π, :Wr]
	titlevec = ["TFP", "Output", "Consumption", "Bonds", "Govt spending", "Lump-sum taxes", "Unemployment", "Default prob", "Welfare in repayment"]

	fvec = Vector{Function}(undef, length(keyvec))
	fvec .= ((x,y) -> y/x - 1)
	gvec = Function[identity for jj in 1:length(keyvec)]
	gvec[7] = x->100*(1-x)
	gvec[[1,8]] .= x->100x
	ytitle = ["" for jj in 1:length(keyvec)]

	fvec[[1,6,7,8,9]] .= (x,y) -> y-x

	fvec .= (x,y) -> x
	# fvec[[2, 5, 6]] .= (x->100*x/meanYb)
	# fvec[3] = (x->100*x/meanCb)
	# fvec[[4]] .= (x->25 * x/meanYb)
	# fvec[7] = (x->100*(1 .- x))
	# ytitle[[4,5,6]] .= "% of mean GDP"
	# ytitle[[2, 3]] .= "% dev from mean"
	ytitle[[1,7,8]] .= "%"

	data = Vector{GenericTrace{Dict{Symbol,Any}}}(undef, 0)
	for (jj, key) in enumerate(keyvec)
		println(key)
		for scat in scats_MIT(pvb, pvh, key, fvec[jj], gvec[jj], axis=jj)
			push!(data, scat)
		end
	end

	a = 1/3
	b = 1/20
	bx = 1/20
	ys = [1, 0.64, 0.3]
	annotations = [
		attr(text=titlevec[jj], x = 2.5, xanchor="center", xref = "x$jj", y = ys[ceil(Int, jj/3)], showarrow=false, font_size=18, yref="paper") for jj in 1:length(titlevec)
		]

	layout = Layout(annotations = annotations,
		height = 1080*yh, width = 1920*0.65, legend = attr(y=0, yref="paper", x=0.5, xanchor="center", xref="paper"),
		xaxis1 = attr(domain = [0a+bx, a-bx], anchor="y1"),
		xaxis2 = attr(domain = [1a+bx, 2a-bx], anchor="y2"),
		xaxis3 = attr(domain = [2a+bx, 3a-bx], anchor="y3"),
		yaxis1 = attr(anchor = "x1", domain = [2a+b, 3a-b], titlefont_size = 14, title=ytitle[1]),
		yaxis2 = attr(anchor = "x2", domain = [2a+b, 3a-b], titlefont_size = 14, title=ytitle[2]),
		yaxis3 = attr(anchor = "x3", domain = [2a+b, 3a-b], titlefont_size = 14, title=ytitle[3]),
		xaxis4 = attr(domain = [0a+bx, a-bx], anchor="y4"),
		xaxis5 = attr(domain = [1a+bx, 2a-bx], anchor="y5"),
		xaxis6 = attr(domain = [2a+bx, 3a-bx], anchor="y6"),
		yaxis4 = attr(anchor = "x4", domain = [1a+b, 2a-b], titlefont_size = 14, title=ytitle[4]),
		yaxis5 = attr(anchor = "x5", domain = [1a+b, 2a-b], titlefont_size = 14, title=ytitle[5]),
		yaxis6 = attr(anchor = "x6", domain = [1a+b, 2a-b], titlefont_size = 14, title=ytitle[6]),
		xaxis7 = attr(domain = [0a+bx, a-bx], anchor="y7"),
		xaxis8 = attr(domain = [1a+bx, 2a-bx], anchor="y8"),
		xaxis9 = attr(domain = [2a+bx, 3a-bx], anchor="y9"),
		yaxis7 = attr(anchor = "x7", domain = [0a+b, 1a-b], titlefont_size = 14, title=ytitle[7]),
		yaxis8 = attr(anchor = "x8", domain = [0a+b, 1a-b], titlefont_size = 14, title=ytitle[8]),
		yaxis9 = attr(anchor = "x9", domain = [0a+b, 1a-b], titlefont_size = 14, title=ytitle[9]),
	)
	data
	plot(data, layout, style=style)
end

panels_defaults(pv::Vector{T}; k=8, style::Style=slides_def, yh = 0.65, indiv=false) where T<:AbstractPath = panels_crises(pv, 0.0, style=style, yh=yh, type="default", indiv=indiv, k=k)
function panels_crises(pv::Vector{T}, πthres::Float64; style::Style=slides_def, yh = 0.65, type="highspreads", indiv=false, k=8) where T<:AbstractPath
	Nc, tvv = get_crises(pv, πthres, k, type=type)
	println("Suggested yh=0.7 for style=paper")
	keyvec = [:z, :Y, :C, :CoY, :B, :ψ, :qg, :π, :L, :mean, :var, :P, :avgω, :p90, :G, :T]

	titlevec = ["TFP", "Output", "Consumption", "<i>C/Y<sup>d</sup>", "Bonds", "Proportion Domestic", "Price of new debt", "Default prob", "Unemployment", "Wealth Dist Mean", "Wealth Dist Variance", "Price of nontradables", "Bondholdings-weighted avg wealth", "top 10% holdings", "Govt spending", "Lump-sum taxes"]

	meanY = mean([mean(series(p, :Y)) for p in pv])
	meanC = mean([mean(series(p, :C)) for p in pv])
	
	funcvec = Vector{Function}(undef, length(keyvec))
	funcvec .= identity
	ytitle = ["" for jj in 1:length(keyvec)]

	funcvec[[2, 15, 16]] .= (x->100*x/meanY)
	funcvec[3] = (x->100*x/meanC)
	funcvec[[5, 10, 13]] .= (x->25 * x/meanY)
	funcvec[[1, 4, 6, 8]] .= (x->100*x)
	funcvec[9] = (x->100*(1 .- x))
	ytitle[[5, 10, 13, 15, 16]] .= "% of mean GDP"
	ytitle[[2, 3]] .= "% dev from mean"
	ytitle[[1,4,6,8,9]] .= "%"

	if haskey(first(pv).data, :Gini)
		titlevec[keyvec.==:var] .= "Wealth Gini"
		funcvec[keyvec.==:var] .= x->100*x
		keyvec[keyvec.==:var] .= :Gini
	end

	
	data = Vector{GenericTrace{Dict{Symbol,Any}}}(undef, 0)
	for (jj, key) in enumerate(keyvec)
		for scat in scats_crises(pv, tvv, key, funcvec[jj], axis=jj, indiv=indiv, k=k)
			push!(data, scat)
		end
	end

	a = 1/4
	b = 1/20
	bx = 1/40
	shapes = [
		# vline(0, line_width=1, marker_color="black")
		]

	ys = [1, 0.75, 0.475, 0.21]
	annotations = [
		attr(text=titlevec[jj], x = -k/2+(3k)/8, xanchor="center", xref = "x$jj", y = ys[ceil(Int, jj/4)], font_size=18, showarrow=false, yref="paper") for jj in 1:length(titlevec)
		]

	layout = Layout(shapes=shapes, annotations = annotations,
		height = 1080*yh, width = 1920*0.65, legend = attr(y=0, yref="paper", x=0.5, xanchor="center", xref="paper"),
		xaxis1 = attr(domain = [0a, a-2bx], range=[-k/2, k/4]),
		xaxis2 = attr(domain = [1a+bx, 2a-bx], range=[-k/2, k/4]),
		xaxis3 = attr(domain = [2a+bx, 3a-bx], range=[-k/2, k/4]),
		xaxis4 = attr(domain = [3a+2bx, 4a], range=[-k/2, k/4]),
		yaxis1 = attr(anchor = "x1", domain = [3a+b, 4a-b], titlefont_size = 16, title=ytitle[1]),
		yaxis2 = attr(anchor = "x2", domain = [3a+b, 4a-b], titlefont_size = 16, title=ytitle[2]),
		yaxis3 = attr(anchor = "x3", domain = [3a+b, 4a-b], titlefont_size = 16, title=ytitle[3]),
		yaxis4 = attr(anchor = "x4", domain = [3a+b, 4a-b], titlefont_size = 16, title=ytitle[4]),
		xaxis5 = attr(domain = [0a, a-2bx], anchor="y5", range=[-k/2, k/4]),
		xaxis6 = attr(domain = [1a+bx, 2a-bx], anchor="y6", range=[-k/2, k/4]),
		xaxis7 = attr(domain = [2a+bx, 3a-bx], anchor="y7", range=[-k/2, k/4]),
		xaxis8 = attr(domain = [3a+2bx, 4a], anchor="y8", range=[-k/2, k/4]),
		yaxis5 = attr(anchor = "x5", domain = [2a+b, 3a-b], titlefont_size = 16, title=ytitle[5]),
		yaxis6 = attr(anchor = "x6", domain = [2a+b, 3a-b], titlefont_size = 16, title=ytitle[6]),
		yaxis7 = attr(anchor = "x7", domain = [2a+b, 3a-b], titlefont_size = 16, title=ytitle[7]),
		yaxis8 = attr(anchor = "x8", domain = [2a+b, 3a-b], titlefont_size = 16, title=ytitle[8]),
		xaxis9 = attr(domain = [0a, a-2bx], anchor="y9", range=[-k/2, k/4]),
		xaxis10 = attr(domain = [1a+bx, 2a-bx], anchor="y10", range=[-k/2, k/4]),
		xaxis11 = attr(domain = [2a+bx, 3a-bx], anchor="y11", range=[-k/2, k/4]),
		xaxis12 = attr(domain = [3a+2bx, 4a], anchor="y12", range=[-k/2, k/4]),
		yaxis9 = attr(anchor = "x9", domain = [1a+b, 2a-b], titlefont_size = 16, title=ytitle[9]),
		yaxis10 = attr(anchor = "x10", domain = [1a+b, 2a-b], titlefont_size = 16, title=ytitle[10]),
		yaxis11 = attr(anchor = "x11", domain = [1a+b, 2a-b], titlefont_size = 16, title=ytitle[11]),
		yaxis12 = attr(anchor = "x12", domain = [1a+b, 2a-b], titlefont_size = 16, title=ytitle[12]),
		xaxis13 = attr(domain = [0a, a-2bx], anchor="y13", range=[-k/2, k/4]),
		xaxis14 = attr(domain = [1a+bx, 2a-bx], anchor="y14", range=[-k/2, k/4]),
		xaxis15 = attr(domain = [2a+bx, 3a-bx], anchor="y15", range=[-k/2, k/4]),
		xaxis16 = attr(domain = [3a+2bx, 4a], anchor="y16", range=[-k/2, k/4]),
		yaxis13 = attr(anchor = "x13", domain = [0a+b, a-b], titlefont_size = 16, title=ytitle[13]),
		yaxis14 = attr(anchor = "x14", domain = [0a+b, a-b], titlefont_size = 16, title=ytitle[14]),
		yaxis15 = attr(anchor = "x15", domain = [0a+b, a-b], titlefont_size = 16, title=ytitle[15]),
		yaxis16 = attr(anchor = "x16", domain = [0a+b, a-b], titlefont_size = 16, title=ytitle[16]),
		)

	plot(data, layout, style=style)
end

function panels_comp(pv_bench::Vector{T}, pv_nodef::Vector{T}, πthres::Float64; style::Style=slides_def, yh = 0.65, k=8) where T<:AbstractPath
	Nc, tvv = get_crises(pv_bench, πthres, k)
	println("Suggested yh=0.7 for style=paper")
	keyvec = [:z, :Y, :C, :B, :G, :T, :L, :qg, :Wr]

	titlevec = ["TFP", "Output", "Consumption", "Bonds", "Govt spending", "Lump-sum taxes", "Unemployment", "Price of new debt", "Welfare in repayment"]

	meanYb = mean([mean(series(p, :Y)) for p in pv_bench])
	meanCb = mean([mean(series(p, :C)) for p in pv_bench])
	meanYn = mean([mean(series(p, :Y)) for p in pv_nodef])
	meanCn = mean([mean(series(p, :C)) for p in pv_nodef])
	
	f1vec = Vector{Function}(undef, length(keyvec))
	f1vec .= identity
	ytitle = ["" for jj in 1:length(keyvec)]

	f1vec[1] = x->100*x
	f1vec[[2, 5, 6]] .= (x->100*x/meanYb)
	f1vec[3] = (x->100*x/meanCb)
	f1vec[[4]] .= (x->25 * x/meanYb)
	f1vec[7] = (x->100*(1 .- x))
	ytitle[[4,5,6]] .= "% of mean GDP"
	ytitle[[2, 3]] .= "% dev from mean"
	ytitle[[1,7]] .= "%"

	f2vec = copy(f1vec)
	f2vec[[2,5,6]] .= (x->100*x/meanYn)
	f2vec[3] = (x->100*x/meanCn)
	f2vec[[4]] .= (x->25 * x/meanYn)

	data = Vector{GenericTrace{Dict{Symbol,Any}}}(undef, 0)
	for (jj, key) in enumerate(keyvec)
		for scat in scats_comp(pv_bench, pv_nodef, tvv, key, f1vec[jj], f2vec[jj], axis=jj, CI=true, k=k)
			push!(data, scat)
		end
	end

	a = 1/3
	b = 1/20
	bx = 1/30
	shapes = [
		# vline(0, line_width=1, marker_color="black")
		]

	ys = [1, 0.64, 0.3]
	annotations = [
		attr(text=titlevec[jj], x = -1, xanchor="center", xref = "x$jj", y = ys[ceil(Int, jj/3)], showarrow=false, font_size = 18, yref="paper") for jj in 1:length(titlevec)
		]

	layout = Layout(shapes=shapes, annotations = annotations,
		height = 1080*yh, width = 1920*0.65, legend = attr(y=0, yref="paper", x=0.5, xanchor="center", xref="paper"),
		xaxis1 = attr(domain = [0a, a-2bx], anchor="y1"),
		xaxis2 = attr(domain = [1a+bx, 2a-bx], anchor="y2"),
		xaxis3 = attr(domain = [2a+2bx, 3a], anchor="y3"),
		yaxis1 = attr(anchor = "x1", domain = [2a+b, 3a-b], titlefont_size = 14, title=ytitle[1]),
		yaxis2 = attr(anchor = "x2", domain = [2a+b, 3a-b], titlefont_size = 14, title=ytitle[2]),
		yaxis3 = attr(anchor = "x3", domain = [2a+b, 3a-b], titlefont_size = 14, title=ytitle[3]),
		xaxis4 = attr(domain = [0a, a-2bx], anchor="y4"),
		xaxis5 = attr(domain = [1a+bx, 2a-bx], anchor="y5"),
		xaxis6 = attr(domain = [2a+2bx, 3a], anchor="y6"),
		yaxis4 = attr(anchor = "x4", domain = [1a+b, 2a-b], titlefont_size = 14, title=ytitle[4]),
		yaxis5 = attr(anchor = "x5", domain = [1a+b, 2a-b], titlefont_size = 14, title=ytitle[5]),
		yaxis6 = attr(anchor = "x6", domain = [1a+b, 2a-b], titlefont_size = 14, title=ytitle[6]),
		xaxis7 = attr(domain = [0a, a-2bx], anchor="y7"),
		xaxis8 = attr(domain = [1a+bx, 2a-bx], anchor="y8"),
		xaxis9 = attr(domain = [2a+2bx, 3a], anchor="y9"),
		yaxis7 = attr(anchor = "x7", domain = [0a+b, 1a-b], titlefont_size = 14, title=ytitle[7]),
		yaxis8 = attr(anchor = "x8", domain = [0a+b, 1a-b], titlefont_size = 14, title=ytitle[8]),
		yaxis9 = attr(anchor = "x9", domain = [0a+b, 1a-b], titlefont_size = 14, title=ytitle[9]),
	)
	plot(data, layout, style=style)
end

function plot_comp(pv_bench::Vector{T}, pv_nodef::Vector{T}, πthres, key::Symbol, f::Function=identity; style::Style=slides_def) where T <: AbstractPath
	Nc, tvv = get_crises(pv_bench, πthres, 8)
	shapes = [vline(0, line_width=1, marker_color="black")]
	scats_comp(pv_bench, pv_nodef, tvv, key, f)
	plot(scats_comp(pv_bench, pv_nodef, tvv, key, f), style=style, Layout(shapes=shapes))
end

function panels_full_comp(pv_bench::Vector{T}, pv_noΔ::Vector{T}, pv_nob::Vector{T}, pv_nodef::Vector{T}, πthres::Float64; style::Style=slides_def, yh = 0.65) where T<:AbstractPath
	Nc, tvv = get_crises(pv_bench, πthres, 8)
	println("Suggested yh=0.7 for style=paper")
	keyvec = [:z, :Y, :C, :B, :G, :T, :L, :qg, :Wr]

	titlevec = ["TFP", "Output", "Consumption", "Bonds", "Govt spending", "Lump-sum taxes", "Unemployment", "Price of new debt", "Welfare in repayment"]

	meanYb = mean([mean(series(p, :Y)) for p in pv_bench])
	meanCb = mean([mean(series(p, :C)) for p in pv_bench])
	meanYΔ = mean([mean(series(p, :Y)) for p in pv_noΔ])
	meanCΔ = mean([mean(series(p, :C)) for p in pv_noΔ])
	meanYb = mean([mean(series(p, :Y)) for p in pv_nob])
	meanCb = mean([mean(series(p, :C)) for p in pv_nob])
	meanYn = mean([mean(series(p, :Y)) for p in pv_nodef])
	meanCn = mean([mean(series(p, :C)) for p in pv_nodef])
	
	f1vec = Vector{Function}(undef, length(keyvec))
	f1vec .= identity
	ytitle = ["" for jj in 1:length(keyvec)]

	f1vec[1] = x->100*x
	f1vec[[2, 5, 6]] .= (x->100*x/meanYb)
	f1vec[3] = (x->100*x/meanCb)
	f1vec[[4]] .= (x->25 * x/meanYb)
	f1vec[7] = (x->100*(1 .- x))
	ytitle[[4,5,6]] .= "% of mean GDP"
	ytitle[[2, 3]] .= "% dev from mean"
	ytitle[[1,7]] .= "%"

	f2vec = copy(f1vec)
	f2vec[[2,5,6]] .= (x->100*x/meanYΔ)
	f2vec[3] = (x->100*x/meanCΔ)
	f2vec[[4]] .= (x->25 * x/meanYΔ)

	f3vec = copy(f1vec)
	f3vec[[2,5,6]] .= (x->100*x/meanYb)
	f3vec[3] = (x->100*x/meanCb)
	f3vec[[4]] .= (x->25 * x/meanYb)

	f4vec = copy(f1vec)
	f4vec[[2,5,6]] .= (x->100*x/meanYn)
	f4vec[3] = (x->100*x/meanCn)
	f4vec[[4]] .= (x->25 * x/meanYn)

	data = Vector{GenericTrace{Dict{Symbol,Any}}}(undef, 0)
	for (jj, key) in enumerate(keyvec)
		for scat in full_scats_comp(pv_bench, pv_noΔ, pv_nob, pv_nodef, tvv, key, f1vec[jj], f2vec[jj], f3vec[jj], f4vec[jj], axis=jj)
			push!(data, scat)
		end
	end

	a = 1/3
	b = 1/20
	bx = 1/30
	shapes = [
		# vline(0, line_width=1, marker_color="black")
		]

	ys = [1, 0.64, 0.3]
	annotations = [
		attr(text=titlevec[jj], x = -1, xanchor="center", xref = "x$jj", y = ys[ceil(Int, jj/3)], showarrow=false, font_size = 18, yref="paper") for jj in 1:length(titlevec)
		]

	layout = Layout(shapes=shapes, annotations = annotations,
		height = 1080*yh, width = 1920*0.65, legend = attr(y=0, yref="paper", x=0.5, xanchor="center", xref="paper"),
		xaxis1 = attr(domain = [0a, a-2bx], anchor="y1"),
		xaxis2 = attr(domain = [1a+bx, 2a-bx], anchor="y2"),
		xaxis3 = attr(domain = [2a+2bx, 3a], anchor="y3"),
		yaxis1 = attr(anchor = "x1", domain = [2a+b, 3a-b], titlefont_size = 14, title=ytitle[1]),
		yaxis2 = attr(anchor = "x2", domain = [2a+b, 3a-b], titlefont_size = 14, title=ytitle[2]),
		yaxis3 = attr(anchor = "x3", domain = [2a+b, 3a-b], titlefont_size = 14, title=ytitle[3]),
		xaxis4 = attr(domain = [0a, a-2bx], anchor="y4"),
		xaxis5 = attr(domain = [1a+bx, 2a-bx], anchor="y5"),
		xaxis6 = attr(domain = [2a+2bx, 3a], anchor="y6"),
		yaxis4 = attr(anchor = "x4", domain = [1a+b, 2a-b], titlefont_size = 14, title=ytitle[4]),
		yaxis5 = attr(anchor = "x5", domain = [1a+b, 2a-b], titlefont_size = 14, title=ytitle[5]),
		yaxis6 = attr(anchor = "x6", domain = [1a+b, 2a-b], titlefont_size = 14, title=ytitle[6]),
		xaxis7 = attr(domain = [0a, a-2bx], anchor="y7"),
		xaxis8 = attr(domain = [1a+bx, 2a-bx], anchor="y8"),
		xaxis9 = attr(domain = [2a+2bx, 3a], anchor="y9"),
		yaxis7 = attr(anchor = "x7", domain = [0a+b, 1a-b], titlefont_size = 14, title=ytitle[7]),
		yaxis8 = attr(anchor = "x8", domain = [0a+b, 1a-b], titlefont_size = 14, title=ytitle[8]),
		yaxis9 = attr(anchor = "x9", domain = [0a+b, 1a-b], titlefont_size = 14, title=ytitle[9]),
	)
	plot(data, layout, style=style)
end


full_scats_comp(pv_bench::Vector{T}, pv_noΔ::Vector{T}, pv_nob::Vector{T}, pv_nodef::Vector{T}, tvv::Vector{Vector{Int64}}, key::Symbol, x1::Float64, x2::Float64=x1, x3::Float64=x1, x4::Float64=x1; CI::Bool=false, avg::Bool=false) where T <: AbstractPath = full_scats_comp(pv_bench, pv_noΔ, pv_nob, pv_nodef, tvv, key, x->x/x1, x->x/x2, x->x/x3, x->x/x4, CI=CI, avg=avg)
function full_scats_comp(pv_bench::Vector{T}, pv_noΔ::Vector{T}, pv_nob::Vector{T}, pv_nodef::Vector{T}, tvv::Vector{Vector{Int64}}, key::Symbol, f1::Function=identity, f2::Function=f1, f3::Function=f1, f4::Function=f1; CI::Bool=false, avg::Bool=false, axis::Int64=1) where T <: AbstractPath
	k = 8
	ybench, bench_up, bench_me, bench_lo, bench_av = series_crises(pv_bench, tvv, key, k)
	ynoΔ, noΔ_up, noΔ_me, noΔ_lo, noΔ_av = series_crises(pv_noΔ, tvv, key, k)
	ynob, nob_up, nob_me, nob_lo, nob_av = series_crises(pv_nob, tvv, key, k)
	ynodef, nodef_up, nodef_me, nodef_lo, nodef_av = series_crises(pv_nodef, tvv, key, k)

	colbench = get(ColorSchemes.vik, 0.25)
	colnodef = get(ColorSchemes.vik, 0.75)
	colnoΔ   = get(ColorSchemes.cork, 1)
	colnob   = get(ColorSchemes.vikO, 1)

	line_bench = scatter(x=(-2k:k)/4, y=f1.(bench_me), name="Benchmark",  mode="lines", line_color=colbench, showlegend=(axis==1), legendgroup = 1, xaxis="x$axis", yaxis="y$axis")
	# lbench_avg = scatter(x=(-2k:k)/4, y=f1.(bench_av), name="Benchmark", line_color=col[1])
	# lb_up = scatter(x=(-2k:k)/4, y=f1.(bench_up), hoverinfo="skip", mode="lines", showlegend=false, line_width=0.001, line_color=col[1])
	# lb_lo = scatter(x=(-2k:k)/4, y=f1.(bench_lo), hoverinfo="skip", mode="lines", showlegend=false, line_width=0.001, line_color=col[1], fill="tonexty")
	noΔ2 = scatter(x=-2k/4, y=f2.(noΔ_me)[1], legendgroup=2, name="Δ = 0", mode = "lines+markers", marker_symbol="diamond", showlegend=true)
	line_noΔ = scatter(x=(-2k:k)/4, y=f2.(noΔ_me), name="Δ = 0", line_color=colnoΔ, mode="lines", showlegend=false, legendgroup = 2, xaxis="x$axis", yaxis="y$axis")
	markers_noΔ = scatter(x=((-2k:k)/4)[1:4:end], y=f2.(noΔ_me)[1:4:end], name="Δ = 0", line_color=colnoΔ, mode="markers", marker_symbol = "diamond", showlegend=(axis==1), legendgroup = 2, xaxis="x$axis", yaxis="y$axis")
	# lnoΔ_avg = scatter(x=(-2k:k)/4, y=f2.(noΔ_av), name="Δ = 0", line_color=col[2])
	# lΔ_up = scatter(x=(-2k:k)/4, y=f2.(noΔ_up), hoverinfo="skip", mode="lines", showlegend=false, line_width=0.001, line_color=col[2])
	# lΔ_lo = scatter(x=(-2k:k)/4, y=f2.(noΔ_lo), hoverinfo="skip", mode="lines", showlegend=false, line_width=0.001, line_color=col[2], fill="tonexty")
	line_nob = scatter(x=(-2k:k)/4, y=f3.(nob_me), name="No dom. holdings", line_color=colnob, mode="lines", line_dash="dash", marker_symbol = "circle", showlegend=false, legendgroup = 3, xaxis="x$axis", yaxis="y$axis")
	markers_nob = scatter(x=((-2k:k)/4)[1:4:end], y=f3.(nob_me)[1:4:end], name="No dom. holdings", marker_color=colnob, mode="markers", marker_symbol = "circle", showlegend=(axis==1), legendgroup = 3, xaxis="x$axis", yaxis="y$axis")
	# lnob_avg = scatter(x=(-2k:k)/4, y=f3.(nob_av), name="No dom. holdings", line_color=col[3])
	# lnob_up = scatter(x=(-2k:k)/4, y=f3.(nob_up), hoverinfo="skip", mode="lines", showlegend=false, line_width=0.001, line_color=col[3])
	# lnob_lo = scatter(x=(-2k:k)/4, y=f3.(nob_lo), hoverinfo="skip", mode="lines", showlegend=false, line_width=0.001, line_color=col[3], fill="tonexty")
	line_nodef = scatter(x=(-2k:k)/4, y=f4.(nodef_me), name="No default", line_color=colnodef, mode="lines", showlegend=(axis==1), legendgroup = 4, xaxis="x$axis", yaxis="y$axis", line_dash="dashdot")
	# lnodef_avg = scatter(x=(-k:k)/4, y=f4.(nodef_av), name="No default", line_color=col[4])
	# ln_up = scatter(x=(-k:k)/4, y=f4.(nodef_up), hoverinfo="skip", mode="lines", showlegend=false, line_width=0.001, line_color=col[4])
	# ln_lo = scatter(x=(-k:k)/4, y=f4.(nodef_lo), hoverinfo="skip", mode="lines", showlegend=false, line_width=0.001, line_color=col[4], fill="tonexty")

	s1 = [line_bench, line_noΔ, line_nob, line_nodef, markers_noΔ, markers_nob, noΔ2]
	if avg
		s1 = [lbench_avg, lnoΔ_avg, lnob_avg, lnodef_avg]
	elseif CI
		s1 = [line_bench, lb_up, lb_lo, line_noΔ, lΔ_up, lΔ_lo, line_nob, lnob_up, lnob_lo, line_nodef, ln_up, ln_lo]
	end
	s1
end

function twisted_π(sd::SOEdef, jϵ=floor(Int, N(sd,:ϵ)/2), eval_points::Vector{Int64}=[default_eval_points(sd)...])
	itp_qᵍ = make_itp(sd, sd.eq[:qᵍ], agg=true)
	itp_vf = make_itp(sd, sd.v[:v])
	rep_mat = reshape_long_shocks(sd, sd.gov[:repay])

	jb, jμ, jσ, jξ, jζ, jz = eval_points

	js = find_points(sd, [jb, jμ, jσ, jξ, jζ, jz])
	pϵ = sd.prob[:ϵ][jϵ,:]

	Tπ = zeros(N(sd,:ω))
	Tv = zeros(N(sd,:ω))
	actual_prob = 0.0
	Jgrid = agg_grid(sd)
	for (jω, ωv) in enumerate(sd.gr[:ω])
		ζv = sd.gr[:ζ][jζ]
		jdef = (ζv == 0)

		pz = sd.prob[:z][jz,:]
		pξ = sd.prob[:ξ][jξ,:]

		ϕa = sd.ϕ[:a][jω, jϵ, jb, jμ, jσ, jξ, jζ, jz]
		ϕb = sd.ϕ[:b][jω, jϵ, jb, jμ, jσ, jξ, jζ, jz]

		for (jξp, ξpv) in enumerate(sd.gr[:ξ]), (jζp, ζpv) in enumerate(sd.gr[:ζ]), (jzp, zpv) in enumerate(sd.gr[:z])
			bpv = sd.eq[:issuance][js]
			μpv = sd.LoM[:μ][js, jξp, jzp][jζp]
			σpv = sd.LoM[:σ][js, jξp, jzp][jζp]

			prob = pz[jzp] * pξ[jξp]
			jdefp = (ζpv == 0)
			if jdefp && jdef # Default in both
				prob *= 1-sd.pars[:θ]
				rep = itp_qᵍ(bpv, μpv, σpv, ξpv, ζpv, zpv)
			elseif jdef # Reentering markets
				prob *= sd.pars[:θ]
				rep = sd.pars[:κ] + (1-sd.pars[:ρ]) * itp_qᵍ(bpv, μpv, σpv, ξpv, ζpv, zpv)
			elseif jdefp # Default at t+1
				bpv *= (1-sd.pars[:ℏ])
				prob *= 1-rep_mat[jb, jμ, jσ, jξ, jζ, jz, jξp, jzp]
				rep = (1-sd.pars[:ℏ]) * itp_qᵍ(bpv, μpv, σpv, ξpv, ζpv, zpv)
			else # repayment
				prob *= rep_mat[jb, jμ, jσ, jξ, jζ, jz, jξp, jzp]
				rep = sd.pars[:κ] + (1-sd.pars[:ρ]) * itp_qᵍ(bpv, μpv, σpv, ξpv, ζpv, zpv)
			end

			ωpv = ϕa + rep * ϕb

			rep_prob = rep_mat[jb, jμ, jσ, jξ, jζ, jz, jξp, jzp]

			for (jϵp, ϵpv) in enumerate(sd.gr[:ϵ])
				Tπ[jω] += prob * pϵ[jϵp] * itp_vf(ωpv, ϵpv, bpv, μpv, σpv, ξpv, ζpv, zpv)^(1-sd.pars[:γ]) * jdefp
				Tv[jω] += prob * pϵ[jϵp] * itp_vf(ωpv, ϵpv, bpv, μpv, σpv, ξpv, ζpv, zpv)^(1-sd.pars[:γ])
				actual_prob += prob * pϵ[jϵp] * jdefp * (jω==1)
			end
		end
	end

	Tπ = Tπ ./ Tv

	return Tπ, actual_prob
end

function make_twisted(sd::SOEdef, eval_points::Vector{Int64}=[default_eval_points(sd)...]; style::Style=slides_def, custom_points::Bool=true, leg::Bool=false)
	Tπ = zeros(N(sd,:ω), N(sd,:ϵ))

	if custom_points
		eval_points = [13,4,2,2,2,7]
	end

	if leg
		legattr = attr(orientation="v", x=0.75, y=0.0)
	else
		legattr = attr()
	end

	for jϵ in 1:N(sd,:ϵ)
		Tπ[:,jϵ], _ = twisted_π(sd, jϵ, eval_points)
	end
	_, actual_prob = twisted_π(sd, 1, eval_points)

	plot([
		scatter(x=extrema(sd.gr[:ω]), y=ones(2) * 100*actual_prob, line_dash="dash", mode="lines", line_color=col[3], name="Actual")
		[scatter(x=sd.gr[:ω], y=100*Tπ[:,jϵ], name="<i>ϵ = $(@sprintf("%0.3g", ϵv))", showlegend=leg, line_color=get(ColorSchemes.lajolla, 1-(jϵ-1)/N(sd,:ϵ))) for (jϵ, ϵv) in enumerate(sd.gr[:ϵ])]
		], style=style, Layout(xaxis_title="<i>ω", yaxis_title="<i>%", title="Twisted default probabilities", legend=legattr))
end

function Wr_Wd(sd::SOEdef, type::String)
	B′ = sd.eq[:issuance]	

	if type == "Wr-Wd"
		itp_obj = make_itp(sd, sd.eq[:welfare]; agg=true)
	elseif type == "T"
		itp_obj = make_itp(sd, -sd.eq[:T]; agg=true)
	end

	# More μ means default more often
	μ_gov = 0.001 * 0.0
	σ_gov = 0.004

	Wr = zeros(size(sd.eq[:wage])..., N(sd,:z))
	Wd = zeros(size(sd.eq[:wage])..., N(sd,:z))

	Jgrid = agg_grid(sd);
	rep_prob = zeros(N(sd,:b), N(sd,:μ), N(sd,:σ), N(sd,:ξ), N(sd,:ζ), N(sd,:z), N(sd,:ξ), N(sd,:z))
	for js in 1:size(Jgrid, 1)
		μ′_arr = sd.LoM[:μ][js,:,:]
		σ′_arr = sd.LoM[:σ][js,:,:]

		jb, jμ, jσ, jξ, jζ, jz = Jgrid[js, :]

		ζv = sd.gr[:ζ][Jgrid[js, 5]]
		
		bpv = B′[js]
		for (jzp, zpv) in enumerate(sd.gr[:z])
			for (jξp, ξpv) in enumerate(sd.gr[:ξ])
				prob = sd.prob[:ξ][jξ,jξp]

				jζp = 1 # Default at t+1
				ζpv = sd.gr[:ζ][jζp]
				μpv = μ′_arr[jξp, jzp][jζp]
				σpv = σ′_arr[jξp, jzp][jζp]
				Wd[js, jzp] += prob * itp_obj((1.0-sd.pars[:ℏ])*bpv, μpv, σpv, ξpv, ζpv, zpv)

				jζp = 2 # No default at t+1
				ζpv = sd.gr[:ζ][jζp]
				μpv = μ′_arr[jξp, jzp][jζp]
				σpv = σ′_arr[jξp, jzp][jζp]
				Wr[js, jzp] += prob * itp_obj(bpv, μpv, σpv, ξpv, ζpv, zpv)
			end
		end
	end

	return Wr, Wd
end

function earnings_default(sd::SOEdef)
	rep_mat = reshape_long_shocks(sd, sd.gov[:repay])
	itp_wL = make_itp(sd, sd.eq[:wage] .* sd.eq[:Ld], agg=true)
	itp_qᵍ = make_itp(sd, sd.eq[:qᵍ], agg=true)

	Ey = zeros(size(sd.eq[:wage])..., N(sd,:z))
	Eret = zeros(size(sd.eq[:wage])..., N(sd,:z))
	Eyd = zeros(size(sd.eq[:wage])..., N(sd,:z))
	Jgrid = agg_grid(sd)

	for js in 1:size(Jgrid,1)
		jξ = Jgrid[js, 4]
		jζ = Jgrid[js, 5]
		jz = Jgrid[js, 6]

		ζv = sd.gr[:ζ][jζ]
		jdef = (ζv == 0)

		pz = sd.prob[:z][jz,:]
		pξ = sd.prob[:ξ][jξ,:]

		sumprob = 0.0
		for (jξp, ξpv) in enumerate(sd.gr[:ξ]), (jzp, zpv) in enumerate(sd.gr[:z]), (jζp, ζpv) in enumerate(sd.gr[:ζ])
			prob = pξ[jξp]
			jdefp = (ζpv == 0)
			μpv = sd.LoM[:μ][js, jξp, jzp][jζp]
			σpv = sd.LoM[:σ][js, jξp, jzp][jζp]

			bpv = sd.eq[:issuance][js]
			
			if jdefp && jdef # Default in both
				prob *= 1-sd.pars[:θ]
				rep = itp_qᵍ(bpv, μpv, σpv, ξpv, ζpv, zpv)
			elseif jdef # Reentering markets
				prob *= sd.pars[:θ]
				rep = sd.pars[:κ] + (1-sd.pars[:ρ]) * itp_qᵍ(bpv, μpv, σpv, ξpv, ζpv, zpv)
			elseif jdefp # Default at t+1
				bpv *= (1-sd.pars[:ℏ])
				prob *= 1-rep_mat[Jgrid[js,:]..., jξp, jzp]
				rep = (1-sd.pars[:ℏ]) * itp_qᵍ(bpv, μpv, σpv, ξpv, ζpv, zpv)
			else # repayment
				prob *= rep_mat[Jgrid[js,:]..., jξp, jzp]
				rep = sd.pars[:κ] + (1-sd.pars[:ρ]) * itp_qᵍ(bpv, μpv, σpv, ξpv, ζpv, zpv)
			end

			# Ey[js, jzp] += prob * itp_wL(bpv, μpv, σpv, ξpv, ζpv, zpv)
			Eret[js, jzp] += prob * rep / sd.eq[:qᵍ][js]
			Ey[js, jzp] += prob * itp_wL(bpv, μpv, σpv, ξpv, 1, zpv)
			Eyd[js, jzp] += prob * itp_wL(bpv, μpv, σpv, ξpv, 0, zpv)
		end
		# !isapprox(sumprob, 1) && println("WARNING: SUM OF PROBS = $sumprob")
	end

	return Ey, Eret, Eyd
end

function find_points(sd::SOEdef, jlong)
	Jgrid = agg_grid(sd)

	js = findfirst([Jgrid[jjs, :] == jlong for jjs in 1:size(Jgrid,1)])
end

function get_earnings_default(sd::SOEdef, eval_points::Vector{Int64}=[default_eval_points(sd)...])
	Ey, Eret, Eyd = earnings_default(sd)

	jb, jμ, jσ, jξ, jζ, jz = eval_points

	js = find_points(sd, [jb, jμ, jσ, jξ, jζ, jz])

	return Ey[js,:], Eret[js,:], Eyd[js,:]
end

function get_Wr_Wd(sd::SOEdef, eval_points::Vector{Int64}=[default_eval_points(sd)...]; type::String="Wr-Wd")
	Wr, Wd = Wr_Wd(sd, type)

	jb, jμ, jσ, jξ, jζ, jz = eval_points

	js = find_points(sd, [jb, jμ, jσ, jξ, jζ, jz])

	return Wr[js,:], Wd[js,:]
end

function make_panels(sd::SOEdef, type::String; style::Style=slides_def, leg::Bool=true)
	Jgrid = agg_grid(sd)

	jbv = [floor(Int, N(sd,:b)*0.25), 10, N(sd,:b)]
	_, jμ, jσ, jξ, jζ, jz = default_eval_points(sd)

	Eyv = Vector{Vector{Float64}}(undef, 0)
	Eydv = Vector{Vector{Float64}}(undef, 0)
	Erv = Vector{Vector{Float64}}(undef, 0)
	Ey2 = Vector{Vector{Float64}}(undef, 0)
	Eyd2 = Vector{Vector{Float64}}(undef, 0)
	Er2 = Vector{Vector{Float64}}(undef, 0)
	for (jj,jb) in enumerate(jbv)
		if type == "Earnings_Default"
			Ey, Er, Eyd = get_earnings_default(sd, [jb,jμ, jσ, 1, jζ, jz])
			Eyy, Err, Eyyd = get_earnings_default(sd, [jb,jμ, jσ, 2, jζ, jz])
			cor = "(y,ret)"
		elseif type == "Wr-Wd"
			Ey, Er = get_Wr_Wd(sd, [jb,jμ, jσ, 1, jζ, jz])
			Eyy, Err = get_Wr_Wd(sd, [jb,jμ, jσ, 2, jζ, jz])
			Eyd = Ey.*NaN
			Eyyd = Ey.*NaN
			cor = "(Wr, Wd)"
		elseif type == "T"
			Ey, Er = get_Wr_Wd(sd, [jb,jμ, jσ, 1, jζ, jz], type="T")
			Eyy, Err = get_Wr_Wd(sd, [jb,jμ, jσ, 2, jζ, jz], type="T")
			Eyd = Ey.*NaN
			Eyyd = Ey.*NaN
			cor = "(Tr, Td)"
		end

		co = [cov(Ey, Er) / (std(Ey)*std(Er)) for js in 1:size(Jgrid,1)]
		println("Corr"*cor*" = $(mean(co))")
		co = [cov(Eyy, Err) / (std(Eyy)*std(Err)) for js in 1:size(Jgrid,1)]
		println("Corr"*cor*" = $(mean(co))")

		push!(Eyv, Ey)
		push!(Eydv, Eyd)
		push!(Erv, Er)
		push!(Ey2, Eyy)
		push!(Eyd2, Eyyd)
		push!(Er2, Err)
	end

	col1 = get(ColorSchemes.romaO, 0.8)
	col2 = get(ColorSchemes.vikO, 0.75)
	col3 = get(ColorSchemes.romaO, 0.3)

	if type == "Earnings_Default"
		title = "Income and returns on government debt"
		ytitle1 = "𝔼[<i>y<sup>r</sup></i>]"
		ytitle1d = "𝔼[<i>y<sup>d</sup></i>]"
		ytitle2 = "𝔼[<i>R<sup>b</sup></i>]"
		modelines = "lines"
	elseif type == "Wr-Wd"
		title = "Value of repayment and default"
		ytitle1 = "<i>W<sup>r"
		ytitle1d = ""
		ytitle2 = "<i>W<sup>d"
		col2, col3 = col3, col2
		modelines = "lines+markers"
	elseif type == "T"
		title = "Transfers"
		ytitle1 = "<i>T<sup>r"
		ytitle1d = ""
		ytitle2 = "<i>T<sup>d"
		col2, col3 = col3, col2
		modelines = "lines+markers"
	end


	minyy, maxyy = extrema(vcat(vcat(Eyv[:]...), vcat(Ey2[:]...)))
	minyy = minyy - 0.05 * (maxyy-minyy)
	maxyy = maxyy + 0.05 * (maxyy-minyy)

	minyr, maxyr = extrema(vcat(vcat(Erv[:]...), vcat(Er2[:]...)))
	minyr = minyr - 0.05 * (maxyr-minyr)
	maxyr = maxyr + 0.05 * (maxyr-minyr)

	maxy = max(maxyy, maxyr)
	miny = min(minyy, minyr)

	if type == "Earnings_Default"
		minyd, maxyd = extrema(vcat(vcat(Eydv[:]...), vcat(Eyd2[:]...)))
		minyd = minyd - 0.05 * (maxyd-minyd)
		maxyd = maxyd + 0.05 * (maxyd-minyd)

		maxy = max(maxy, maxyd)
		miny = min(miny, minyd)
	end

	ξv = sd.gr[:ξ]

	data = [
		[scatter(x=sd.gr[:z], y=Eyv[jb], xaxis="x$jb",yaxis="y$jb", name=ytitle1, marker_color=col1, showlegend=((jb==1)&(leg)), line_width=2.5) for jb in 1:length(Eyv)]
		[scatter(x=sd.gr[:z], y=Erv[jb], xaxis="x$jb",yaxis="y$jb", name=ytitle2, marker_color=col3, mode=modelines, line_dash="dot", showlegend=((jb==1)&(leg)), line_width=2.5) for jb in 1:length(Eyv)]
		[scatter(x=sd.gr[:z], y=Ey2[jb], xaxis="x$(jb+3)",yaxis="y$(jb+3)", name=ytitle1, marker_color=col1, showlegend=false, line_width=2.5) for jb in 1:length(Ey2)]
		[scatter(x=sd.gr[:z], y=Er2[jb], xaxis="x$(jb+3)",yaxis="y$(jb+3)", name=ytitle2, marker_color=col3, mode=modelines, line_dash="dot", showlegend=false, line_width=2.5) for jb in 1:length(Ey2)]
	]

	if type == "Earnings_Default"
		for jb in 1:length(Eyv)
			push!(data, scatter(x=sd.gr[:z], y=Eydv[jb], xaxis="x$jb",yaxis="y$jb", name=ytitle1d, marker_color=col2, line_dash="dashdot", showlegend=((jb==1)&(leg)), line_width=2.5) )
		end
		for jb in 1:length(Ey2)
			push!(data, scatter(x=sd.gr[:z], y=Eyd2[jb], xaxis="x$(jb+3)",yaxis="y$(jb+3)", name=ytitle1d, marker_color=col2, line_dash="dashdot", showlegend=false, line_width=2.5))
		end
	end

	annotations = [
		[attr(x=0, xanchor="center", xref="x$jb", y=1.01, yref="paper", yanchor="bottom", text="<i>B = $(@sprintf("%0.3g", sd.gr[:b][bv])), ξ = $(@sprintf("%0.3g",100*ξv[1]))%", showarrow=false) for (jb,bv) in enumerate(jbv)]
		[attr(x=0, xanchor="center", xref="x$jb", y=0.44, yref="paper", yanchor="bottom", text="<i>B = $(@sprintf("%0.3g", sd.gr[:b][bv])), ξ = $(@sprintf("%0.3g",100*ξv[2]))%", showarrow=false) for (jb,bv) in enumerate(jbv)]
		]

	layout = Layout(
		yaxis1=attr(anchor="x1", domain = [0.575,1]),
		yaxis2=attr(anchor="x2", domain = [0.575,1]),
		yaxis3=attr(anchor="x3", domain = [0.575,1]),
		yaxis4=attr(anchor="x4", domain = [0,0.425]),
		yaxis5=attr(anchor="x5", domain = [0,0.425]),
		yaxis6=attr(anchor="x6", domain = [0,0.425]),
		xaxis1=attr(zeroline=false, domain=[0,0.3], anchor="y1"),
		xaxis2=attr(zeroline=false, domain=[0.33, 0.67], anchor="y2"),
		xaxis3=attr(zeroline=false, domain=[0.7, 1], anchor="y3"),
		xaxis4=attr(title="<i>z′", zeroline=false, domain=[0,0.3], anchor="y4"),
		xaxis5=attr(title="<i>z′", zeroline=false, domain=[0.33, 0.67], anchor="y5"),
		xaxis6=attr(title="<i>z′", zeroline=false, domain=[0.7, 1], anchor="y6"),
		title=ifelse(leg,title,""), height = 1080*0.45,
		legend = attr(orientation="v", x=0.95, xanchor="right", y=0.95),
		annotations = annotations
		)
	plot(data, layout, style=style)
end




























function lines(h::Hank, y, x_dim, name=""; custom_w::Int=0)
	jshow_b, jshow_μ, jshow_σ, jshow_ξ, jshow_ζ, jshow_z = ceil(Int, h.Nb/2), ceil(Int, h.Nμ/2), ceil(Int, h.Nσ/2), 2, 1, ceil(Int, h.Nz/2)

	if custom_w != 0
		jshow_ξ = custom_w
	end

	x = h.bgrid
	xlabel = "B"
	if x_dim == 1
		y = y[:, jshow_μ, jshow_σ, jshow_ξ, jshow_ζ, jshow_z]
	elseif x_dim == 2
		x, xlabel = h.μgrid, "μ"
		y = y[jshow_b, :, jshow_σ, jshow_ξ, jshow_ζ, jshow_z]
	elseif x_dim == 3
		x, xlabel = h.σgrid, "σ"
		y = y[jshow_b, jshow_μ, :, jshow_ξ, jshow_ζ, jshow_z]
	elseif x_dim == 4
		x, xlabel = h.ξgrid, "w"
		y = y[jshow_b, jshow_μ, jshow_σ, :, jshow_ζ, jshow_z]
	elseif x_dim == 6
		x, xlabel = h.zgrid, "z"
		y = y[jshow_b, jshow_μ, jshow_σ, jshow_ξ, jshow_ζ, :]
	else
		print_save("x_dim wrong")
	end


	layout = Layout(;	xaxis=attr(title=xlabel, zeroline=true),
						yaxis=attr(zeroline=true),
						font_family="Fira Sans Light")

	l = scatter(;x=x, y=y, showlegend=false)
	p = plot(l, layout)
	if name == ""
	else
		relayout!(p, title=name)
	end
	return p
end

function plot_gov_welf(h::Hank; remote::Bool=false)
	itp_vf = make_itp(h, h.vf; agg=false)

	B′_mat = reshape(h.issuance, h.Nb, h.Nμ, h.Nσ, h.Nξ, h.Nζ, h.Nz)
	μ′_mat = reshape(h.μ′, h.Nb, h.Nμ, h.Nσ, h.Nξ, h.Nζ, h.Nz, h.Nz, 2)
	σ′_mat = reshape(h.σ′, h.Nb, h.Nμ, h.Nσ, h.Nξ, h.Nζ, h.Nz, h.Nz, 2)

	Wr_vec = zeros(size(h.Jgrid, 1))
	Wd_vec = zeros(size(h.Jgrid, 1))
	for js in 1:length(Wr_vec)
		jb = h.Jgrid[js, 1]
		jμ = h.Jgrid[js, 2]
		jσ = h.Jgrid[js, 3]
		jξ = h.Jgrid[js, 4]
		jζ = 1
		jz = h.Jgrid[js, 6]

		EWr, EWd = 0., 0.

		bvp = B′_mat[jb, jμ, jσ, jξ, jζ, jz]
		for (jξp, ϵpv) in enumerate(h.ξgrid), jzp in 1:h.Nz
			prob = h.Pξ[jξ, jξp] * h.Pz[jz, jzp]
			μvp = μ′_mat[jb, jμ, jσ, jξ, jζ, jz, jzp, 1]
			σvp = σ′_mat[jb, jμ, jσ, jξ, jζ, jz, jzp, 1]
			EWr += prob * integrate_itp(h, bvp, μvp, σvp, ξvp, 1, jzp, itp_vf)
			μvp = μ′_mat[jb, jμ, jσ, jξ, jζ, jz, jzp, 2]
			σvp = σ′_mat[jb, jμ, jσ, jξ, jζ, jz, jzp, 2]
			EWd += prob * integrate_itp(h, (1.0 .-h.ℏ)*bvp, μvp, σvp, ξvp, 2, jzp, itp_vf)
		end

		Wr_vec[js] = EWr
		Wd_vec[js] = EWd
	end

	Wr_mat = reshape(Wr_vec, h.Nb, h.Nμ, h.Nσ, h.Nξ, h.Nζ, h.Nz)
	Wd_mat = reshape(Wd_vec, h.Nb, h.Nμ, h.Nσ, h.Nξ, h.Nζ, h.Nz)
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
		path = pwd() * "/../Graphs/"
		# savefig(p, path * "objfunc.pdf")
		return p
	end
	Void
end

function plot_govt_reaction(h::Hank; Wdiff::Bool=false, Ts::Bool=false, Bs::Bool=false, Ws::Bool=false, qs::Bool=false, remote::Bool=false)
	jμ, jσ, jξ = ceil(Int, h.Nμ*0.75), ceil(Int, h.Nσ*0.25), 2
	# jμ, jσ, jξ = ceil(Int, h.Nμ*0.25), ceil(Int, h.Nσ*0.75), 2
	μv, σv, ξv = h.μgrid[jμ], h.σgrid[jσ], h.ξgrid[jξ]
	jζ = 1

	itp_vf = make_itp(h, h.vf; agg=false)

	wL = h.Ld .* h.wage .* (1.0-h.τ)
	T_mat = govt_bc(h, h.wage .* h.Ld)
	itp_T = make_itp(h, -T_mat; agg=true)
	wL_mat = reshape(wL, h.Nb, h.Nμ, h.Nσ, h.Nξ, h.Nζ, h.Nz)


	B′_mat = reshape(h.issuance, h.Nb, h.Nμ, h.Nσ, h.Nξ, h.Nζ, h.Nz)
	μ′_mat = reshape(h.μ′, h.Nb, h.Nμ, h.Nσ, h.Nξ, h.Nζ, h.Nz, h.Nξ, h.Nz, 2)
	σ′_mat = reshape(h.σ′, h.Nb, h.Nμ, h.Nσ, h.Nξ, h.Nζ, h.Nz, h.Nξ, h.Nz, 2)

	rep_mat = reshape(h.repay, h.Nb, h.Nμ, h.Nσ, h.Nξ, h.Nζ, h.Nz, h.Nξ, h.Nz)

	itp_B = make_itp(h, B′_mat; agg=true)
	itp_W = make_itp(h, wL_mat; agg=true)
	itp_qᵍ = make_itp(h, h.qᵍ; agg=true)

	midb = ceil(Int, h.Nb/2)
	midb = h.Nb-3
	states = gridmake([3; midb; h.Nb], [1; h.Nξ])
	jz = ceil(Int, h.Nz/2)
	jξ = 1
	p_vec = Array{PlotlyJS.SyncPlot}(undef, size(states,1))
	for js in 1:size(states,1)
		Wr = zeros(h.Nz)
		Wd = zeros(h.Nz)
		Tr, Td = zeros(h.Nz), zeros(h.Nz)
		Br, Bd = zeros(h.Nz), zeros(h.Nz)
		Yr, Yd = zeros(h.Nz), zeros(h.Nz)
		qr, qd = zeros(h.Nz), zeros(h.Nz)
		exp_rep = zeros(h.Nz)
		jb, jξp = states[js, :]
		bvp = B′_mat[jb, jμ, jσ, jξ, jζ, jz]
		ξpv = h.ξgrid[jξp]
		for jzp in 1:h.Nz
			μvp = μ′_mat[jb, jμ, jσ, jξ, jζ, jz, jξp, jzp, 1]
			σvp = σ′_mat[jb, jμ, jσ, jξ, jζ, jz, jξp, jzp, 1]
			exp_rep[jzp] = (rep_mat[jb, jμ, jσ, jξ, jζ, jz, jξp, jzp] * (h.κ + (1.0-h.ρ) * itp_qᵍ[bvp, μvp, σvp, ξpv, 1, jzp] ) + (1.0-rep_mat[jb, jμ, jσ, jξ, jζ, jz, jξp, jzp]) * (1.0-h.ρ) * (1.0-h.ℏ) * itp_qᵍ[(1.0 .-h.ℏ)*bvp, μvp, σvp, ξpv, 2, jzp])	/ h.qᵍ[js]
			# exp_rep[jzp] = rep_mat[jb, jμ, jσ, jξ, jζ, jz, jξp, jzp]
			qr[jzp] = itp_qᵍ[bvp, μvp, σvp, ξpv, 1, jzp]
			Wr[jzp] = integrate_itp(h, bvp, μvp, σvp, ξpv, 1, jzp, itp_vf) + 0.0075
			Tr[jzp] = itp_T[bvp, μvp, σvp, ξpv, 1, jzp]
			Br[jzp] = itp_B[bvp, μvp, σvp, ξpv, 1, jzp]
			Yr[jzp] = itp_W[bvp, μvp, σvp, ξpv, 1, jzp]
			μvp = μ′_mat[jb, jμ, jσ, jξ, jζ, jz, jξp, jzp, 2]
			σvp = σ′_mat[jb, jμ, jσ, jξ, jζ, jz, jξp, jzp, 2]
			qd[jzp] = itp_qᵍ[(1.0 .-h.ℏ)*bvp, μvp, σvp, ξpv, 2, jzp]
			Wd[jzp] = integrate_itp(h, (1.0 .-h.ℏ)*bvp, μvp, σvp, ξpv, 2, jzp, itp_vf)
			Td[jzp] = itp_T[(1.0 .-h.ℏ)*bvp, μvp, σvp, ξpv, 2, jzp]
			Bd[jzp] = itp_B[(1.0 .-h.ℏ)*bvp, μvp, σvp, ξpv, 2, jzp]
			Yd[jzp] = itp_W[(1.0 .-h.ℏ)*bvp, μvp, σvp, ξpv, 2, jzp]
		end
		if Wdiff 
			p_vec[js] = plot(scatter(;x=h.zgrid, y=Wd-Wr, marker_color=col[1], showlegend=false, line_width=2), Layout(;title="𝐵=$(h.bgrid[jb]), ξ'=$(round((h.ξgrid[jξp]),digits=2))", titlefont_size=32))
		elseif Ts
			p_vec[js] = plot([
				scatter(;x=h.zgrid, y=Tr, marker_color=col[1], showlegend=false, line_width = 2),
				scatter(;x=h.zgrid, y=Td, marker_color=col[4], showlegend=false, line_dash="dashdot", line_width = 2),
				], Layout(;title="𝐵=$(h.bgrid[jb]), ξ'=$(round((h.ξgrid[jξp]),digits=2))", titlefont_size=32)
				)
		elseif Bs
			p_vec[js] = plot([
				scatter(;x=h.zgrid, y=Br, marker_color=col[1], showlegend=false, line_width = 2),
				scatter(;x=h.zgrid, y=Bd, marker_color=col[4], showlegend=false, line_dash="dashdot", line_width = 2),
				], Layout(;title="𝐵=$(h.bgrid[jb]), ξ'=$(round((h.ξgrid[jξp]),digits=2))", titlefont_size=32)
				)
		elseif qs
			p_vec[js] = plot([
				scatter(;x=h.zgrid, y=qr, marker_color=col[1], showlegend=false, line_width = 2),
				scatter(;x=h.zgrid, y=qd, marker_color=col[4], showlegend=false, line_dash="dashdot", line_width = 2),
				], Layout(;title="𝐵=$(h.bgrid[jb]), ξ'=$(round((h.ξgrid[jξp]),digits=2))", titlefont_size=32)
				)
		elseif Ws
			p_vec[js] = plot([
				scatter(;x=h.zgrid, y=Yr, marker_color=col[1], showlegend=false, line_width = 2),
				scatter(;x=h.zgrid, y=exp_rep, marker_color=col[2], showlegend=false, line_dash="dot", line_width = 2, mode ="lines"),
				scatter(;x=h.zgrid, y=Yd, marker_color=col[4], showlegend=false, line_dash="dashdot", line_width = 2),
				], Layout(;title="𝐵=$(h.bgrid[jb]), ξ'=$(round((h.ξgrid[jξp]),digits=2))", titlefont_size=32)
				)
			# , cov = $(round(cov(Yr, exp_rep)/std(Yr)/std(exp_rep),3))
		else

			p_vec[js] = plot(  [scatter(;x=h.zgrid, y=Wr, marker_color=col[1], showlegend=false, line_width=2),
						scatter(;x=h.zgrid, y=Wd, marker_color=col[4], showlegend=false, line_dash="dashdot", line_width=2)],
						Layout(;title="𝐵=$(h.bgrid[jb]), ξ'=$(round((h.ξgrid[jξp]),digits=2))", titlefont_size=32))
		end
	end

	p_paper = [p_vec[1] p_vec[2] p_vec[3]; p_vec[4] p_vec[5] p_vec[6]]
	if Ws
		p_paper = [p_vec[1] p_vec[2] p_vec[3]]
	end
	p_paper.plot.layout["font_family"] = "STIX Two Text"
	p_slides = [p_vec[1] p_vec[3]; p_vec[4] p_vec[6]]
	p_slides = [p_vec[1] p_vec[2] p_vec[3]; p_vec[4] p_vec[5] p_vec[6]]
	p_slides.plot.layout["font_family"] = "Fira Sans Light"
	p_slides.plot.layout["plot_bgcolor"] = "rgba(250, 250, 250, 1.0)"
	p_slides.plot.layout["paper_bgcolor"] = "rgba(250, 250, 250, 1.0)"

	if remote
		path = pwd() * "/../../Graphs/"
		save(path * "p_reactions.jld", "p", p)
	else
		path = pwd() * "/../Graphs/"
		if is_linux()
			return p_paper, p_slides
		else
			# savefig(p, path * "reactions.pdf")
		end
	end
	Void
end

function plot_debtprice(h::Hank; remote::Bool=false)

	qʰ_mat, qᵍ_mat, wL_mat, T_mat, pC_mat, Π_mat = _unpackstatefs(h)
	T_vec  = reshape(T_mat, length(T_mat))
	Π_vec  = reshape(Π_mat, length(Π_mat))
	wL_vec = reshape(wL_mat, length(wL_mat))

	ϕc_mat = h.ϕc
	yd_mat = zeros(size(h.ϕc))
	pC_big = zeros(size(h.ϕc))

	adj = sum(h.λϵ.*exp.(h.ϵgrid))
	agg_income = wL_vec + Π_vec / adj

	def_prob = zeros(h.Nb, h.Nμ, h.Nσ, h.Nξ, h.Nζ, h.Nz)
	rep_mat = reshape(h.repay, h.Nb, h.Nμ, h.Nσ, h.Nξ, h.Nζ, h.Nz, h.Nz)

	jpn = ceil(Int, length(h.pngrid)/2)
	pnv = h.pngrid[jpn]
	N = size(h.Jgrid, 1)
	wage_pn, labor_pn, profits_pn = Array{Float64, 1}(undef, N), Array{Float64, 1}(undef, N), Array{Float64, 1}(undef, N)
	for js in 1:N
		jζ = h.Jgrid[js, 5]
		jz = h.Jgrid[js, 6]

		wv = h.wbar
		ζv = h.ζgrid[jζ]
		zv = h.zgrid[jz]

		labor_pn[js], wage_pn[js], profits_pn[js], _ = labor_market(h, ζv, zv, pnv)
	end

	pC = price_index(h, pnv)
	pC_fix = ones(h.Nb, h.Nμ, h.Nσ, h.Nξ, h.Nζ, h.Nz) * pC

	T_fix = govt_bc(h, wage_pn.*labor_pn)# - reshape(profits_pn - h.profits, h.Nb, h.Nμ, h.Nσ, h.Nξ, h.Nζ, h.Nz)
	Π_fix = reshape(profits_pn, h.Nb, h.Nμ, h.Nσ, h.Nξ, h.Nζ, h.Nz)

	wL_fix  = reshape(wage_pn.*labor_pn, h.Nb, h.Nμ, h.Nσ, h.Nξ, h.Nζ, h.Nz) * (1.0 - h.τ)
	yd_fix = zeros(size(h.ϕc))
	pC_bigfix = zeros(size(h.ϕc))
	for js in 1:size(h.Jgrid, 1)
		jb = h.Jgrid[js, 1]
		jμ = h.Jgrid[js, 2]
		jσ = h.Jgrid[js, 3]
		jξ = h.Jgrid[js, 4]
		jζ = h.Jgrid[js, 5]
		jz = h.Jgrid[js, 6]
		for (jϵ, ϵv) in enumerate(h.ϵgrid), (jω, ωv) in enumerate(h.ωgrid)
			yd_mat[jω, jϵ, jb, jμ, jσ, jξ, jζ, jz] = ωv + agg_income[js] * exp(ϵv) - T_vec[js]
			pC_big[jω, jϵ, jb, jμ, jσ, jξ, jζ, jz] = pC_mat[js]

			yd_fix[jω, jϵ, jb, jμ, jσ, jξ, jζ, jz] = ωv + (wL_fix[jb, jμ, jσ, jξ, jζ, jz] + Π_fix[jb, jμ, jσ, jξ, jζ, jz]/adj) * exp(ϵv) - T_fix[jb, jμ, jσ, jξ, jζ, jz]
			pC_bigfix[jω, jϵ, jb, jμ, jσ, jξ, jζ, jz] = pC_fix[jb, jμ, jσ, jξ, jζ, jz]
		end
		for jzp in 1:h.Nz
			def_prob[jb, jμ, jσ, jξ, jζ, jz] += h.Pz[jz, jzp] * (1.0 .-rep_mat[jb, jμ, jσ, jξ, jζ, jz, jzp])
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

	jω1, jω2 = 2, 4
	jϵ_show = ceil(Int, h.Nϵ/2)
	pc1p = lines(h, Srate[jω1, jϵ_show,:,:,:,:,:,:],  1, "Saving rate at ω = $(round(h.ωgrid[jω1],digits=2))")
	pc2p = lines(h, Srate[jω1, jϵ_show,:,:,:,:,:,:],  2)
	pc3p = lines(h, Srate[jω1, jϵ_show,:,:,:,:,:,:],  3)
	pc4p = lines(h, Srate[jω1, jϵ_show,:,:,:,:,:,:],  4)
	pc6p = lines(h, Srate[jω1, jϵ_show,:,:,:,:,:,:],  6)

	pc1r = lines(h, Srate[jω2, jϵ_show,:,:,:,:,:,:],  1, "Saving rate at ω = $(round(h.ωgrid[jω2],digits=2))")
	pc2r = lines(h, Srate[jω2, jϵ_show,:,:,:,:,:,:],  2)
	pc3r = lines(h, Srate[jω2, jϵ_show,:,:,:,:,:,:],  3)
	pc4r = lines(h, Srate[jω2, jϵ_show,:,:,:,:,:,:],  4)
	pc6r = lines(h, Srate[jω2, jϵ_show,:,:,:,:,:,:],  6)
	
	pc1pf = lines(h, Sratef[jω1, jϵ_show,:,:,:,:,:,:],  1, "S/Y at ω = $(round(h.ωgrid[jω1],digits=2)), fixed pN")
	pc2pf = lines(h, Sratef[jω1, jϵ_show,:,:,:,:,:,:],  2)
	pc3pf = lines(h, Sratef[jω1, jϵ_show,:,:,:,:,:,:],  3)
	pc4pf = lines(h, Sratef[jω1, jϵ_show,:,:,:,:,:,:],  4)
	pc6pf = lines(h, Sratef[jω1, jϵ_show,:,:,:,:,:,:],  6)

	pc1rf = lines(h, Sratef[jω2, jϵ_show,:,:,:,:,:,:],  1, "S/Y at ω = $(round(h.ωgrid[jω2],digits=2)), fixed pN")
	pc2rf = lines(h, Sratef[jω2, jϵ_show,:,:,:,:,:,:],  2)
	pc3rf = lines(h, Sratef[jω2, jϵ_show,:,:,:,:,:,:],  3)
	pc4rf = lines(h, Sratef[jω2, jϵ_show,:,:,:,:,:,:],  4)
	pc6rf = lines(h, Sratef[jω2, jϵ_show,:,:,:,:,:,:],  6)


	p = [pq1 pq2 pq3 pq4 pq6; pd1 pd2 pd3 pd4 pd6; pc1p pc2p pc3p pc4p pc6p; pc1r pc2r pc3r pc4r pc6r; pc1pf pc2pf pc3pf pc4pf pc6pf; pc1rf pc2rf pc3rf pc4rf pc6rf]
	p.plot.layout["width"] = 800
	p.plot.layout["height"] = 800*1.15
	p.plot.layout["font_family"] = "Fira Sans Light"
	if remote
		path = pwd() * "/../../Graphs/"
		save(path * "p_debtprice.jld", "p", p)
	else
		path = pwd() * "/../Graphs/"
		# savefig(p, path * "debtprice.pdf")
		return p
	end
	Void
end

function plot_contour_unemp(h::Hank, savedir::String)
	p1, p2, p3, p4, pu, pY, pT, pTY = plot_state_funcs(h, MV=true)

	pl = style_contour(pu, 2, slides = true)
	savejson(pl, savedir * "unemp_slides.json")
	try
		savefig(pl, savedir*"unemp_slides.pdf")
	catch
	end

	pl = style_contour(pu, 2, slides = false)
	savejson(pl, savedir * "unemp_paper.json")

	nothing
end

function plot_contour_debtprice(h::Hank, savedir::String)
	p, p1, p2 = contour_debtprice(h, MV=true)

	pl = style_contour(p, 2, slides = true)
	savejson(pl, savedir * "debtprice_slides.json")
	try
		savefig(pl, savedir*"debtprice_slides.pdf")
	catch
	end

	pl = style_contour(p, 2, slides = false)
	savejson(pl, savedir * "debtprice_paper.json")

	nothing
end

function contour_debtprice(h::Hank; run_number::Int64=0, MV::Bool=true)
	qᵍ_mat  = reshape(h.qᵍ, h.Nb, h.Nμ, h.Nσ, h.Nξ, h.Nζ, h.Nz)

	jshow_b, jshow_μ, jshow_σ, jshow_ξ, jshow_ζ, jshow_z = ceil(Int, h.Nb*0.8), ceil(Int, h.Nμ*0.1), ceil(Int, h.Nσ*1), 2, 1, ceil(Int, h.Nz*0.1)

	debtcolors = [ ["0.0", "rgb(165,0,38)"],
					   ["0.2", "rgb(215,48,39)"],
					   ["0.4", "rgb(244,109,67)"],
					   ["0.6", "rgb(253,174,97)"],
					   ["0.8", "rgb(254,224,144)"],
					   ["1.0", "rgb(254,254,254)"]]

	tickmin = minimum(qᵍ_mat[:,:,:,jshow_ξ,1,:])
	tickmax = maximum(qᵍ_mat[:,:,:,jshow_ξ,1,:])
	ctbz = contour(;
		x = h.bgrid, y = exp.(h.zgrid),
		z = qᵍ_mat[:, jshow_μ, jshow_σ, jshow_ξ, jshow_ζ, :],
		# contours_coloring="heatmap",
		contours_start=tickmin, contours_end=tickmax,
		colorbar_tick0 = 0., #colorbar_dtick=floor(Int, 1.0./5),
		# colorscale = debtcolors, 
		colorscale = "Reds", reversescale = true,
		colorbar_dtick=0.1, colorbar_xpad=14
		)
	pbz = plot(ctbz, Layout(;xaxis_title="𝐵", yaxis_title="𝑧"))

	if MV 
		itp_qᵍ = make_itp(h, h.qᵍ; agg=true)
		qg_mat, xgrid, ygrid = reeval_mat_MV(h, itp_qᵍ, jshow_b, jshow_ξ, jshow_z, lb = 0)
		xax, yax = "Mean", "Variance"
	else
		qg_mat = qᵍ_mat[jshow_b, :, :, jshow_ξ, jshow_ζ, jshow_z]
		xgrid, ygrid = h.μgrid, h.σgrid
		xax, yax = "μ", "σ"
	end

	ctμσ = contour(;
		x = xgrid, y = ygrid,
		z = qg_mat,
		# contours_coloring="heatmap",
		contours_start=tickmin, contours_end=tickmax,
		colorbar_tick0 = 0., #colorbar_dtick=floor(Int, 1.0./5),
		# colorscale = debtcolors, 
		colorscale = "Reds", reversescale = true,
		colorbar_dtick=0.1, colorbar_xpad=14
		)

	pμσ = plot(ctμσ, Layout(;xaxis_title=xax, yaxis_title=yax))

	p = [pbz pμσ]
	relayout!(p, title = "Price of Debt")

	return p, pbz, pμσ
end

function reeval_mat_MV(h::Hank, itp_obj, jb, jξ, jz; lb=-Inf, ub=Inf)
	lb < ub || throw(error("Must specify upper bound greater than lower bound"))
	m_min, v_min = unmake_logN(h.μgrid[1], h.σgrid[1])
	m_max, v_max = unmake_logN(h.μgrid[end], h.σgrid[end])

	itp_obj = extrapolate(itp_obj, Interpolations.Flat())
	
	N = max(4*h.Nμ, 4*h.Nσ)

	mgrid = range(m_min, m_max, length=N)
	vgrid = range(v_min, v_max, length=N)
	
	show_b, show_ξ, jζ = h.bgrid[jb], h.ξgrid[jξ], 1

	mat = zeros(N, N)
	for (jm, m) in enumerate(mgrid)
		for (jv, v) in enumerate(vgrid)
			μv, σv = make_logN(m, v)
			Y = itp_obj(show_b, μv, σv, show_ξ, jζ, jz)
			mat[jm, jv] = max(lb, min(ub, Y))
		end
	end
	
	return mat, mgrid, vgrid
end

function plot_eulereq(h::Hank; remote::Bool=false)
	ExpRealRet = zeros(h.Ns, h.Nξ, h.Nz, 2)
	ExpExpRealRet = zeros(h.Ns, h.Nξ, h.Nz)
	probs = zeros(h.Ns, h.Nξ, h.Nz)
	ExpTRet = zeros(h.Ns, h.Nξ, h.Nz, 2)
	Exp_pC = zeros(h.Ns, h.Nξ, h.Nz, 2)
	EZ = zeros(h.Nω, h.Nϵ, 1, h.Nϵ, h.Nξ, h.Nz, 2)
	EIS = zeros(h.Nω, h.Nϵ, 1, h.Nϵ, h.Nξ, h.Nz, 2)

	jshow_b, jshow_μ, jshow_σ, jshow_ξ, jshow_ζ, jshow_z = ceil(Int, h.Nb*0.2), ceil(Int, h.Nμ/2), ceil(Int, h.Nσ/2), 2, 1, ceil(Int, h.Nz/2)

	jshow_ω, jshow_ϵ = ceil(Int, h.Nω*0.75), ceil(Int, h.Nϵ*0.5)

	jshow_s = findfirst((h.Jgrid[:,1].==jshow_b) .*
						(h.Jgrid[:,2].==jshow_μ) .* 
						(h.Jgrid[:,3].==jshow_σ) .*
						(h.Jgrid[:,4].==jshow_ξ) .*
						(h.Jgrid[:,5].==jshow_ζ) .*
						(h.Jgrid[:,6].==jshow_z)
						)

	pC_vec = price_index(h, h.pN)

	itp_qᵍ = make_itp(h, h.qᵍ, agg=true)
	itp_pC = make_itp(h, pC_vec, agg=true)
	itp_ϕc = make_itp(h, h.ϕc, agg=false)
	itp_vf = make_itp(h, h.vf, agg=false)

	rep_mat = reshape(h.repay, h.Nb, h.Nμ, h.Nσ, h.Nξ, h.Nζ, h.Nz, h.Nξ, h.Nz)

	for (js, js_show) in enumerate(1:size(h.Jgrid, 1))
		jb = h.Jgrid[js_show, 1]
		jμ = h.Jgrid[js_show, 2]
		jσ = h.Jgrid[js_show, 3]
		jξ = h.Jgrid[js_show, 4]
		jζ = h.Jgrid[js_show, 5]
		jz = h.Jgrid[js_show, 6]

		bp = h.issuance[js_show]
		μp = h.μ′[js_show,:,:,:]
		σp = h.σ′[js_show,:,:,:]

		pCv = price_index(h, h.pN[js_show])
		for (jξp, ξpv) in enumerate(h.ξgrid), jzp in 1:h.Nz
			prob = h.Pξ[jξ, jξp] * h.Pz[jz, jzp]
			# In repayment
			bpv = bp
			μpv = μp[jξp, jzp, 1]
			σpv = σp[jξp, jzp, 1]
			Rb = h.κ + (1.0 .-h.ρ) * itp_qᵍ[bpv, μpv, σpv, ξpv, 1, jzp]
			ExpRealRet[js, jzp, 1] = Rb * pCv / itp_pC[bpv, μpv, σpv, ξpv, 1, jzp] * prob / h.qᵍ[js]
			ExpExpRealRet[js, jzp] += ExpRealRet[js, jzp, 1] * rep_mat[jb, jμ, jσ, jξ, jζ, jz, jξp, jzp]
			ExpTRet[js, jzp, 1] = Rb * prob
			Exp_pC[js, jzp, 1] = pCv / itp_pC[bpv, μpv, σpv, ξpv, 1, jzp] * prob

			# In default
			haircut = (1.0 .-h.ℏ*(jζ==1))
			bpv = haircut * bp
			μpv = μp[jξp, jzp, 2]
			σpv = σp[jξp, jzp, 2]
			Rb = (1.0 .-h.ρ) * haircut * itp_qᵍ[bpv, μpv, σpv, ξpv, 2, jzp]
			ExpRealRet[js, jzp, 2] = Rb * pCv / itp_pC[bpv, μpv, σpv, ξpv, 2, jzp] * prob / h.qᵍ[js]
			ExpExpRealRet[js, jzp] += ExpRealRet[js, jzp, 2] * (1.0 .-rep_mat[jb, jμ, jσ, jξ, jζ, jz, jξp, jzp])
			ExpTRet[js, jzp, 2] = Rb * prob
			Exp_pC[js, jzp, 2] = pCv / itp_pC[bpv, μpv, σpv, ξpv, 2, jzp] * prob
		end
	end
	for (js, js_show) in enumerate(jshow_s)
		jb = h.Jgrid[js_show, 1]
		jμ = h.Jgrid[js_show, 2]
		jσ = h.Jgrid[js_show, 3]
		jξ = h.Jgrid[js_show, 4]
		jζ = h.Jgrid[js_show, 5]
		jz = h.Jgrid[js_show, 6]
		
		bp = h.issuance[js_show]
		μp = h.μ′[js_show,:,:,:]
		σp = h.σ′[js_show,:,:,:]

		for (jω, ωv) in enumerate(h.ωgrid)
			for jϵ in 1:h.Nϵ
				Tvf = 0.
				V = zeros(h.Nϵ, h.Nξ, h.Nz, 2)
				Cv = h.ϕc[jω, jϵ, jb, jμ, jσ, jξ, jζ, jz]
				Vf = h.vf[jω, jϵ, jb, jμ, jσ, jξ, jζ, jz]

				A = h.ϕa[jω, jϵ, jb, jμ, jσ, jξ, jζ, jz]
				B = h.ϕb[jω, jϵ, jb, jμ, jσ, jξ, jζ, jz]
				for (jξp, ξpv) in enumerate(h.ξgrid), jzp in 1:h.Nz
					rep_prob = rep_mat[jb, jμ, jσ, jξ, jζ, jz, jξp, jzp] * (jζ == 1) + h.θ * (jζ == 2)

					# First in repayment
					bpv = bp
					μpv = μp[jξp, jzp, 1]
					σpv = σp[jξp, jzp, 1]
					R = h.κ + (1.0 .-h.ρ) * itp_qᵍ[bpv, μpv, σpv, ξpv, 1, jzp]
					
					ωpv = A + R * B
					for jϵp in 1:h.Nϵ
						V_t = itp_vf[ωpv, jϵp, bpv, μpv, σpv, ξpv, 1, jzp]
						V[jϵp, jξp, jzp, 1] = V_t
						EIS[jω, jϵ, js, jϵp, jzp, 1] = (itp_ϕc[ωpv, jϵp, bpv, μpv, σpv, ξpv, 1, jzp] / Cv)^(-1.0./h.ψ)
						Tvf += V_t^(1.0 .-h.γ) * h.Pϵ[jϵ, jϵp] * h.Pz[jz, jzp] * rep_prob
					end
					
					# Then in default
					haircut = (1.0 .-h.ℏ*(jζ==1))
					bpv = haircut * bp
					μpv = μp[jξp, jzp, 2]
					σpv = σp[jξp, jzp, 2]
					R = h.κ + (1.0 .-h.ρ) * haircut * itp_qᵍ[bpv, μpv, σpv, ξpv, 2, jzp]
					
					ωpv = A + R * B
					for jϵp in 1:h.Nϵ
						V_t = itp_vf[ωpv, jϵp, bpv, μpv, σpv, ξpv, 2, jzp]
						V[jϵp, jξp, jzp, 2] = V_t
						EIS[jω, jϵ, js, jϵp, jzp, 2] = (itp_ϕc[ωpv, jϵp, bpv, μpv, σpv, ξpv, 2, jzp] / Cv)^(-1.0./h.ψ)
						Tvf += V_t^(1.0 .-h.γ) * h.Pϵ[jϵ, jϵp] * h.Pz[jz, jzp] * (1.0 .-rep_prob)
						
						EZ[jω, jϵ, js, jϵp, jξp, jzp, 1] = (V[jϵp, jξp, jzp, 1] ./ Tvf).^(1.0./h.ψ - h.γ)
						EZ[jω, jϵ, js, jϵp, jξp, jzp, 2] = (V[jϵp, jξp, jzp, 2] ./ Tvf).^(1.0./h.ψ - h.γ)
					end
				end
			end
		end
	end
	SDF = EZ .* EIS

	ESDF = zeros(h.Nω, h.Nϵ, h.Nξ, h.Nz, 2)
	for jζp in 1:2, jzp in 1:h.Nz, jϵp in 1:h.Nϵ, js in 1:1, jϵ in 1:h.Nϵ, jω in 1:h.Nω
		jz = h.Jgrid[jshow_s, 6]
		ESDF[jω, jϵ, jξp, jzp, jζp] += SDF[jω, jϵ, js, jϵp, jξp, jzp, jζp] * h.Pϵ[jϵ,jϵp] * h.Pξ[jξ, jξp] * h.Pz[jz, jzp]
	end

	p = plot([
		scatter(;x=h.zgrid, y=ExpRealRet[jshow_s, 1, :, 1], name="Ret in rep")
		scatter(;x=h.zgrid, y=ExpRealRet[jshow_s, 1, :, 2], name="Ret in def", line_dash = "dashdot")
		scatter(;x=h.zgrid, y=ExpExpRealRet[jshow_s, 1, :], name="Avg ret", line_dash = "dot", opacity = 0.75, line_width = 1)
		scatter(;x=h.zgrid, y=probs[jshow_s, 1, :], name="prob of z'", opacity=0.75)
		# scatter(;x=h.zgrid, y=ESDF[jshow_ω, jshow_ϵ, :, 1], name="SDF in rep")
		# scatter(;x=h.zgrid, y=ESDF[jshow_ω, jshow_ϵ, :, 2], name="SDF in def")
		])

	p2 = plot([
		scatter(;x=h.zgrid, y=ExpTRet[jshow_s, 1, :, 1], name="Ret_T in rep")
		scatter(;x=h.zgrid, y=ExpTRet[jshow_s, 1, :, 2], name="Ret_T in def")
		scatter(;x=h.zgrid, y=Exp_pC[jshow_s, 1, :, 1], line_dash="dashdot", name="pC/pC' in rep")
		scatter(;x=h.zgrid, y=Exp_pC[jshow_s, 1, :, 2], line_dash="dashdot", name="pC/pC' in def")
		])

	p = [p; p2]
end

function plot_defprob(h::Hank; remote::Bool=false)
	twisted_π = zeros(h.Nω, h.Nϵ, h.Nb, h.Nμ, h.Nσ, h.Nξ, h.Nζ, h.Nz)
	π_vec = zeros(h.Nω, h.Nϵ, h.Nb, h.Nμ, h.Nσ, h.Nξ, h.Nζ, h.Nz)

	itp_vf = make_itp(h, h.vf; agg=false)
	itp_qᵍ = make_itp(h, h.qᵍ, agg=true)

	for js in 1:size(h.Jgrid, 1)
		jb = h.Jgrid[js, 1]
		jμ = h.Jgrid[js, 2]
		jσ = h.Jgrid[js, 3]
		jξ = h.Jgrid[js, 4]
		jζ = h.Jgrid[js, 5]
		jz = h.Jgrid[js, 6]

		bp = h.issuance[js]
		μp = h.μ′[js,:,:,:]
		σp = h.σ′[js,:,:,:]
	
		rep_mat = reshape(h.repay, h.Nb, h.Nμ, h.Nσ, h.Nξ, h.Nζ, h.Nz, h.Nξ, h.Nz)

		for (jω, ωv) in enumerate(h.ωgrid), (jϵ, ϵv) in enumerate(h.ϵgrid)
			
			Tv = 0.0
			# First compute the denominator of the SDF
			ϕa = h.ϕa[jω, jϵ, jb, jμ, jσ, jξ, jζ, jz]
			ϕb = h.ϕb[jω, jϵ, jb, jμ, jσ, jξ, jζ, jz]

			for (jξp, ξpv) in enumerate(h.ξgrid), jzp in 1:h.Nz
		
				rep_prob = rep_mat[jb, jμ, jσ, jξ, jζ, jz, jξp, jzp]

				# In repayment
				bpv = bp
				μpv = μp[jξp, jzp, 1]
				σpv = σp[jξp, jzp, 1]

				Rb = h.κ + (1.0 .-h.ρ) * itp_qᵍ[bpv, μpv, σpv, ξpv, 1, jzp]


				ωpv = ϕa + Rb * ϕb

				for (jϵp, ϵpv) in enumerate(h.ϵgrid)
					prob = h.Pξ[jξ, jξp] * h.Pz[jz, jzp] * rep_prob * h.Pϵ[jϵ, jϵp]
					
					vfp = itp_vf[ωpv, jϵp, bpv, μpv, σpv, ξpv, 1, jzp]

					Tv += vfp^(1.0-h.γ) * prob
				end
				# In default
				bpv = (1.0-h.ℏ) * bp
				μpv = μp[jξp, jzp, 2]
				σpv = σp[jξp, jzp, 2]

				Rb = (1.0-h.ℏ) * itp_qᵍ[bpv, μpv, σpv, ξpv, 2, jzp]

				ωpv = ϕa + Rb * ϕb

				for (jϵp, ϵpv) in enumerate(h.ϵgrid)
					prob = h.Pξ[jξ, jξp] * h.Pz[jz, jzp] * (1.0-rep_prob) * h.Pϵ[jϵ, jϵp]
					
					vfp = itp_vf[ωpv, jϵp, bpv, μpv, σpv, ξpv, 2, jzp]

					Tv += vfp^(1.0-h.γ) * prob
				end
			end
			Tv = Tv^(1.0/(1.0-h.γ))

			# Now take the whole expectation to construct the likelihood ratio thing
			ϕa = h.ϕa[jω, jϵ, jb, jμ, jσ, jξ, jζ, jz]
			ϕb = h.ϕb[jω, jϵ, jb, jμ, jσ, jξ, jζ, jz]

			for (jξp, ξpv) in enumerate(h.ξgrid), jzp in 1:h.Nz
		
				rep_prob = rep_mat[jb, jμ, jσ, jξ, jζ, jz, jξp, jzp]

				# In default
				bpv = (1.0-h.ℏ) * bp
				μpv = μp[jξp, jzp, 2]
				σpv = σp[jξp, jzp, 2]

				Rb = (1.0-h.ℏ) * itp_qᵍ[bpv, μpv, σpv, ξpv, 2, jzp]

				ωpv_def = ϕa + Rb * ϕb

				bpv = (1.0-h.ℏ) * bp
				μpv = μp[jξp, jzp, 2]
				σpv = σp[jξp, jzp, 2]

				Rb = (1.0-h.ℏ) * itp_qᵍ[bpv, μpv, σpv, ξpv, 2, jzp]

				ωpv = ϕa + Rb * ϕb

				for (jϵp, ϵpv) in enumerate(h.ϵgrid)
					prob = h.Pξ[jξ, jξp] * h.Pz[jz, jzp] * (1.0-rep_prob) * h.Pϵ[jϵ, jϵp]
					
					vfp = itp_vf[ωpv_def, jϵp, bpv, μpv, σpv, ξpv, 1, jzp]

					twisted_π[jω, jϵ, jb, jμ, jσ, jξ, jζ, jz] += prob * (vfp/Tv)^(1.0/h.ψ - h.γ)
					π_vec[jω, jϵ, jb, jμ, jσ, jξ, jζ, jz] += prob
				end
			end
		end

	end

	return twisted_π, π_vec

end



function plot_aggcons(h::Hank; remote::Bool=false)
	jμ, jσ, jξ = ceil(Int, h.Nμ/2), ceil(Int, h.Nσ/2), ceil(Int, h.Nξ/2)
	μv, σv, ξv = h.μgrid[jμ], h.σgrid[jσ], h.ξgrid[jξ]
	jζ = 1

	itp_ϕc = make_itp(h, h.ϕc; agg=false)
	itp_ϕc2 = make_itp(h, h.ϕc.^2; agg=false)

	B′_mat = reshape(h.issuance, h.Nb, h.Nμ, h.Nσ, h.Nξ, h.Nζ, h.Nz)
	μ′_mat = reshape(h.μ′, h.Nb, h.Nμ, h.Nσ, h.Nξ, h.Nζ, h.Nz, h.Nz, 2)
	σ′_mat = reshape(h.σ′, h.Nb, h.Nμ, h.Nσ, h.Nξ, h.Nζ, h.Nz, h.Nz, 2)

	states = gridmake([1; h.Nb], [1; h.Nz])
	# p_vec = Array{PlotlyJS.SyncPlot{PlotlyJS.ElectronDisplay}}(size(states,1))
	p_vec = Array{PlotlyJS.SyncPlot}(undef, size(states,1))
	# p2_vec = Array{PlotlyJS.SyncPlot{PlotlyJS.ElectronDisplay}}(size(states,1))
	p2_vec = Array{PlotlyJS.SyncPlot}(undef, size(states,1))
	for js in 1:size(states,1)
		C_r = zeros(h.Nz)
		VarCr = zeros(h.Nz)
		C_d = zeros(h.Nz)
		VarCd = zeros(h.Nz)
		jb, jz = states[js, :]
		bvp = B′_mat[jb, jμ, jσ, jξ, jζ, jz]
		ξpv = h.ξgrid[1]
		for jzp in 1:h.Nz
			μvp = μ′_mat[jb, jμ, jσ, jξ, jζ, jz, jξp, jzp, 1]
			σvp = σ′_mat[jb, jμ, jσ, jξ, jζ, jz, jξp, jzp, 1]
			C_r[jzp] = integrate_itp(h, bvp, μvp, σvp, ξvp, 1, jzp, itp_ϕc)
			VarCr[jzp] = integrate_itp(h, bvp, μvp, σvp, ξvp, 1, jzp, itp_ϕc2) - C_r[jzp]^2
			μvp = μ′_mat[jb, jμ, jσ, jξ, jζ, jz, jξp, jzp, 2]
			σvp = σ′_mat[jb, jμ, jσ, jξ, jζ, jz, jξp, jzp, 2]
			C_d[jzp] = integrate_itp(h, (1.0 .-h.ℏ)*bvp, μvp, σvp, ξvp, 2, jzp, itp_ϕc)
			VarCd[jzp] = integrate_itp(h, (1.0 .-h.ℏ)*bvp, μvp, σvp, ξvp, 2, jzp, itp_ϕc2) - C_d[jzp]^2
		end
		p_vec[js] = plot(  [scatter(;x=h.zgrid, y=C_r, marker_color=col[1], showlegend=false),
						scatter(;x=h.zgrid, y=C_d, marker_color=col[4], showlegend=false, line_dash="dashdot")],
						Layout(;title="𝐵=$(h.bgrid[jb]), 𝑧=$(exp(h.zgrid[jz]))"))
		p2_vec[js] = plot(  [scatter(;x=h.zgrid, y=VarCr, marker_color=col[1], showlegend=false),
						scatter(;x=h.zgrid, y=VarCd, marker_color=col[4], showlegend=false, line_dash="dashdot")],
						Layout(;title="𝐵=$(h.bgrid[jb]), 𝑧=$(exp(h.zgrid[jz]))"))
	end

	p = [p_vec[1] p_vec[2]; p_vec[3] p_vec[4]]
	p.plot.layout["width"] = 800
	p.plot.layout["height"] = 800
	p.plot.layout["font_family"] = "Fira Sans Light"
	if remote
		path = pwd() * "/../../Graphs/"
		save(path * "p_aggcons.jld", "p", p)
	else
		path = pwd() * "/../Graphs/"
		# savefig(p, path * "aggcons.pdf")
		return p
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

function plot_state_funcs(h::Hank; remote::Bool=false, MV::Bool=true)

	pN_mat = reshape(h.pN,     h.Nb, h.Nμ, h.Nσ, h.Nξ, h.Nζ, h.Nz)
	w_mat  = reshape(h.wage,   h.Nb, h.Nμ, h.Nσ, h.Nξ, h.Nζ, h.Nz)
	u_mat  = 100*max.(1.0 .- reshape(h.Ld,     h.Nb, h.Nμ, h.Nσ, h.Nξ, h.Nζ, h.Nz),0)
	Y_mat  = reshape(h.output, h.Nb, h.Nμ, h.Nσ, h.Nξ, h.Nζ, h.Nz)
	Π_mat  = reshape(h.profits,h.Nb, h.Nμ, h.Nσ, h.Nξ, h.Nζ, h.Nz)
	g_mat  = reshape(h.spending,h.Nb, h.Nμ, h.Nσ, h.Nξ, h.Nζ, h.Nz)

	b_mat = copy(h.issuance)
	for js in 1:size(h.Jgrid, 1)
		bv = h.bgrid[h.Jgrid[js,1]]
		b_mat[js] = b_mat[js] - (1.0-h.ρ)*bv
	end
	b_mat  = reshape(b_mat, h.Nb, h.Nμ, h.Nσ, h.Nξ, h.Nζ, h.Nz)

	T_mat  = govt_bc(h, h.wage.*h.Ld)


	for (jp, jξ) in enumerate([1; h.Nξ])
		ppN1 = lines(h, pN_mat, 1, "Price of nontradables"; custom_w = jξ)
		pw1  = lines(h, w_mat, 1, "Wage"; custom_w = jξ)
		pY1  = lines(h, Y_mat, 1, "Output"; custom_w = jξ)
		pΠ1  = lines(h, Π_mat, 1, "Profits"; custom_w = jξ)
		pT1  = lines(h, T_mat, 1, "Taxes"; custom_w = jξ)
		pg1	 = lines(h, g_mat, 1, "Govt spending"; custom_w = jξ)
		pb1	 = lines(h, b_mat, 1, "Issuance"; custom_w = jξ)

		ppN2 = lines(h, pN_mat, 2; custom_w = jξ)
		pw2  = lines(h, w_mat, 2; custom_w = jξ)
		pY2  = lines(h, Y_mat, 2; custom_w = jξ)
		pΠ2  = lines(h, Π_mat, 2; custom_w = jξ)
		pT2  = lines(h, T_mat, 2; custom_w = jξ)
		pg2	 = lines(h, g_mat, 2; custom_w = jξ)
		pb2	 = lines(h, b_mat, 2; custom_w = jξ)


		ppN3 = lines(h, pN_mat, 3; custom_w = jξ)
		pw3  = lines(h, w_mat, 3; custom_w = jξ)
		pY3  = lines(h, Y_mat, 3; custom_w = jξ)
		pΠ3  = lines(h, Π_mat, 3; custom_w = jξ)
		pT3  = lines(h, T_mat, 3; custom_w = jξ)
		pg3	 = lines(h, g_mat, 3; custom_w = jξ)
		pb3	 = lines(h, b_mat, 3; custom_w = jξ)


		ppN4 = lines(h, pN_mat, 4; custom_w = jξ)
		pw4  = lines(h, w_mat, 4; custom_w = jξ)
		pY4  = lines(h, Y_mat, 4; custom_w = jξ)
		pΠ4  = lines(h, Π_mat, 4; custom_w = jξ)
		pT4  = lines(h, T_mat, 4; custom_w = jξ)
		pg4	 = lines(h, g_mat, 4; custom_w = jξ)
		pb4	 = lines(h, b_mat, 4; custom_w = jξ)


		ppN6 = lines(h, pN_mat, 6; custom_w = jξ)
		pw6  = lines(h, w_mat, 6; custom_w = jξ)
		pY6  = lines(h, Y_mat, 6; custom_w = jξ)
		pΠ6  = lines(h, Π_mat, 6; custom_w = jξ)
		pT6  = lines(h, T_mat, 6; custom_w = jξ)
		pg6	 = lines(h, g_mat, 6; custom_w = jξ)
		pb6	 = lines(h, b_mat, 6; custom_w = jξ)


		p1 = [ppN1 ppN2 ppN3 ppN4 ppN6; pw1 pw2 pw3 pw4 pw6] 
		p2 = [pΠ1 pΠ2 pΠ3 pΠ4 pΠ6; pT1 pT2 pT3 pT4 pT6]
		p3 = [pY1 pY2 pY3 pY4 pY6; pΠ1 pΠ2 pΠ3 pΠ4 pΠ6]
		p4 = [pg1 pg2 pg3 pg4 pg6; pb1 pb2 pb3 pb4 pb6]

		jshow_b, jshow_μ, jshow_σ, jshow_ξ, jshow_ζ, jshow_z = ceil(Int, h.Nb*0.5), ceil(Int, h.Nμ*0.75), ceil(Int, h.Nσ*0.5), 1, 1, ceil(Int, h.Nz*0.25)

		jshow_ξ = 2
		jshow_σ = 1

		function make_contour(y::Array; MV::Bool=true, cscale::String="Reds", perc::Bool=true)
			if MV 
				itp_y = make_itp(h, y; agg=true)
				lb, ub = minimum(y), maximum(y)
				yn_mat, xgrid, ygrid = reeval_mat_MV(h, itp_y, jshow_b, jshow_ξ, jshow_z, lb=lb, ub=ub)
				xax, yax = "Mean", "Variance"
			else
				yn_mat = y[jshow_b, :, :, jshow_ξ, jshow_ζ, jshow_z]
				xgrid, ygrid = h.μgrid, h.σgrid
				xax, yax = "μ", "σ"
			end
			
			tickmax = maximum(y[:,:,:,jshow_ξ,jshow_ζ,:])
			tickmin = minimum(y[:,:,:,jshow_ξ,jshow_ζ,:])

			# tickmin = 5.
			tickmax = 27.

			perc ? suffix = "%" : suffix = ""

			ctbz = contour(;
				x=h.bgrid, y=exp.(h.zgrid),
				z = y[:, jshow_μ, jshow_σ, jshow_ξ, jshow_ζ, :],
				# contours_coloring="heatmap",
				colorscale=cscale, contours_start=tickmin+0.01, contours_end=tickmax,
				colorbar_tick0 = 0., colorbar_dtick=floor(Int, tickmax/5),
				colorbar_ticksuffix=suffix, colorbar_showticksuffix="all"
				)
			ctμσ = contour(;
				x = xgrid, y = ygrid,
				z = yn_mat,
				# contours_coloring="heatmap",
				colorscale = cscale, contours_start=tickmin+0.01, contours_end=tickmax,
				colorbar_tick0 = 0., colorbar_dtick=floor(Int, tickmax/5),
				colorbar_ticksuffix=suffix, colorbar_showticksuffix="all"
				)
			pbz = plot(ctbz, Layout(;xaxis_title="𝐵", yaxis_title="𝑧"))	
			pμσ = plot(ctμσ, Layout(;xaxis_title=xax, yaxis_title=yax))
			
			p = [pbz pμσ]
		end

		pu = make_contour(u_mat)
		relayout!(pu, title = "Unemployment")

		pY = make_contour(Y_mat; perc=false)
		relayout!(pY, title = "Output")

		pT = make_contour(T_mat; perc=false)
		relayout!(pT, title = "Taxes")
		
		pTY = make_contour(T_mat./Y_mat * 100; perc=true)
		relayout!(pTY, title = "Taxes to GDP")
		
		return p1, p2, p3, p4, pu, pY, pT, pTY
	end
	nothing
end

function plot_LoM(h::Hank; remote::Bool=false)
	jz = ceil(Int, h.Nz/2)

	μ′_mat = zeros(h.Nb*h.Nμ*h.Nσ*h.Nξ*h.Nζ*h.Nz)
	σ′_mat = zeros(h.Nb*h.Nμ*h.Nσ*h.Nξ*h.Nζ*h.Nz)

	for js in 1:size(h.Jgrid, 1)
		jz = h.Jgrid[js, 6]

		μ′_mat[js] = h.μ′[js,jz,1]
		σ′_mat[js] = h.σ′[js,jz,1]
	end

	μ′_mat = reshape(μ′_mat, h.Nb, h.Nμ, h.Nσ, h.Nξ, h.Nζ, h.Nz)
	σ′_mat = reshape(σ′_mat, h.Nb, h.Nμ, h.Nσ, h.Nξ, h.Nζ, h.Nz)

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
		path = pwd() * "/../Graphs/"
		# savefig(p, path * "LoMs.pdf")
		return p
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
	l = scatter(;y=h.ξgrid, x=ones(h.ξgrid), line_dash="dashdot", marker_color="black", showlegend=false, mode="lines", title="Labor market")
	for (jpN, pNv) in enumerate(h.pngrid)
		Ld = labor_demand(h, h.ξgrid, exp(z_show), pNv)
		label = "pₙ = $(round(pNv,digits=2))"
		l = hcat(l, scatter(;y=h.ξgrid, x=Ld, name=label, marker_color=col[jpN], line_shape="spline"))
		if minimum(Ld) < vl
			vl = minimum(Ld)
		end
	end
	shapes = [hline(minimum(h.ξgrid), line_width=1)]
	layout = Layout(;	xaxis=attr(title="L", zeroline=true, range=[0., 3.]),
						yaxis=attr(title="w", zeroline=true),
						title="Labor Market",
						annotations=[attr(x=1, y=maximum(h.ξgrid),text="Lˢ", xanchor="center", yanchor="bottom", showarrow=false, font_size=18)],
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
		# savefig(p, path * "labordemand.pdf")
		return p
	end
	Void
end

function plot_nontradables(h::Hank; wrt::String="", remote::Bool=false)
	if wrt != "B" && wrt != "z"
		throw(error("Must specify B or z as wrt"))
	end
	jb = ceil(Int, h.Nb/2)
	jμ = ceil(Int, h.Nμ/2)
	jσ = ceil(Int, h.Nσ/2)
	jξ = ceil(Int, h.Nξ/2)
	jζ = ceil(Int, h.Nζ/2)
	jz = ceil(Int, h.Nz/2)

	bv, μv, σv, wv, ζv, zv = h.bgrid[jb], h.μgrid[jμ], h.σgrid[jσ], h.ξgrid[jξ], h.ζgrid[jζ], h.zgrid[jz]

	G_mat = reshape(h.spending, h.Nb, h.Nμ, h.Nσ, h.Nξ, h.Nζ, h.Nz)
	B_mat = reshape(h.issuance, h.Nb, h.Nμ, h.Nσ, h.Nξ, h.Nζ, h.Nz)

	itp_ϕc = make_itp(h, h.ϕc_ext; agg = false)

	pNmin, pNmax = minimum(h.pngrid), maximum(h.pngrid)

	l = Array{PlotlyBase.GenericTrace{Dict{Symbol,Any}}}(undef, 2*h.Nb)
	maxq = 0.
	minq = 10.
	for (jb, bv) in enumerate(h.bgrid)
		sup = zeros(size(h.pngrid))
		dem = zeros(size(h.pngrid))
		G   = G_mat[jb, jμ, jσ, jξ, jζ, jz]
		Bpv = B_mat[jb, jμ, jσ, jξ, jζ, jz]
		for (jpn, pnv) in enumerate(h.pngrid)
			sup[jpn], dem[jpn] = mkt_clearing(h, itp_ϕc, G, Bpv, pnv, pNmin, pNmax, bv, μv, σv, wv, jζ, jz, (jζ!=1); get_both=true)
		end
		l[jb] = scatter(; y=h.pngrid, x=sup, marker_color=col[jb], name="B = $(round(bv, digits=2))")
		l[h.Nb+jb] = scatter(; y=h.pngrid, x=dem, marker_color=col[jb], name="B = $(round(bv, digits=2))", showlegend=false)
		maxq = max(max(maximum(dem), maximum(sup)), maxq)
		minq = min(min(minimum(dem), minimum(sup)), minq)
	end
	maxq = min(maxq * 1.10, 3.)
	minq = minq * 0.9

	p = plot([l[jb] for jb in 1:2*h.Nb], Layout(; yaxis_title="pₙ", xaxis_title="Q", xaxis_range=[0., maxq]))
	if wrt == "B"
		if remote
			path = pwd() * "/../../Graphs/"
			save(path * "p_nontradables_B.jld", "p", p)
		else
			path = pwd() * "/../Graphs/"
			# savefig(p, path * "nontradables_B.pdf")
			return p
		end
	end

	jb = ceil(Int, h.Nb/2)
	bv, μv, σv, wv, ζv, zv = h.bgrid[jb], h.μgrid[jμ], h.σgrid[jσ], h.ξgrid[jξ], h.ζgrid[jζ], h.zgrid[jz]
	l = Array{PlotlyBase.GenericTrace{Dict{Symbol,Any}}}(undef, 2*h.Nz,2)
	maxq = 0.
	minq = 10.
	for (jz, zv) in enumerate(h.zgrid)
		sup = zeros(size(h.pngrid))
		dem = zeros(size(h.pngrid))
		supN = zeros(size(h.pngrid))
		G   = G_mat[jb, jμ, jσ, jξ, jζ, jz]
		Bpv = B_mat[jb, jμ, jσ, jξ, jζ, jz]
		for (jpn, pnv) in enumerate(h.pngrid)
			sup[jpn], dem[jpn] = mkt_clearing(h, itp_ϕc, G, Bpv, pnv, pNmin, pNmax, bv, μv, σv, wv, jζ, jz, (jζ!=1); get_both=true)

			zv = h.zgrid[jz]
			Ld, w_new, profits, output = labor_market(h, jζ, zv, pnv)
			Ld_N, _  = labor_demand(h, w_new, zv, jζ, pnv; get_both=true)
			supN[jpn] = TFP_N(zv, h.Δ, jζ) * Ld_N^(h.α_N)
		end
		l[jz,1] = scatter(; y=h.pngrid, x=sup, marker_color=col[ceil(Int,10*jz/h.Nz)], name="z = $(round(exp(zv), digits=2))")
		l[h.Nz+jz,1] = scatter(; y=h.pngrid, x=dem, marker_color=col[ceil(Int,10*jz/h.Nz)], name="z = $(round(exp(zv), digits=2))", showlegend=false)
		l[jz,2] = scatter(; x=supN, y=h.pngrid, marker_color=col[ceil(Int,10*jz/h.Nz)], name="z = $(round(exp(zv), digits=2))")
		maxq = max(max(maximum(dem), maximum(sup)), maxq)
		minq = min(min(minimum(dem), minimum(sup)), minq)
	end
	maxq = min(maxq * 1.10, 3.)
	minq = minq * 0.9

	p = plot([l[jz,1] for jz in 1:2*h.Nz], Layout(; yaxis_title="pₙ", xaxis_title="Q", xaxis_range=[0., maxq]))

	if wrt == "z"
		if remote
			path = pwd() * "/../../Graphs/"
			save(path * "p_nontradables_z.jld", "p", p)
			p = plot([l[jz,2] for jz in 1:h.Nz], Layout(;xaxis_title="Q", yaxis_title="pₙ", xaxis_range=[0., maxq]))
			save(path * "p_nontradables_z2.jld", "p", p)
		else
			path = pwd() * "/../Graphs/"
			# savefig(p, path * "nontradables_z.pdf")
			return p
		end
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
		# savefig(p, path * "conv.pdf")
		return p
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
		# savefig(p, path * "outconv.pdf")
		return p
	end
	Void
end

function plot_simul(path_entry::Path; remote::Bool=false, trim::Int=0)
	name = ""
	path = trim_path(path_entry, trim)
	if trim > 0
	else
		name = "_full"
	end

	T = size(path.data, 1)

	B_vec = series(path,:B)
	μ_vec = series(path,:μ)
	σ_vec = series(path,:σ)
	w_vec = series(path,:w)
	ζ_vec = series(path,:ζ) .- 1
	z_vec = exp.(series(path,:z))
	Y_vec = 4 * series(path,:Y)
	L_vec = series(path,:L)
	π_vec = series(path,:π)
	P_vec = series(path,:P)
	Pe_vec= series(path,:Pe)
	ψ_vec = series(path,:ψ)
	Bh_vec= series(path,:Bh)
	A_vec = series(path,:A)
	Bf_vec= series(path,:Bf)
	Wr_vec= series(path,:Wr)
	Wd_vec= series(path,:Wd)
	qᵍ_vec= series(path,:qg)
	ξ_vec = series(path,:ξ)
	NX_vec= series(path,:NX)

	m_vec, v_vec = unmake_logN(μ_vec, σ_vec)

	shiftζ = [0; ζ_vec[1:end-1]]

	defaults = findall((ζ_vec.==1) .* (shiftζ.==0))./4
	exits    = findall((ζ_vec.==0) .* (shiftζ.==1))./4

	times = (1:T)./4

	default_shades = rect(defaults, exits, 0, 1; fillcolor="#d3d3d3", opacity=0.5, line_width=0, xref="x", yref="paper")

	pB = plot([	scatter(; x=times, y=100*B_vec./Y_vec, marker_color=col[1], showlegend=false)
				#, scatter(; x=times, y=ones(times)*minimum(h.bgrid), showlegend=false, line_dash="dashdot", marker_color="black", line_width=0.5),
				#, scatter(; x=times, y=ones(times)*maximum(h.bgrid), showlegend=false, line_dash="dashdot", marker_color="black", line_width=0.5)
				],
						Layout(; title="Debt-to-GDP", xaxis=attr(title="𝑡"), yaxis_title="% of GDP"));
	pμ = plot([ scatter(; x=times, y=m_vec, marker_color=col[1], showlegend=false),
				# scatter(; x=times, y=ones(times)*minimum(h.μgrid), showlegend=false, line_dash="dashdot", marker_color="black", line_width=0.5),
				# scatter(; x=times, y=ones(times)*maximum(h.μgrid), showlegend=false, line_dash="dashdot", marker_color="black", line_width=0.5)
				],
						Layout(; title="Mean", yaxis_zeroline=false, xaxis=attr(title="𝑡")));
	pσ = plot([ scatter(; x=times, y=v_vec, marker_color=col[1], showlegend=false)
				# scatter(; x=times, y=ones(times)*maximum(h.σgrid), showlegend=false, line_dash="dashdot", marker_color="black", line_width=0.5),
				# scatter(; x=times, y=ones(times)*minimum(h.σgrid), showlegend=false, line_dash="dashdot", marker_color="black", line_width=0.5)
				],
						Layout(; title="Variance", yaxis_zeroline=false, xaxis=attr(title="𝑡")));
	pw = plot([ scatter(; x=times, y=w_vec, marker_color=col[1], showlegend=false)
				#, scatter(; x=times, y=ones(times)*minimum(h.ξgrid), showlegend=false, line_dash="dashdot", marker_color="black", line_width=0.5)
				#, scatter(; x=times, y=ones(times)*maximum(h.ξgrid), showlegend=false, line_dash="dashdot", marker_color="black", line_width=0.5)
				],
						Layout(; title="Wage", xaxis=attr(title="𝑡")));
	pz = plot(scatter(; x=times, y=z_vec, marker_color=col[1], showlegend=false), Layout(; title="TFP", xaxis=attr(title="𝑡")));
	pY = plot([ scatter(; x=times, y=Y_vec, z=100*(Y_vec .- mean(Y_vec)) / mean(Y_vec), marker_color=col[1], showlegend=false)
				#, scatter(; x=times, y=100 * (1.0 - L_vec), marker_color=col[2], showlegend=false, line_dash="dashdot")
			],
			Layout(; title="Output", yaxis_title="", xaxis=attr(title="𝑡")));
	pu = plot([
		scatter(; x=times, y=100*(1.0 .- L_vec), marker_color=col[1], showlegend=false)
		], Layout(; title="Unemployment", yaxis_title="%", xaxis_title="𝑡"));
	pπ = plot([scatter(; x=times, y=ζ_vec, marker_color=col[1], showlegend=false),
				scatter(; x=times, y=π_vec, marker_color=col[2], showlegend=false, line_dash="dashdot")],
			Layout(; title="Default prob", xaxis=attr(title="𝑡")));
	pP = plot([ scatter(; x=times, y=P_vec, marker_color=col[1], showlegend=false),
				# scatter(; x=times, y=ones(times)*maximum(h.pngrid), showlegend=false, line_dash="dashdot", marker_color="black", line_width=0.5),
				# scatter(; x=times, y=ones(times)*minimum(h.pngrid), showlegend=false, line_dash="dashdot", marker_color="black", line_width=0.5),
				scatter(; x=times, y=Pe_vec,marker_color=col[4], showlegend=false, line_dash="dashdot")],
			Layout(; title="Price of nontradables", xaxis=attr(title="𝑡")));
	pψ = plot(scatter(; x=times, y=100*ψ_vec, marker_color=col[1],  showlegend=false), Layout(; title="Fraction domestic", xaxis=attr(title="𝑡"), yaxis_title="% of total"));
	pA = plot(scatter(; x=times, y=A_vec, marker_color=col[1],  showlegend=false), Layout(; title="Domestic risk-free debt", xaxis_title="𝑡"));
	pBf= plot(scatter(; x=times, y=100*Bf_vec./Y_vec, marker_color=col[1], showlegend=false), Layout(; title="Foreign debt", xaxis_title="𝑡", yaxis_title="% of GDP"));
	pW = plot([ scatter(;x=times, y=Wr_vec, marker_color=col[1], showlegend=false),
				scatter(;x=times, y=Wd_vec, marker_color=col[2], showlegend=false, line_dash="dashdot")], Layout(;title="Welfare", xaxis_title="𝑡"));
	pBh = plot(scatter(; x=times, y=Bh_vec, marker_color=col[1], showlegend=false), Layout(;title="Domestic debt", xaxis_title="𝑡"))
	pqᵍ = plot(scatter(; x=times, y=qᵍ_vec, marker_color=col[1], showlegend=false), Layout(;title="Price of debt", xaxis_title="𝑡"))
	pξ = plot(scatter(; x=times, y=ξ_vec, marker_color=col[1], showlegend=false), Layout(;title="ξ", xaxis_title="𝑡"))
	pNX = plot(scatter(; x=times, y=NX_vec, marker_color=col[1], showlegend=false), Layout(;title="Net Exports", xaxis_title="𝑡"))

	p = [pB pw pz pY; pμ pσ pA pu; pψ pπ pW pP; pBh pqᵍ pξ pNX]
	# p.plot.layout["shapes"] = default_shades
	# p.plot.layout["width"] = 900
	# p.plot.layout["height"] = 700
	# p.plot.layout["font_family"] = "Fira Sans Light"

	name = "simul"*name
	if remote
		path = pwd() * "/../../Graphs/"
		save(path * "p_"*name*".jld", "p", p)
	else
		path = pwd() * "/../Graphs/"
		# savefig(p, path*name*".pdf")
		# savefig(p, path*name*".png")
		return p
	end

	nothing
end

function stats_sample(p::Path, sample; relative::Bool=false)
	println(size(sample))
	sample_stats = zeros(size(sample,1), size(sample,2), 4)
	Ndef = size(sample,3)

	if relative
		for jvar in 1:size(sample,1)
			for jep in 1:size(sample,3)
				if sample[jvar, 1, jep] > 0
					sample[jvar, :, jep] = sample[jvar, :, jep] / sample[jvar, 1, jep]
				end
			end
		end
	end
	
	for jvar in 1:size(sample,1)
		for jtr in 1:size(sample,2)
			sample_stats[jvar, jtr, 1] = quantile(vec(sample[jvar, jtr, :]), 0.25)
			sample_stats[jvar, jtr, 2] = quantile(vec(sample[jvar, jtr, :]), 0.50)
			sample_stats[jvar, jtr, 3] = quantile(vec(sample[jvar, jtr, :]), 0.75)
			sample_stats[jvar, jtr, 4] = mean(sample[jvar, jtr, :])
		end
	end

	return sample_stats
end

function stack_sample(p::Path, sample)
	Nepi = size(sample,3)
	Nt = size(sample,2)

	stack = zeros(size(sample,1), Nt*Nepi)
	for jepi in 1:Nepi
		stack[:, 1+(jepi-1)*Nt:jepi*Nt] = sample[:, :, jepi]
	end
	return stack
end

function volYC(p::Path; episode_type::String="default", πthres::Float64=0.975)
	sample, N, πv = find_episodes(p, episode_type=episode_type, πthres=πthres)

	stack = stack_sample(p, sample)

	σY = get_AR1(stack[p.n[:Y],:])[2]
	σC = get_AR1(stack[p.n[:C],:])[2]

	return σY, σC, σC/σY
end

function plot_episodes(p::Path; episode_type::String="default", slides::Bool=true, πthres::Float64=0.975, onlystats::Bool=true)

	sample, N, πv = find_episodes(p, episode_type=episode_type, πthres=πthres)
	iter = 0
	maxiter = 25
	while N == 0 && iter < maxiter
		iter += 1
		πthres *= 0.95
		sample, N, πv = find_episodes(p, episode_type=episode_type, πthres=πthres)
	end
	if iter == maxiter
		print_save("WARNING: No episodes of $(episode_type) found")
		return plot(), πthres
	else
		print_save("\nπthres = $πv, quantile $πthres")
	end

	sample_stats = stats_sample(p, sample)
	rel_sample_stats = stats_sample(p, sample; relative=true)

	function plot_sample(sym::Symbol, sampstats=sample_stats, sample=sample; onlystats::Bool=true, f::Function=identity, title::String="", yaxis_title::String="")
		ylow = sampstats[p.n[sym],:, 1]
		ymed = sampstats[p.n[sym],:, 2]
		yhig = sampstats[p.n[sym],:, 3]
		yavg = sampstats[p.n[sym],:, 4]

		data = Array{GenericTrace{Dict{Symbol,Any}},1}(undef,0)
		if onlystats
			push!(data, scatter(;x = -2.5:0.25:2.5, y = f(ylow), marker_color=col[1], line_dash="dot", opacity=0.25, showlegend=false, name="q=0.25", hoverinfo="skip"))
			push!(data, scatter(;x = -2.5:0.25:2.5, y = f(yhig), marker_color=col[1], line_dash="dot", opacity=0.25, fill="tonexty", fillcolor="rgba(31,119,180,0.1)", showlegend=false, name="q=0.75", hoverinfo="skip"))
				# scatter(;x = -2.5:0.25:2.5, y = f(yavg), marker_color=col[3], line_dash="dashdot", showlegend=false, opacity=0.4, name="mean")
		else
			for jj in 1:size(sample, 3)
				y = sample[p.n[sym], :, jj]
				scat = scatter(x = -2.5:0.25:2.5, y=f(y), marker_color=col[1], opacity = 0.25, showlegend=false, hoverinfo="skip")
				push!(data, scat)
			end
		end
		push!(data, scatter(;x = -2.5:0.25:2.5, y = f(ymed), marker_color=col[4], line_dash="solid", showlegend=false, name="q=0.5"))
		
		p1 = plot(data, Layout(;title=title, yaxis_title=yaxis_title, yaxis_zeroline=false, xaxis_range=[-2.5; 2.5]))
		return p1 
	end

	meanY = mean(p.data[:,p.n[:Y]])
	meanC = mean(p.data[:,p.n[:C]])
	meanμ = mean(p.data[:,p.n[:μ]]) 
	meanσ = mean(p.data[:,p.n[:σ]]) 
	pY = plot_sample(:Y, f=x->100*(x.-meanY)./meanY, title="Output", yaxis_title="% dev from mean", onlystats=onlystats)
	pu = plot_sample(:L, f=x->100*(1.0.-x), title="Unemployment", yaxis_title="%", onlystats=onlystats)
	pB = plot_sample(:B, f=x->100*x/(4*meanY), title="Bonds", yaxis_title="% of mean GDP", onlystats=onlystats)
	pG = plot_sample(:G, f=x->100*x/(meanY), title="Govt spending", yaxis_title="% of mean GDP", onlystats=onlystats)
	pT = plot_sample(:T, f=x->100*x/(meanY), title="Lump-sum taxes", yaxis_title="% of mean GDP", onlystats=onlystats)
	pμ = plot_sample(:mean, f=x->100*x/(4*meanY), title="Wealth Dist Mean", yaxis_title="% of mean GDP", onlystats=onlystats)
	pσ = plot_sample(:var, title="Wealth Dist Variance", onlystats=onlystats)
	pz = plot_sample(:z, f=x->100*(exp.(x).-1), title="TFP", onlystats=onlystats)
	pw = plot_sample(:w, title="Wage rate", onlystats=onlystats)
	pψ = plot_sample(:ψ, f=x->100*x, title="Proportion Domestic", yaxis_title="%", onlystats=onlystats)
	pπ = plot_sample(:π, f=x->100*x, title="Default prob", yaxis_title="%", onlystats=onlystats)
	pP = plot_sample(:P, title="Price of nontradables", onlystats=onlystats)
	pq = plot_sample(:qg, title="Price of new debt", onlystats=onlystats)
	pCf = plot_sample(:CoYd, f=x->100*x, title="C/Yᵈ", yaxis_title="%", onlystats=onlystats)
	pCl = plot_sample(:C, f=x->100*(x.-meanC)./meanC, title="Consumption", yaxis_title="% dev from mean", onlystats=onlystats)
	meanC = mean(rel_sample_stats[p.n[:C], 1, 4])
	pCs = plot_sample(:C, rel_sample_stats, f=x->100*(x.-meanC)./meanC, title="Consumption", yaxis_title="% dev from mean", onlystats=onlystats)
	meanY = mean(rel_sample_stats[p.n[:Y], 1, 4])
	pYs = plot_sample(:Y, rel_sample_stats, f=x->100*(x.-meanY)./meanY, title="Output", yaxis_title="% dev from mean", onlystats=onlystats)
	pp25 = plot_sample(:p25, rel_sample_stats, f=x->100x, title="bottom 25% holdings", onlystats=onlystats)
	pavgω = plot_sample(:avgω, title="Bondholdings-weighted avg wealth", onlystats=onlystats)
	pp90 = plot_sample(:p90, f=x->100*x, title="top 10% holdings", onlystats=onlystats)
	pNX = plot_sample(:NX, title="Net Exports", onlystats=onlystats)

	# p = [pz pY pCl pCf; pB pψ pq pπ; pu pμ pσ pP; pavgω pp90 pG pT]
	p = [pz pY pCl pP; pB pψ pπ pp90; pμ pσ pG pT]
	p = [pz pY pCl pP; pB pψ pπ pp90; pμ pσ pG pNX]
	slides ? font = "Fira Sans Light" : font = "STIX Two Text"
	relayout!(p, font_family = font, font_size = 18, height = 600, width = 1100)
	if slides
		relayout!(p, plot_bgcolor = "rgba(250, 250, 250, 1.0)", paper_bgcolor = "rgba(250, 250, 250, 1.0)")
	else
		relayout!(p, title="")
	end
	return p, πthres
end

function plot_comparison_episodes(path_bench::Path, path_nodef::Path, path_nodelta::Path=path_nodef; episode_type::String="default", slides::Bool=true, πthres::Float64=0.95, levels::Bool=true)

	t_epi, N, πv = find_times_episodes(path_bench; episode_type=episode_type, πthres=πthres)

	sample_bench = collect_episodes(path_bench, t_epi, N)
	sample_nodef = collect_episodes(path_nodef, t_epi, N)
	sample_nodelta = collect_episodes(path_nodelta, t_epi, N)

	samplestats_bench = stats_sample(path_bench, sample_bench)
	rel_samplestats_bench = stats_sample(path_bench, sample_bench; relative=true)

	samplestats_nodef = stats_sample(path_nodef, sample_nodef)
	rel_samplestats_nodef = stats_sample(path_nodef, sample_nodef; relative=true)

	samplestats_nodelta = stats_sample(path_nodelta, sample_nodelta)
	rel_samplestats_nodelta = stats_sample(path_nodelta, sample_nodelta; relative=true)

	function plot_sample(sym::Symbol, samplestats_bench=samplestats_bench, samplestats_nodef=samplestats_nodef, samplestats_nodelta=samplestats_nodelta; fb::Function=identity, fd::Function=fb, fΔ::Function=fb, title::String="", yaxis_title::String="")
		# ylow = sampstats[p.n[sym],:, 1]
		ybench = samplestats_bench[path_bench.n[sym],:, 2]
		ynodef = samplestats_nodef[path_nodef.n[sym],:, 2]
		ynodelta = samplestats_nodelta[path_nodelta.n[sym],:, 2]
		# yhig = sampstats[p.n[sym],:, 3]
		# yavg = sampstats[p.n[sym],:, 4]

		lines = [
			# scatter(;x = -2.5:0.25:2.5, y = f(ylow), marker_color=col[1], line_dash="dot", opacity=0.25, showlegend=false, name="q=0.1")
			# scatter(;x = -2.5:0.25:2.5, y = f(yhig), marker_color=col[1], line_dash="dot", opacity=0.25, fill="tonexty", fillcolor="rgba(31,119,180,0.1)", showlegend=false, name="q=0.9")
			scatter(;x = -2.5:0.25:2.5, y = fb(ybench), marker_color=col[1], line_dash="solid", showlegend=false, name="Benchmark")
			scatter(;x = -2.5:0.25:2.5, y = fd(ynodef), marker_color=col[2], line_dash="dashdot", showlegend=false, name="No default")
			# scatter(;x = -2.5:0.25:2.5, y = f(yavg), marker_color=col[3], line_dash="dashdot", showlegend=false, opacity=0.4, name="mean")
			]

		if ynodelta != ynodef
			push!(lines, 
				scatter(;x = -2.5:0.25:2.5, y = fΔ(ynodelta), marker_color=col[3], line_dash="dot", showlegend=false, name="Δ = 0")
				)
		end


		p1 = plot(lines, Layout(;title=title, yaxis_title=yaxis_title, yaxis_zeroline=false, xaxis_range=[-2.5; 2.5]))
		return p1 
	end

	meanY = mean(path_bench.data[:,path_bench.n[:Y]])
	meanYn = mean(path_nodef.data[:,path_nodef.n[:Y]])	
	meanYΔ = mean(path_nodelta.data[:,path_nodelta.n[:Y]])	
	meanC = mean(path_bench.data[:,path_bench.n[:C]])
	meanCd = mean(path_nodef.data[:,path_nodef.n[:C]])	
	meanCΔ = mean(path_nodelta.data[:,path_nodelta.n[:C]])	
	meanμ = mean(path_bench.data[:,path_bench.n[:μ]])
	meanσ = mean(path_bench.data[:,path_bench.n[:σ]])
	pY = plot_sample(:Y, fb=x->100*(x.-meanY)./meanY, fd = x->100*(x.-meanYn)/meanYn, fΔ = x->100*(x.-meanYΔ)/meanYΔ, title="Output", yaxis_title="% dev from mean")
	pu = plot_sample(:L, fb=x->100*(1.0.-x), title="Unemployment", yaxis_title="%")
	pB = plot_sample(:B, fb=x->100*x/(4*meanY), fd = x->100*x/(4*meanYn), fΔ = x->100*x/(4*meanYΔ), title="Bonds", yaxis_title="% of mean GDP")
	pG = plot_sample(:G, fb=x->100*x/(meanY), fd = x->100*x/(1*meanYn), fΔ = x->100*x/(1*meanYΔ), title="Govt spending", yaxis_title="% of mean GDP")
	pT = plot_sample(:T, fb=x->100*x/(meanY), fd = x->100*x/(1*meanYn), fΔ = x->100*x/(1*meanYΔ), title="Lump-sum taxes", yaxis_title="% of mean GDP")
	pz = plot_sample(:z, fb=x->100*(exp.(x).-1), title="TFP")
	pπ = plot_sample(:π, fb=x->100*x, title="Default prob", yaxis_title="%")
	pCl = plot_sample(:C, fb=x->100*(x.-meanC)./meanC, fd = x->100*(x.-meanCd)/meanCd, fΔ = x->100*(x.-meanCΔ)/meanCΔ, title="Consumption", yaxis_title="% dev from mean")
	pq = plot_sample(:qg, title="Price of new debt")
	pWr = plot_sample(:Wr, title="Welfare in repayment")

	if !levels
		pY = plot_sample(:Y, rel_samplestats_bench, rel_samplestats_nodef, rel_samplestats_nodelta, title="Output", yaxis_title="% dev from start")
		pCl = plot_sample(:C, rel_samplestats_bench, rel_samplestats_nodef, rel_samplestats_nodelta, title="Consumption", yaxis_title="% dev from start")
	end

	println(mean(samplestats_bench[path_bench.n[:Wr],:, 2]))
	println(mean(samplestats_nodef[path_nodef.n[:Wr],:, 2]))

	println((mean(samplestats_bench[path_bench.n[:Wr],:, 2]) - mean(samplestats_nodef[path_bench.n[:Wr],:, 2]))/mean(samplestats_bench[path_bench.n[:Wr],:, 2])*100)

	p1 = [pz pY pCl; pB pG pT; pu pπ pWr]
	slides ? font = "Fira Sans Light" : font = "STIX Two Text"
	# p1.plot.layout["font_family"] = font
	# p1.plot.layout["fontsize"] = 18
	# p1.plot.layout["height"] = 600
	# p1.plot.layout["width"] = 1100
	if slides
		p1.plot.layout["plot_bgcolor"] = "rgba(250, 250, 250, 1.0)"
		p1.plot.layout["paper_bgcolor"] = "rgba(250, 250, 250, 1.0)"
	else
		p1.plot.layout["title"] = ""
	end

	return p1

end


function make_IRF_plots(p::Path; slides::Bool=true, response::String="Y", impulse::String="z",
        verbose::Bool=false, create_plots::Bool=false, savedir::String="")
    zvec = series(p, :z)
    ξvec = series(p, :ξ)
    Wvec = series(p, :Wr) - series(p, :Wd)
    Cvec = log.(series(p, :C))
    Yvec = log.(series(p, :Y))
    πvec = series(p, :π)*100
    Bvec = series(p, :B)./series(p, :Y)
    t0 = 11
    T = length(zvec)

    if impulse == "z"
        shock_vec = zvec
        ρ_shock = 0.9
    elseif impulse == "ξ"
        shock_vec = ξvec
        ρ_shock = 0.95
    else
        throw(error("Keyword impulse has to be 'z' or 'ξ'"))
    end
    
    if response == "Y"
    	print_save("\nResponse variable Y")
        yvec = Yvec
    elseif response == "C"
    	print_save("\nResponse variable C")
        yvec = Cvec
    else
    	throw(error("No response variable selected"))
    end

    E_shock = shock_vec * ρ_shock
    y = yvec[t0:T]
    X = Bvec[t0:T]

    eps_z = [shock_vec[jt] - E_shock[jt-1] for jt in t0:T]

    Bmed = quantile(Bvec, 0.9)
    zlow = quantile(eps_z, 0.1)
    # eps_z = min.(eps_z, zlow)
    # eps_z = max.(eps_z, quantile(eps_z, 0.9))

    p1 = plot(
        scatter(;x=1:T, y=eps_z)
        , Layout(;shapes=[hline(minimum(eps_z)); hline(quantile(eps_z,0.01))])
    )

    H = 20
    β = Matrix{Float64}(undef, H+1,3)
    βhB = Matrix{Float64}(undef, H+1,6)
    βlz = Matrix{Float64}(undef, H+1,6)
    βboth = Matrix{Float64}(undef, H+1,3)
    for jh in 1:H+1
        yh = y[jh:end]
        # print_save("\nVariance of y = $(var(yh))")
        # print_save("\n$(yh)")
        Bh = X[1:end-jh+1]
        e_z = eps_z[1:end-jh+1]

        dummy_highB = (Bh.>Bmed)
        dummy_lowz = (e_z.<zlow)
        data = DataFrame(yh=yh, eps_z=e_z, X = Bh, ind_B = dummy_highB, ind_z = dummy_lowz)

        #data = data[data[:eps_z].<= zlow,:]
        # println(size(data[(data[:eps_z].<=zlow) .* (data[:X].>=Bmed),:]))
        OLS = glm(@formula(yh ~ eps_z), data, Normal(), IdentityLink())
        if verbose
            println(OLS)
        end
        β[jh,1] = coef(OLS)[2]
        β[jh,2] = coef(OLS)[2] + stderror(OLS)[2]*1.96
        β[jh,3] = coef(OLS)[2] - stderror(OLS)[2]*1.96
        #data_cond = data[(data[:X].>=Bmed),:]
        #OLS = glm(@formula(yh ~ eps_z), data_cond, Normal(), IdentityLink())    
        OLS = glm(@formula(yh ~ eps_z + eps_z&ind_B), data, Normal(), IdentityLink())
        if verbose
            println(OLS)
        end
        βhB[jh,1] = coef(OLS)[2]
        βhB[jh,2] = coef(OLS)[2] + stderror(OLS)[2]*1.96
        βhB[jh,3] = coef(OLS)[2] - stderror(OLS)[2]*1.96
        βhB[jh,4] = coef(OLS)[2] + coef(OLS)[3]
        βhB[jh,5] = coef(OLS)[2] + coef(OLS)[3] + stderror(OLS)[3]*1.96 + stderror(OLS)[2]*1.96
        βhB[jh,6] = coef(OLS)[2] + coef(OLS)[3] - stderror(OLS)[3]*1.96 - stderror(OLS)[2]*1.96
        #data_cond = data[(data[:X].<=Bmed),:]
        #OLS = glm(@formula(yh ~ eps_z), data_cond, Normal(), IdentityLink())    
        OLS = glm(@formula(yh ~ eps_z + eps_z&ind_z), data, Normal(), IdentityLink())
        if verbose
            println(OLS)
        end
        βlz[jh,1] = coef(OLS)[2]
        βlz[jh,2] = coef(OLS)[2] + stderror(OLS)[2]*1.96
        βlz[jh,3] = coef(OLS)[2] - stderror(OLS)[2]*1.96
        βlz[jh,4] = coef(OLS)[2] + coef(OLS)[3]
        βlz[jh,5] = coef(OLS)[2] + coef(OLS)[3] + stderror(OLS)[3]*1.96 #+ stderror(OLS)[2]*1.96
        βlz[jh,6] = coef(OLS)[2] + coef(OLS)[3] - stderror(OLS)[3]*1.96 #- stderror(OLS)[2]*1.96
        if verbose
            println(coef(OLS))
            println(βlz[jh,:])
        end
        #data_cond = data[(data[:eps_z].<zlow) .* (data[:X].>=Bmed),:]
        #OLS = glm(@formula(yh ~ eps_z), data_cond, Normal(), IdentityLink())    
        # OLS = glm(@formula(yh ~ eps_z + ind_B + eps_z*ind_B + ind_z + eps_z*ind_z + eps_z*ind_B*ind_z), data, Normal(), IdentityLink())
        # if verbose
        #     println(OLS)
        # end
        # βboth[jh,1] = coef(OLS)[2]
        # βboth[jh,2] = coef(OLS)[2] + stderror(OLS)[2]*1.96
        # βboth[jh,3] = coef(OLS)[2] - stderror(OLS)[2]*1.96
    end
    #β_h *= sqrt(var(zvec))
    
    yaxistitle = "∂log " * response * " / ∂log z"
    
    pYz = plot([
            scatter(;x=0:H, y=β[:,1], line_color=col[1], name="βₕ")        
            scatter(;x=0:H, y=β[:,3], line_width=0, showlegend=false, name="βl")
            scatter(;x=0:H, y=β[:,2], fill="tonexty", line_color="#bbd6e8", line_width=0, showlegend=false, name="βh")
        ], Layout(;xaxis_zeroline=false, yaxis_zeroline=false, xaxis_title="Quarters", yaxis_title=yaxistitle,
            legend_orientation="h", legend_x = 0.1, width = 600, height = 250, font_family = "STIX Two Text",
            shapes=[hline(0, line_dash="dot", line_width=1)]))

    pYzz = plot([
            scatter(;x=0:H, y=βlz[:,3], line_width=0, showlegend=false, name="βl")
            scatter(;x=0:H, y=βlz[:,2], fill="tonexty", line_color="#bbd6e8", line_width=0, showlegend=false, name="βh")
            scatter(;x=0:H, y=βlz[:,6], line_width=0, showlegend=false, name="βl")
            scatter(;x=0:H, y=βlz[:,5], fill="tonexty", line_color="#bfe2bf", line_width=0, showlegend=false, name="βh")
            scatter(;x=0:H, y=βlz[:,4], line_color=col[3], name="βₕ (low z)", line_dash="dashdot")        
            scatter(;x=0:H, y=βlz[:,1], line_color=col[1], name="βₕ")        
        ], Layout(;xaxis_zeroline=false, yaxis_zeroline=false, xaxis_title="Quarters", yaxis_title=yaxistitle,
            legend=attr(;orientation="h", x = 0.1, traceorder="reversed"), width = 600, height = 400, font_family = "STIX Two Text",
            shapes=[hline(0, line_dash="dot", line_width=1)]))

    pYzB = plot([
            scatter(;x=0:H, y=β[:,3], line_width=0, showlegend=false, name="βl")
            scatter(;x=0:H, y=β[:,2], fill="tonexty", line_color="#bbd6e8", line_width=0, showlegend=false, name="βh")
            scatter(;x=0:H, y=βhB[:,6], line_width=0, showlegend=false, name="βl")
            scatter(;x=0:H, y=βhB[:,5], fill="tonexty", line_color="#bfe2bf", line_width=0, showlegend=false, name="βh")
            scatter(;x=0:H, y=βhB[:,4], line_color=col[3], name="βₕ (high B)", line_dash="dashdot")        
            scatter(;x=0:H, y=β[:,1], line_color=col[1], name="βₕ")        
        ], Layout(;xaxis_zeroline=false, yaxis_zeroline=false, xaxis_title="Quarters", yaxis_title=yaxistitle,
            legend=attr(;orientation="h", x = 0.1, traceorder="reversed"), width = 600, height = 400, font_family = "STIX Two Text",
            shapes=[hline(0, line_dash="dot", line_width=1)]))
    
    if slides
        function style_plot_IRF!(pl)
        	relayout!(pl, font_family = "Fira Sans Light", font_size=16, plot_bgcolor="rgb(250, 250, 250, 1.0)", paper_bgcolor="rgb(250, 250, 250, 1.0")
        	nothing
        end
        style_plot_IRF!(pYz)
        style_plot_IRF!(pYzz)
        style_plot_IRF!(pYzB)
    end
    
    slides ? name = "_slides" : name = ""

    if create_plots
        savejson(pYz, savedir * "elast_" * response * "z" * name * ".json")
        savejson(pYzz, savedir * "elast_" * response * "zz" * name * ".json")
        savejson(pYzB, savedir * "elast_" * response * "zB" * name * ".json")
        try
	        savefig(pYz, savedir * "elast_" * response * "z" * name * ".pdf")
	        savefig(pYzz, savedir * "elast_" * response * "zz" * name * ".pdf")
	        savefig(pYzB, savedir * "elast_" * response * "zB" * name * ".pdf")
	    catch
	    end
    end
    nothing
end
