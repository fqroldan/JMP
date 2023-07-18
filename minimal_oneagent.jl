using QuantEcon, Distributions, ForwardDiff, Optim, PlotlyJS, ColorSchemes, ProgressBars

""" Define styles """

sand() = "#F5F3F1"
darkbgd() = "#272929"
lightgrid() = "#353535"
darkgrid() = "#e2e2e2"
gridcol(dark=false) = ifelse(dark, lightgrid(), darkgrid())

q_axis(dark) = attr(showgrid = true, gridcolor=gridcol(dark), gridwidth = 0.5, zeroline=false)
bgcol(slides, dark) = ifelse(slides, ifelse(dark, darkbgd(), sand()), "white")
qleg() = attr(orientation = "h", x=0.05, xanchor="left")

qwidth(slides) = 864
qheight(slides) = ceil(Int, qwidth(slides) * ifelse(slides, 10/16, 7/16))

function qtemplate(;dark=false, slides=!dark)
    axis = q_axis(dark)
    width = 864 #1920 * 0.45
    l = Layout(
        xaxis = axis, yaxis = axis,
        width = width,
        height = width * ifelse(slides, 10/16, 7/16),
        font = attr(
            family = ifelse(slides, "Lato", "Linux Libertine"),
            size = 16, color = ifelse(dark, sand(), darkbgd())
        ),
        paper_bgcolor = bgcol(slides, dark), plot_bgcolor = bgcol(slides, dark),
        legend = qleg(),
    )
    return Template(layout = l)
end

abstract type SOE
end

mutable struct SOEmin <: SOE
	pars::Dict{Symbol, Float64}

	zgrid::Vector{Float64}
	pz::Vector{Float64}
end

function SOEmin(;
	β=1.06^-1,	# Discount factor
	r=1.02-1,	# Risk-free rate
	
	γ=6,		# Risk aversion

	α=0.75, 		# Curvature of production function
	ϖ=0.55,			# Relative weight of nontradables
	η=1/0.83-1,		# Elasticity of substitution btw T and N
	Δ=.075,			# Productivity loss in default
	wbar=0.35,		# Wage rigidity

	Nz = 500,

	ξ_d = 1, 
)

	pars = Dict(:β => β, :r => r, :γ => γ, :α=>α, :ϖN => ϖ, :ϖT => 1-ϖ, :η => η, :wbar => wbar, :ξ_d => ξ_d)
	zgrid = range(0.01, 1.2, length=Nz)

	pz = pdf.(LogNormal(1, 0.75), zgrid)
	pz *= 1/sum(pz)

	return SOEmin(pars, zgrid, pz)
end

function eval_f_period2(sm::SOEmin, s1, d, Δ, f::Function)
	ξ_d = sm.pars[:ξ_d]

	u = 0.0
	sum_prob = 0.0
	for (jz, zv) in enumerate(sm.zgrid)

		c2 = 0.0
		# y2 = exp(zv)
		y2 = zv
		prob = sm.pz[jz]

		if y2 < d/Δ
			def = true
			c2 = y2 * (1-Δ) + s1
		else
			def = false
			c2 = y2 - ξ_d * d + s1
		end

		u += f(sm, c2) * prob
		sum_prob += prob
	end

	return sm.pars[:β] * u / sum_prob
end

eval_period2(sm::SOEmin, s1, d, Δ) = eval_f_period2(sm, s1, d, Δ, utility)
Euler_eq_RHS(sm::SOEmin, s1, d, Δ) = eval_f_period2(sm, s1, d, Δ, uprime) * (1+sm.pars[:r])


function CES_aggregator(sr::SOE, cT, cN)
	ϖN, ϖT, η = [sr.pars[key] for key in [:ϖN, :ϖT, :η]]

	return (ϖN^(1/η) * cN^(1-1/η) + ϖT^(1/η) * cT^(1-1/η))^(η/(η-1))
end

utility(sr::SOE, c) = utility(c, sr.pars[:γ])
function utility(c, γ)
	cmin = 1e-8
	if c < cmin
		return utility(cmin,γ) + (c-cmin) * (cmin)^-γ
	else
		γ == 1 && return log(c)
		return c^(1-γ)/(1-γ)
	end
end
uprime(sr::SOE, c) = ForwardDiff.derivative(x->utility(sr, x), c)

utility(sr::SOE, cT, cN) = utility(sr, CES_aggregator(sr, cT, cN))
∇u(sr::SOE, cT, cN) = ForwardDiff.gradient(x->utility(sr, x[1], x[2]), [cT,cN])

uT(sr::SOE, cT, cN) = ∇u(sr, cT, cN)[1]
uN(sr::SOE, cT, cN) = ∇u(sr, cT, cN)[2]


function prod_N(sr::SOE, h)
	""" Computes production of nontradables at input h """
	yN = h^sr.pars[:α]
	return yN
end
prod_N_prime(sr::SOE, h) = ForwardDiff.derivative(x->prod_N(sr,x), h)

function H(sr::SOE, cT, w)
	""" Computes labor supply consistent with consumption of tradables + wage """
	α, η, ϖN, ϖT = [sr.pars[key] for key in [:α, :η, :ϖN, :ϖT]]

	# return (ϖN/ϖT * α/w)^(1/(1+α*η)) * cT^(((1+η)/(1+α*η)))
	return (ϖN/ϖT * α/w)^(1/(1+α*η)) * cT^(1+η)
end

function eq_h(sr::SOE, cT)
	""" Computes labor supply consistent with consumption of tradables """
	Ls = 1

	h = H(sr, cT, sr.pars[:wbar])
	labor = min(h, Ls)
	return labor
end
H_prime(sr::SOE, cT) = ForwardDiff.derivative(x->eq_h(sr, x), cT)

function outcomes_period1(sm::SOEmin, s1, ω1; planner=false)
	cT = ω1 - s1 / (1+sm.pars[:r])
	h = eq_h(sm, cT)
	cN = prod_N(sm, h)

	if planner #&& h < 0.99
		μ = prod_N_prime(sm, h) * uN(sm, cT, cN)
	else
		μ = 0.0
	end
	return cT, h, cN, μ
end


function Euler_eq_LHS(sm::SOEmin, s1, ω1; planner=false)
	cT, h, cN, μ = outcomes_period1(sm, s1, ω1, planner=planner)

	LHS = uT(sm, cT, cN) + H_prime(sm, cT) * μ
end	

function def_prob(sm::SOEmin, d, Δ)

	Eid = 0.0
	sum_prob = 0.0
	for (jz, zv) in enumerate(sm.zgrid)
		# y2 = exp(zv)
		y2 = zv
		prob = sm.pz[jz]

		if y2 < d/Δ
			def = true
		else
			def = false
		end

		Eid += def * prob
		sum_prob += prob
	end

	Eid / sum_prob
end
q1(sm::SOEmin, d, Δ) = (1 - def_prob(sm, d, Δ)) / (1+sm.pars[:r])

function eval_period1(sm::SOEmin; s1, debt, Δ, y1=1, planner)
	# q1v = q1(sm, debt-s1, Δ)
	q1v = q1(sm, debt, Δ) * 0.0

	π = def_prob(sm, debt, Δ)

	cT, h, cN, μ = outcomes_period1(sm, s1, y1+q1v*debt, planner=planner)
	# cT, h, cN, μ = outcomes_period1(sm, 0.0, y1-q1v*s1, planner=planner)

	h = eq_h(sm, cT)
	h′ = H_prime(sm, cT)

	u1 = utility(sm, cT, cN)
	u2 = eval_period2(sm, s1, debt, Δ)
	# u2 = eval_period2(sm, 0.0, debt-s1, Δ)

	v = u1 + u2
	return (cT=cT, h=h, cN=cN, s_opt=s1, v=v, μ=μ*h′, q1v=q1v, debt=debt, π=π)
end

function solve_period1(sm::SOEmin; debt=0.2, Δ=0.4, y1=1, planner, repurchase=false)

	# obj_f(s1) = s1

	# if repurchase
	# 	obj_f(s1) = (Euler_eq_LHS(sm, 0.0, y1+q1(sm, debt-s1, Δ)*(-s1), planner=planner) - Euler_eq_RHS(sm, 0.0, debt-s1, Δ))^2
	# else
	# end
	obj_f(s1) = (Euler_eq_LHS(sm, s1, y1+q1v*debt; planner=planner) - Euler_eq_RHS(sm, s1, debt, Δ))^2
	q1v = q1(sm, debt, Δ) * 0.0

	res = Optim.optimize(obj_f, 0, y1/(1+sm.pars[:r]) - 1e-4, Brent())

	s_opt = res.minimizer
	dist = res.minimum
	
	dist < 1e-2 || print("WARNING: prob at debt $debt, Δ $Δ, planner $planner. Dist = $dist\n")
	
	eval_period1(sm, s1=s_opt, debt=debt, Δ=Δ, planner=planner)
end

function solve_opt_period1(sm::SOEmin; Δ=0.2, y1=1, planner=true, verbose=planner)

	# obj_f(X) = -eval_period1(sm, s1=X[2], debt=X[1], Δ=Δ, planner=planner).v + (Euler_eq_LHS(sm, X[2], y1+q1(sm, X[1], Δ)*X[1]; planner=planner) - Euler_eq_RHS(sm, X[2], X[1], Δ))^2

	# res = Optim.optimize(obj_f, [0.0, 0.0], [0.2, 0.15], [0.1, 0.1], Fminbox(NelderMead()))

	# return eval_period1(sm, debt=res.minimizer[1], s1 = res.minimizer[2], Δ=Δ, planner=planner)

	obj_f(debt) = -solve_period1(sm, debt=first(debt), Δ=Δ, planner=planner).v

	res = Optim.optimize(obj_f, 0, 0.15, Brent())

	# res = Optim.optimize(obj_f, [0.0], [0.9], [Δ], Fminbox(GradientDescent()))
	return solve_period1(sm, debt=res.minimizer, Δ=Δ, planner=planner)
end

function makeplots_minimal(sm::SOEmin; move, optim_debt = false, slides = true, dark = slides, template::Template=qtemplate(slides=slides, dark=dark), yh=1)

	Δvec = range(0,0.4,length=151)
	# Δvec = Δvec[2:end]

	cT, h, cN, s, v, μ, q, d, π = [zeros(length(Δvec), 2) for jj in 1:9]
	for ii in ProgressBar(eachindex(Δvec))
		for (jj, pl) in enumerate([true, false])
			if move == "d"
				cT[ii,jj], h[ii,jj], cN[ii,jj], s[ii,jj], v[ii,jj], μ[ii,jj], q[ii,jj], d[ii,jj], π[ii,jj] = solve_period1(sm, debt = Δvec[ii], planner = pl)
			elseif move == "Δ"
				cT[ii,jj], h[ii,jj], cN[ii,jj], s[ii,jj], v[ii,jj], μ[ii,jj], q[ii,jj], d[ii,jj], π[ii,jj] = solve_period1(sm, Δ = Δvec[ii], planner = pl)
			end
		end
	end

	annotations = [
		# attr(text="Value function", x=mean(Δvec), xref="x1", xanchor="center", y=1.1, font_size=18, yref="paper", showarrow=false)
		attr(text="Employment", x=mean(Δvec), xref="x2", xanchor="center", y=1.1, font_size=18, yref="paper", showarrow=false)
		attr(text="Savings", x=mean(Δvec), xref="x3", xanchor="center", y=1.1, font_size=18, yref="paper", showarrow=false)
		attr(text="Default probability", x=mean(Δvec), xref="x4", xanchor="center", y=1.1, font_size=18, yref="paper", showarrow=false)
		# [attr(text="<i>Δ", x=mean(Δvec), xref=xv, xanchor="center", y=-0.15, yref="paper", font_size=18, showarrow=false) for xv in ["x3", "x4"]]
	]

	data = [
		# [scatter(x=Δvec, y=v[:, jj], name=ifelse(pl, "Planner", "Decentralized"), legend_group=jj, line_color=get(ColorSchemes.southwest, (jj-1)), xaxis="x1", yaxis="y1") for (jj, pl) in enumerate([true, false])]
		[scatter(x=Δvec, y=cN[:, jj].^(1/sm.pars[:α]), name=ifelse(pl, "Planner", "Decentralized"), legend_group=jj, showlegend=true, line_color=get(ColorSchemes.southwest, (jj-1), :extrema), xaxis="x2", yaxis="y2") for (jj, pl) in enumerate([true, false])]
		[scatter(x=Δvec, y=s[:, jj], name=ifelse(pl, "Planner", "Decentralized"), legend_group=jj, showlegend=false, line_color=get(ColorSchemes.southwest, (jj-1), :extrema), xaxis="x3", yaxis="y3") for (jj, pl) in enumerate([true, false])]
		[scatter(x=Δvec, y=π[:, jj], name=ifelse(pl, "Planner", "Decentralized"), legend_group=jj, showlegend=false, line_color=get(ColorSchemes.southwest, (jj-1), :extrema), xaxis="x4", yaxis="y4") for (jj, pl) in enumerate([true, false])]
	]

	layout = Layout(
		template = template,
		annotations = annotations,
		# height = style.layout[:height] * yh,
		legend = attr(orientation="v", x=0.65, xanchor="right", y=0.1),
		# yaxis1=attr(domain=[0.575, 1], anchor="x1"),
		yaxis2=attr(domain=[0,1], anchor="x2"),
		yaxis3=attr(domain=[0,1], anchor="x3"),
		yaxis4=attr(domain=[0,1], anchor="x4"),
		# xaxis1 = attr(title="<i>Δ", anchor="y1", domain=[0, 0.45]),
		xaxis2 = attr(title="<i>"*move, anchor="y2", domain=[0, 0.3]),
		xaxis3 = attr(title="<i>"*move, anchor="y3", domain=[0.35, 0.65]),
		xaxis4 = attr(title="<i>"*move, anchor="y4", domain=[0.7, 1]),
		)

	if optim_debt 
		layout[:yaxis2_range] = [0.75, 1]
	end

	plot(data, layout)
end



function minimal_twoagents(γ = 4; slides = true, dark = slides, template::Template=qtemplate(slides=slides, dark=dark))

	πvec = range(0, 0.12, length=101);
	kvec = range(0, 0.16, length=101);

	c(k,π,γ) = 1 - 1/γ * log(1-π + π * exp(γ*k))

	layout = Layout(
		template = template,
		title = "Output in first period",
		xaxis_title = "Probability of transfer (%)",
		yaxis = attr(title = "Size of transfer (% of agg. consumption)", tick_padding=200),
		)

	plot(contour(x=100πvec, y=100kvec, z=100 * [c(k,π,γ)/c(0,0,γ) - 1 for π in πvec, k in kvec], line_width=0.1, contours=Dict(:start=>0,:end=>-2.5, :size=>0.125, :coloring=>"fill"), colorbar=attr(tick0=0, dtick=0.5)), layout)
end

function equil_period1(sm::SOEmin, f::Function; Δv, dv, planner)

	cT, h, cN, s, v, μ, q, d = solve_period1(sm, Δ=Δv, debt=dv, planner=false)

	if planner && h < 1
		cT, h, cN, s, v, μ, q, d = solve_period1(sm, Δ=Δv, debt=dv, planner=true)
	end

	eval_period1(sm, s1=s, debt=dv, Δ=Δv, planner=planner)
end