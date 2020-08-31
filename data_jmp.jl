using QuantEcon, CSV, DataFrames, Dates, GLM, PlotlyJS, ColorSchemes, ORCA, FixedEffectModels, RegressionTables, SparseArrays, ExcelReaders

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

function load_all()

	gdp_base = CSV.read("../Data/Eurostat aggregates/gdp_all.csv", missingstring=":")
	rates_base = CSV.read("../Data/Eurostat aggregates/rates_all.csv", missingstring=":")	
	debt_base = CSV.read("../Data/Eurostat aggregates/debt_all.csv", missingstring=":")
	# debt_base.Value[is.na(debt_base$Value)] = 0
	unemp_base = CSV.read("../Data/Eurostat aggregates/unemp_all.csv", missingstring=":")
	# unemp_base$Value[is.na(unemp_base$Value)] = 0
	exports_base = CSV.read("../Data/Eurostat aggregates/exports_all_namq_10_gdp.csv", missingstring=":")
	imports_base = CSV.read("../Data/Eurostat aggregates/imports_all_namq_10_gdp.csv", missingstring=":")
	g_base = CSV.read("../Data/Eurostat aggregates/g_all_namq_10_gdp.csv", missingstring=":")
	sav_base = CSV.read("../Data/Eurostat aggregates/nasq_10_ki_1_Data.csv", missingstring=":")
	c_base = CSV.read("../Data/Eurostat aggregates/HhC_all_namq_10_gdp_1_Data.csv", missingstring=":")
	coy_base = CSV.read("../Data/Eurostat aggregates/CoY_namq_10_gdp_1_Data.csv", missingstring=":")
	DeltaC_base = CSV.read("../Data/Eurostat aggregates/HHDC_namq_10_gdp_1_Data.csv", missingstring=":")

	for x in [gdp_base, rates_base, debt_base, unemp_base, exports_base, imports_base, g_base, sav_base, c_base, coy_base, DeltaC_base]
		
		x.TIME = Date.(["$(parse(Int, x.TIME[jt][1:4]))-$(parse(Int,x.TIME[jt][end])*3-2)" for jt in 1:length(x.TIME)], "yyyy-mm")
		
		sort!(x, [:GEO, :TIME])
	end

	TIME = gdp_base.TIME
	GEO = gdp_base.GEO
	gdp = gdp_base.Value
	rates = rates_base.Value
	debt = debt_base.Value
	unemp = unemp_base.Value
	exports = exports_base.Value
	imports = imports_base.Value
	g_spend = g_base.Value
	sav = sav_base.Value
	coy = coy_base.Value
	cons = parse.(Float64, replace.(c_base.Value, ',' => ""))

	data = DataFrame(TIME = TIME, GEO = GEO, debt=debt, gdp=gdp, rates=rates, unemp=unemp, exports=exports, cons = cons, imports=imports, coy=coy, g_spend=g_spend, sav=sav)

	data = data[data.GEO .!= "Greece",:]


	temp = CSV.read("../Data/fred_data.csv", missingstring=":")
	temp.TIME = Date.([temp.TIME[jt][1:6] * ifelse(parse(Int,temp.TIME[jt][7:8]) < 60, "20","19") * temp.TIME[jt][7:8] for jt in 1:length(temp.TIME)], "mm/dd/yyyy")
	temp = stack(temp, Not(:TIME))
	rename!(temp, :variable=>:GEO, :value=>:consnom)

	data = join(data, temp, on = [:GEO, :TIME])

	temp = CSV.read("../Data/fred_cpi.csv", missingstring=":")
	temp = stack(temp, Not(:TIME))
	rename!(temp, :variable=>:GEO, :value=>:cpi)

	data = join(data, temp, on = [:GEO, :TIME])

	data.cons = data.cons ./ data.cpi
	data.cons2 = data.consnom ./ data.cpi
	data.lcons = log.(data.cons)
	data.lgdp = log.(data.gdp)

	data.debt_level = data.debt .* data.gdp

	temp = CSV.read("../Data/fred_gdp.csv", missingstring=":")
	temp = stack(temp, Not(:TIME))
	rename!(temp, :variable=>:GEO, :value=>:gdp2)

	data = join(data, temp, on = [:GEO, :TIME])

	data.GEO[data.GEO.=="Germany (until 1990 former territory of the FRG)"] .= "Germany"

	# data.sav = 100 * (1 .- data.cons ./ data.gdp)

	data.NX = data.exports - data.imports
	data.NXsq = data.NX.^2

	data.year = year.(data.TIME)
	data.postDraghi = data.TIME .>= Date("2012-09-01")

	sort!(data, [:GEO, :TIME])
	data.debt_lag = by(data, :GEO, x = :debt => Base.Fix2(lag, 1)).x
	data.debt_lead = by(data, :GEO, x = :debt => Base.Fix2(lead, 1)).x
	data.debt_level_lead = by(data, :GEO, x = :debt_level => Base.Fix2(lead, 1)).x


	data.x = data.debt .* (data.TIME .== Date("2007-01-01"))
	temp = by(data, :GEO, b0=:x => x->maximum(skipmissing(x)))
	data = join(data, temp, on = [:GEO])

	data.x = data.rates .* (data.GEO .== "Germany (until 1990 former territory of the FRG)")
	temp = by(data, :TIME, GERint=:x => x->maximum(skipmissing(x)))
	data = join(data, temp, on = [:TIME])
	data.spread = data.rates - data.GERint


	data
end

function regs_sec2(data::DataFrame)

	t0 = Date("2008-01-01") # Initial debt
	t0 = Date("2008-01-01") 
	t1 = Date("2010-01-01") # Initial output
	t2 = Date("2013-01-01") # Final output

	df = data[(data.TIME .>= t1) .& (data.TIME .<= t2),:]
	# df = data	

	df.lcons_ = df.lcons
	# fs = reg(df, @formula(spread ~ b0 + fe(TIME)), save=true)
	# df.spr_hat = coef(fs)[1] * df.b0 + fe(fs).fe_TIME

	regY1 = reg(df, @formula(lgdp ~ spread + fe(GEO) + fe(TIME)))
	regC1 = reg(df, @formula(lcons ~ spread + fe(GEO) + fe(TIME)))

	regY2 = reg(df, @formula(lgdp ~ spread + debt + fe(GEO) + fe(TIME)))
	regC2 = reg(df, @formula(lcons ~ spread + debt + fe(GEO) + fe(TIME)))

	res1 = reg(df, @formula(lcons_ ~ spread + lgdp + fe(GEO) + fe(TIME)))
	res2 = reg(df, @formula(lcons_ ~ spread + debt + lgdp + fe(GEO) + fe(TIME)))

	regtable(regY1, regY2, regC1, regC2, res1, res2, renderSettings = latexOutput("../Output/current_best/table1.txt"), regression_statistics = [:nobs, :r2_within], print_estimator_section = false, labels = Dict("__LABEL_FE_YES__" => "\\checkmark", "__LABEL_STATISTIC_adjr2__" => "Adj.~\$R^2\$", "TIME"=>"Time FE", "GEO"=>"Country FE", "lgdp"=>"\$\\log Y_t\$", "lcons"=>"\$\\log C_t\$", "lcons_"=>"\$\\log C_{t}\$", "spread"=>"Spread\$_t\$", "debt"=>"\$B_t/Y_t\$"))
	regtable(regY1, regY2, regC1, regC2, res1, res2, regression_statistics = [:nobs, :adjr2, :r2_within], print_estimator_section = false, labels=Dict("TIME"=>"Time FE","GEO"=>"Country FE", "lgdp"=>"log Y_t", "lcons"=>"log C_t", "lcons_"=>" log C_t"))
end

function regs_fiscalrules(df::DataFrame; style::Style=slides_def, yh = 1)

	df = df[df.TIME .>= Date("2000-01-01"),:]

	df.net_iss = (df.debt_level_lead - (1-0.05) * df.debt_level) ./ df.gdp

	df.unemp2 = df.unemp.^2
	df.BoverY2 = df.debt.^2

	reg1G = reg(df, @formula(g_spend ~ unemp + unemp2 + debt + BoverY2 + NX + NXsq + fe(GEO) + fe(TIME)), save=true)
	reg2G = reg(df, @formula(g_spend ~ unemp + debt + NX + fe(GEO) + fe(TIME)), save=true)
	reg1B = reg(df, @formula(net_iss ~ unemp + unemp2 + debt + BoverY2 + NX + NXsq + fe(GEO) + fe(TIME)), save=true)
	reg2B = reg(df, @formula(net_iss ~ unemp + debt + NX + fe(GEO) + fe(TIME)), save=true)

	iss_hat = reg1B.fe.fe_GEO + reg1B.fe.fe_TIME
	g_hat   = reg1G.fe.fe_GEO + reg1G.fe.fe_TIME
	for jj in 1:length(reg1B.coefnames)
		iss_hat += reg1B.coef[jj] * df[!,reg1B.coefnames[jj]]
		g_hat += reg1G.coef[jj]  * df[!,reg1G.coefnames[jj]]
	end

	reg3G = lm(@formula(g_spend ~ unemp + debt + NX), df[df.GEO.=="Spain",:])
	reg4G = lm(@formula(g_spend ~ unemp + unemp2 + debt + BoverY2 + NX + NXsq), df[df.GEO.=="Spain",:])
	reg3B = lm(@formula(net_iss ~ unemp + debt + NX), df[df.GEO.=="Spain",:])
	reg4B = lm(@formula(net_iss ~ unemp + unemp2 + debt + BoverY2 + NX + NXsq), df[df.GEO.=="Spain",:])

	g_hat = predict(reg4G)
	iss_hat = predict(reg4B)

	println(regtable(reg1G, reg2G, reg3G, reg4G, reg1B, reg2B, reg3B, reg4B, renderSettings=asciiOutput(),regression_statistics = [:nobs, :r2, :adjr2, :r2_within], print_estimator_section = false, labels=Dict("TIME"=>"Time FE","GEO"=>"Country FE", "g_spend" => "G_t/Y_t", "unemp" => "Unemployment_t", "unemp2" => "Unemployment_t^2", "debt" => "B_t/Y_t", "BoverY2" => "(B_t/Y_t)^2", "NX" => "Net Exports_t", "NXsq" => "Net Exports_t^2", "net_iss" => "(B_t' - (1-ρ)B_t) / Y_t")))

	col = [	
		get(ColorSchemes.davos, 0.2)
		get(ColorSchemes.lajolla, 0.6)
		get(ColorSchemes.cork, 0.9)
	]

	data = [
		scatter(x=df[df.GEO.=="Spain",:].TIME, y=df[df.GEO.=="Spain",:].net_iss, name = "Observed", legendgroup = 1, yaxis="y1", marker_color=col[1])
		scatter(x=df[df.GEO.=="Spain",:].TIME, y=df[df.GEO.=="Spain",:].g_spend, name = "Observed", legendgroup = 1, yaxis="y2", showlegend=false, marker_color=col[1])
		scatter(x=df[df.GEO.=="Spain",:].TIME, y=iss_hat, name = "Fitted", line_dash="dashdot", legendgroup=2, yaxis="y1", marker_color=col[2])
		scatter(x=df[df.GEO.=="Spain",:].TIME, y=g_hat, name = "Fitted", line_dash="dashdot", legendgroup=2, yaxis="y2", showlegend=false, marker_color=col[2])
	]

	annotations = [
		attr(x=0.5, xref="paper", xanchor="center", y=1, yref="paper", text = "Government Spending", font_size=18, showarrow=false)
		attr(x=0.5, xref="paper", xanchor="center", y=0.45, yref="paper", text = "Debt Issuances", font_size=18, showarrow=false)
		]

	layout = Layout(annotations = annotations, height = style.layout[:height]*yh,
		yaxis1 = attr(domain = [0, 0.4], title="% of GDP"),
		yaxis2 = attr(domain = [0.55, 0.95], title="% of GDP"),
		)

	pg = plot([
		scatter(x=df[df.GEO.=="Spain",:].TIME, y=df[df.GEO.=="Spain",:].net_iss, name = "Observed")
		scatter(x=df[df.GEO.=="Spain",:].TIME, y=iss_hat, name = "Fitted")
		]
		)
	pb = plot([
		scatter(x=df[df.GEO.=="Spain",:].TIME, y=df[df.GEO.=="Spain",:].g_spend, name = "Observed")
		scatter(x=df[df.GEO.=="Spain",:].TIME, y=g_hat, name = "Fitted")
		]
		)

	plot(data, layout, style=style)
	# plot([
	# 	scatter(x=df[df.GEO.=="Spain",:].TIME, y=df[df.GEO.=="Spain",:].g_spend, name = "Observed")
	# 	scatter(x=df[df.GEO.=="Spain",:].TIME, y=g_hat, name = "Fitted")
	# 	]
	# 	)

	# plot([scatter(x=df[df.GEO.==country,:].net_iss, y=iss_hat[df.GEO.==country], name=country, text=df[df.GEO.==country,:].TIME, mode="markers") for country in unique(df.GEO)])

	# return coef(reg4G), coef(reg4B)
end

function hp_filter(y::Vector{Float64}, lambda::Float64=1600.)
	n = length(y)
	#@assert n >= 4

	diag2 = lambda*ones(n-2)
	diag1 = [ -2lambda; -4lambda*ones(n-3); -2lambda ]
	diag0 = [ 1+lambda; 1+5lambda; (1+6lambda)*ones(n-4); 1+5lambda; 1+lambda ]

	D = spdiagm(-2=>diag2, -1=>diag1, 0=>diag0, 1=>diag1, 2=>diag2)

	D\y
end

function hp_detrend(y::Vector{Float64}, lambda::Float64=1600.)
	return y - hp_filter(y, lambda)
end

function get_AR1(y::Vector)
	y_lag = y[1:end-1]
	y = y[2:end]

	data = DataFrame(yt = y, ylag = y_lag)
	OLS = glm(@formula(yt ~ ylag), data, Normal(), IdentityLink())

	println(OLS); ρ = coef(OLS)[2]

	ϵ = y - predict(OLS)

	σ = (sum(ϵ.^2)/(length(ϵ)))^.5

	return ρ, σ
end

function load_GDP_SPA()
	data_raw 		= readxl("../Data/Eurostat aggregates/Spain_gdp.xls", "Data2!B12:I103");

	d = readxl("../Data/Eurostat aggregates/Spain_gdp.xls", "Data2!A12:A103");
	dates_GDP 		= Date.(["$(parse(Int, d[jt][1:4]))-$(parse(Int,d[jt][end])*3-2)" for jt in 1:length(d)], "yyyy-mm")

	varnames = readxl("../Data/Eurostat aggregates/Spain_gdp.xls", "Data2!B11:I11")[:]

	varlabels = ["GDP", "C", "G", "C_hh", "C_npish", "Kform", "X", "M"]

	vardict = Dict(varnames[jj] => varlabels[jj] for jj in 1:length(varnames))

	df = DataFrame(convert(Matrix{Float64},data_raw))
	rename!(df, [jj => vardict[vn] for (jj, vn) in enumerate(varnames)])

	df.date = dates_GDP

	for yy in varlabels
		df[!,"$(yy)_hp"] = hp_detrend(log.(df[!,yy]))
	end

	df.CoverY = df.C_hh ./ df.GDP * 100

	return df
end

function load_SPA()
	df = load_GDP_SPA()

	data_rates = readxl("../Data/Eurostat aggregates/Spain_rates.xls", "Data!B10:D101")
	d2 = readxl("../Data/Eurostat aggregates/Spain_rates.xls", "Data!A10:A101")
	dates_rates = Date.(["$(parse(Int, d2[jt][1:4]))-$(parse(Int,d2[jt][end])*3-2)" for jt in 1:length(d2)], "yyyy-mm")
	spread = 100 * convert(Vector{Float64}, data_rates[:,2] - data_rates[:,1]);

	SPR = DataFrame("date" => dates_rates, "spread" => spread)

	df = join(df, SPR, on=["date"])

	unemp_base = CSV.read("../Data/Eurostat aggregates/Unemployment/une_rt_q_2_Data.csv", datarow=6)
	unemp_base.TIME = Date.(["$(parse(Int, unemp_base.TIME[jt][1:4]))-$(parse(Int,unemp_base.TIME[jt][end])*3-2)" for jt in 1:length(unemp_base.TIME)], "yyyy-mm")
	unemp_base = DataFrame(date = unemp_base.TIME, unemp = unemp_base.Value)
	df = join(df, unemp_base, on=["date"])

	debt_base = CSV.read("/home/q/Dropbox/Research/Active/Default_Inequality_Supply/empirics/Debt/data_spain.csv", datarow=2)[1:end-3,:]

	debt_base.Column1 = Date.(debt_base.Column1, "u-yy") .+ Year(2000) .- Month(2)
	rename!(debt_base, "Column1" => "date")

	debt_base = DataFrame([y=>debt_base[y] for y in ["date", "Total Foreign", "Total domestic", "Debt to GDP"]])

	return df
end
	
function SPA_CvY(;style::Style=slides_def, yh = 1)

	df = load_SPA()

	df = df[df.date.>= Date("2000-01-01"),:]

	jT = findfirst(df.date .>= Date("2008-01-01"))

	df.Output = 100 * df.GDP ./ df.GDP[jT]
	df.Consumption = 100 * df.C_hh ./ df.C_hh[jT]

	col = [	
		get(ColorSchemes.davos, 0.2)
		get(ColorSchemes.lajolla, 0.6)
		get(ColorSchemes.cork, 0.9)
	]

	lineattr = Dict("Output" => attr(color=col[1]), "Consumption" => attr(color=col[2], dash="dashdot"), "Spread" => attr(color=col[3], dash="dot"))

	date1Y = Date("2008-01-01")
	date1C = Date("2008-02-01")
	date2Y = Date("2013-04-01")
	date2C = Date("2013-06-01")
	datemid = Date("2010-10-01")

	maxY, maxC = 100, 100
	midY, midC = first(df.Output[df.date .>= datemid]), first(df.Consumption[df.date .>= datemid])
	minY, minC = df.Output[df.date .>= Date("2013-05-01")][1], df.Consumption[df.date .>= Date("2013-05-01")][1]

	shapes = [
			hline([minY], [Date(date1Y)], [Date(date2Y)], line = attr(width = 1.25, dash="dot", color=col[1]))
			hline([minC], [Date(date1C)], [Date(date2C)], line = attr(width = 1.25, dash="dot", color=col[2]))
			vline(date1Y, minY, 100, line = attr(width=1.5, dash="dot", color=col[1]))
			vline(date1C, minC, 100, line = attr(width=1.5, dash="dot", color=col[2]))
			vline(datemid-Dates.Month(1), minY, midY, line = attr(width=1.5, dash="dot", color=col[1]))
			vline(datemid, minC, midC, line = attr(width=1.5, dash="dot", color=col[2]))
		]
	annotations = [
			attr(x = Date(date1Y), y = 95, showarrow=false, xanchor="right", text="$(round(Int,100*(minY-maxY)/maxY))%", font_color=col[1])
			attr(x = Date(date1C), y = 92, showarrow=false, xanchor="left", text="$(round(Int,100*(minC-maxC)/maxC))%", font_color=col[2])
			attr(x = Date(datemid)-Dates.Month(1), y = 93, showarrow=false, xanchor="right", text="$(round(Int,100*(minY-midY)/midY))%", font_color=col[1])
			attr(x = Date(datemid), y = 88, showarrow=false, xanchor="right", text="$(round(Int,100*(minC-midC)/midC))%", font_color=col[2])
		]


	layout = Layout(shapes = shapes, annotations = annotations,
		height = style.layout[:height]*yh, 
		yaxis1 = attr(title="%"),
		yaxis2 = attr(titlefont_size=18, title="bps", overlaying = "y", side="right", zeroline=false),
		legend = attr(x=0.5, xanchor = "center"),
		)

	plot([
		[scatter(x=df.date, y=df[!,y], line = lineattr[y], name=y) for y in ["Output", "Consumption"]]
		scatter(x=df.date, y=df.spread, name="Spread (right axis)", line = lineattr["Spread"], yaxis="y2")
		], layout, style=style)
end

function SPA_nw(; style::Style=slides_def)

	function load_SPA_nw(k::String)
		data_raw = readxl("../Data/Eurostat aggregates/wealth_statistics_nasq_10_f_bs.xlsx", "Data$(k)!A12:F89");

		GEO = data_raw[1,:]
		GEO[1] = "date"

		data = DataFrame(data_raw[2:end, :])
		rename!(data, [jj => geo for (jj, geo) in enumerate(GEO)])
		data.date = Date.(["$(parse(Int, data.date[jt][1:4]))-$(parse(Int,data.date[jt][end])*3-2)" for jt in 1:length(data.date)], "yyyy-mm")

		df = DataFrame(date = data.date, assets = convert(Vector{Float64},data.Spain))
	end

	df1 = load_SPA_nw("")
	df2 = load_SPA_nw("2")
	rename!(df2, "assets" => "liabilities")

	df = join(df1, df2, on = "date")

	col = [	
		get(ColorSchemes.davos, 0.2)
		get(ColorSchemes.lajolla, 0.6)
		get(ColorSchemes.cork, 0.9)
	]

	annot = [
		attr(x = Date("2009-02"), y = mean(df.assets) - 1
			, ax = -40, ay = 40, xanchor = "right", yanchor = "top", text="Mean $(round(Int, mean(df.assets)))%");
		attr(x = Date("2004-01"), y = mean(df.liabilities) + 1
			, ax = 40, ay = -40, xanchor = "left", yanchor = "bottom", text="Mean $(round(Int, mean(df.liabilities)))%");
		]
	shapes = [
		hline(mean(df.assets), line = attr(width = 1, dash = "dot", color=col[1]))
		hline(mean(df.liabilities), line = attr(width = 1, dash = "dot", color=col[2]))
		]

	layout = Layout(annotations = annot, shapes = shapes,
		yaxis = attr(title="% of GDP"),
		)

	pNW = plot([
			scatter(;x = df.date, y=df[!,y], line_color=col[jj], name=uppercasefirst(y)) for (jj,y) in enumerate(["assets", "liabilities"])
			], layout, style=style)
end

