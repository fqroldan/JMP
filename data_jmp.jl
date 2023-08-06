using QuantEcon, CSV, DataFrames, Dates, GLM, PlotlyJS, ColorSchemes, FixedEffectModels, RegressionTables, SparseArrays, ExcelReaders, XLSX
 
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


function load_all(loaddir = "../Data/")

	gdp_base = CSV.read(loaddir * "/Eurostat aggregates/gdp_all.csv", DataFrame, missingstring=":")
	rates_base = CSV.read(loaddir * "/Eurostat aggregates/rates_all.csv", DataFrame, missingstring=":")	
	debt_base = CSV.read(loaddir * "/Eurostat aggregates/debt_all.csv", DataFrame, missingstring=":")
	# debt_base.Value[is.na(debt_base$Value)] = 0
	unemp_base = CSV.read(loaddir * "/Eurostat aggregates/unemp_all.csv", DataFrame, missingstring=":")
	# unemp_base$Value[is.na(unemp_base$Value)] = 0
	exports_base = CSV.read(loaddir * "/Eurostat aggregates/exports_all_namq_10_gdp.csv", DataFrame, missingstring=":")
	imports_base = CSV.read(loaddir * "/Eurostat aggregates/imports_all_namq_10_gdp.csv", DataFrame, missingstring=":")
	g_base = CSV.read(loaddir * "/Eurostat aggregates/g_all_namq_10_gdp.csv", DataFrame, missingstring=":")
	sav_base = CSV.read(loaddir * "/Eurostat aggregates/nasq_10_ki_1_Data.csv", DataFrame, missingstring=":")
	c_base = CSV.read(loaddir * "/Eurostat aggregates/HhC_all_namq_10_gdp_1_Data.csv", DataFrame, missingstring=":")
	coy_base = CSV.read(loaddir * "/Eurostat aggregates/CoY_namq_10_gdp_1_Data.csv", DataFrame, missingstring=":")
	DeltaC_base = CSV.read(loaddir * "/Eurostat aggregates/HHDC_namq_10_gdp_1_Data.csv", DataFrame, missingstring=":")

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


	temp = CSV.read(loaddir * "/fred_data.csv", DataFrame, missingstring=":")
	temp.TIME = Date.([temp.TIME[jt][1:6] * ifelse(parse(Int,temp.TIME[jt][7:8]) < 60, "20","19") * temp.TIME[jt][7:8] for jt in 1:length(temp.TIME)], "mm/dd/yyyy")
	temp = stack(temp, Not(:TIME))
	rename!(temp, :variable=>:GEO, :value=>:consnom)

	data = innerjoin(data, temp, on = [:GEO, :TIME])

	temp = CSV.read(loaddir * "/fred_cpi.csv", DataFrame, missingstring=":")
	temp = stack(temp, Not(:TIME))
	rename!(temp, :variable=>:GEO, :value=>:cpi)

	data = innerjoin(data, temp, on = [:GEO, :TIME])

	data.cons = data.cons ./ data.cpi
	data.cons2 = data.consnom ./ data.cpi
	data.lcons = log.(data.cons)
	data.lgdp = log.(data.gdp)

	data.debt_level = data.debt .* data.gdp

	temp = CSV.read(loaddir * "/fred_gdp.csv", DataFrame, missingstring=":")
	temp = stack(temp, Not(:TIME))
	rename!(temp, :variable=>:GEO, :value=>:gdp2)

	data = innerjoin(data, temp, on = [:GEO, :TIME])

	data.GEO[data.GEO.=="Germany (until 1990 former territory of the FRG)"] .= "Germany"

	# data.sav = 100 * (1 .- data.cons ./ data.gdp)

	data.NX = data.exports - data.imports
	data.NXsq = data.NX.^2

	data.year = year.(data.TIME)
	data.postDraghi = data.TIME .>= Date("2012-09-01")

	sort!(data, [:GEO, :TIME])
	data.debt0_temp = ifelse.(data.TIME .== Date("2007-01-01"), data.debt, missing)
	data.GERint_temp = ifelse.(data.GEO.=="Germany", data.rates, missing)
	gdf = groupby(data, :GEO)
	data.debt_lag = combine(gdf, :debt => Base.Fix2(lag, 1) => :x).x
	data.debt_lead = combine(gdf, :debt => Base.Fix2(lead, 1) => :x).x
	data.debt_level_lead = combine(gdf, :debt_level => Base.Fix2(lead, 1) => :x).x

	data.cpi_lag = combine(gdf, :cpi => Base.Fix2(lag, 4) => :x).x
	data.inflation = ((data.cpi ./ data.cpi_lag).^1 .- 1) * 100

	temp = combine(gdf, :debt0_temp => (x->maximum(skipmissing(x))) => :debt0)
	data = innerjoin(data, temp, on = [:GEO])

	sort!(data, [:TIME, :GEO])
	gdf = groupby(data, :TIME)
	temp = combine(gdf, :GERint_temp => (x->maximum(skipmissing(x))) => :GERint)

	data = innerjoin(data, temp, on = [:TIME])
	data.spread = 100 * ( data.rates - data.GERint )

	data
end

function regs_sec2(df_raw::DataFrame=load_all(); savedir="../Output/current_best/")

	t0 = Date("2008-01-01") # Initial debt
	t0 = Date("2008-01-01") 
	t1 = Date("2010-01-01") # Initial output
	t2 = Date("2013-01-01") # Final output

	df = df_raw[(df_raw.TIME .>= t1) .& (df_raw.TIME .<= t2),:]

	df.spread *= 0.01 # Measured in percent

	df.lcons_ = df.lcons
	# fs = reg(df, @formula(spread ~ b0 + fe(TIME)), save=true)
	# df.spr_hat = coef(fs)[1] * df.b0 + fe(fs).fe_TIME

	regY1 = reg(df, @formula(lgdp ~ spread + fe(GEO) + fe(TIME)))
	regC1 = reg(df, @formula(lcons ~ spread + fe(GEO) + fe(TIME)))
	regu1 = reg(df, @formula(unemp ~ spread + fe(GEO) + fe(TIME)))

	regY2 = reg(df, @formula(lgdp ~ spread + debt + fe(GEO) + fe(TIME)))
	regC2 = reg(df, @formula(lcons ~ spread + debt + fe(GEO) + fe(TIME)))
	regu2 = reg(df, @formula(unemp ~ spread + debt + fe(GEO) + fe(TIME)))

	res1 = reg(df, @formula(lcons_ ~ spread + lgdp + fe(GEO) + fe(TIME)))
	res2 = reg(df, @formula(lcons_ ~ spread + debt + lgdp + fe(GEO) + fe(TIME)))

	regtable(regY1, regY2, regC1, regC2, res1, res2, renderSettings = latexOutput(savedir*"table1.txt"), regression_statistics = [:nobs, :r2_within], print_estimator_section = false, labels = Dict("__LABEL_FE_YES__" => "\\checkmark", "__LABEL_STATISTIC_adjr2__" => "Adj.~\$R^2\$", "TIME"=>"Time FE", "GEO"=>"Country FE", "lgdp"=>"\$\\log Y_t\$", "lcons"=>"\$\\log C_t\$", "lcons_"=>"\$\\log C_{t}\$", "spread"=>"Spread\$_t\$", "debt"=>"\$B_t/Y_t\$"))

	regtable(regY1, regY2, regC1, regC2, regu1, regu2, res1, res2, regression_statistics = [:nobs, :adjr2, :r2_within], print_estimator_section = false, labels=Dict("TIME"=>"Time FE","GEO"=>"Country FE", "lgdp"=>"log Y_t", "lcons"=>"log C_t", "lcons_"=>" log C_t"))
end

function regs_fiscalrules(df::DataFrame; slides = true, dark = slides, template::Template=qtemplate(slides=slides, dark=dark), yh = 1, savedir = "../Output/current_best/")

    df = df[df.TIME.>=Date("2000-01-01"), :]

    # df = df[df.TIME .<= Date("2010-01-01"),:]

    df.net_iss = (df.debt_level_lead - (1 - 0.05) * df.debt_level) ./ df.gdp

    df.unemp2 = df.unemp .^ 2
    df.BoverY2 = df.debt .^ 2

    df.NXsq = df.NX .^ 2 .* sign.(df.NX)

    reg1G = reg(df, @formula(g_spend ~ unemp + unemp2 + debt + BoverY2 + NX + fe(GEO) + fe(TIME)), save = true)
    reg2G = reg(df, @formula(g_spend ~ unemp + debt + NX + fe(GEO) + fe(TIME)), save = true)
    reg1B = reg(df, @formula(net_iss ~ unemp + unemp2 + debt + BoverY2 + NX + fe(GEO) + fe(TIME)), save = true)
    reg2B = reg(df, @formula(net_iss ~ unemp + debt + NX + fe(GEO) + fe(TIME)), save = true)

    iss_hat = reg1B.fe.fe_GEO + reg1B.fe.fe_TIME
    g_hat = reg1G.fe.fe_GEO + reg1G.fe.fe_TIME
    for jj in 1:length(reg1B.coefnames)
        iss_hat += reg1B.coef[jj] * df[!, reg1B.coefnames[jj]]
        g_hat += reg1G.coef[jj] * df[!, reg1G.coefnames[jj]]
    end

    reg3G = lm(@formula(g_spend ~ unemp + debt + NX), df[df.GEO.=="Spain", :])
    reg4G = lm(@formula(g_spend ~ unemp + unemp2 + debt + BoverY2 + NX), df[df.GEO.=="Spain", :])
    reg3B = lm(@formula(net_iss ~ unemp + debt + NX), df[df.GEO.=="Spain", :])
    # reg4B = lm(@formula(net_iss ~ unemp + unemp2 + debt + BoverY2 + NX + NXsq), df[df.GEO.=="Spain",:])
    reg4B = lm(@formula(net_iss ~ unemp + unemp2 + debt + BoverY2 + NX), df[df.GEO.=="Spain", :])

    # g_hat = predict(reg4G)
    # iss_hat = predict(reg4B)
    g_hat = predict(reg3G)
    iss_hat = predict(reg3B)

    println(regtable(reg1G, reg2G, reg3G, reg4G, reg1B, reg2B, reg3B, reg4B, renderSettings = asciiOutput(), regression_statistics = [:nobs, :r2, :adjr2, :r2_within], print_estimator_section = false, labels = Dict("TIME" => "Time FE", "GEO" => "Country FE", "g_spend" => "G_t/Y_t", "unemp" => "Unemployment_t", "unemp2" => "Unemployment_t^2", "debt" => "B_t/Y_t", "BoverY2" => "(B_t/Y_t)^2", "NX" => "Net Exports_t", "NXsq" => "Net Exports_t^2", "net_iss" => "(B_t' - (1-ρ)B_t) / Y_t")))

    println(regtable(reg4G, reg3G, reg4B, reg3B, renderSettings = asciiOutput(), regression_statistics = [:nobs, :r2], print_estimator_section = false, labels = Dict("TIME" => "Time FE", "GEO" => "Country FE", "g_spend" => "G_t/Y_t", "unemp" => "Unemployment_t", "unemp2" => "Unemployment_t^2", "debt" => "B_t/Y_t", "BoverY2" => "(B_t/Y_t)^2", "NX" => "Net Exports_t", "NXsq" => "Net Exports_t^2", "net_iss" => "(B_t' - (1-ρ)B_t) / Y_t")))

    regtable(reg4G, reg3G, reg4B, reg3B, renderSettings = latexOutput(savedir * "table_fiscal.txt"), regression_statistics = [:nobs, :r2], print_estimator_section = false, labels = Dict("__LABEL_STATISTIC_adjr2__" => "Adj.~\$R^2\$", "__LABEL_STATISTIC_N__" => "Observations", "TIME" => "Time FE", "GEO" => "Country FE", "(Intercept)" => "Constant", "g_spend" => "\$G_t/Y_t\$", "unemp" => "Unemployment\$_t\$", "unemp2" => "Unemployment\$_t^2\$", "debt" => "\$B_t/Y_t\$", "BoverY2" => "\$(B_t/Y_t)^2\$", "NX" => "Net Exports\$_t\$", "NXsq" => "Net Exports\$_t^2\$", "net_iss" => "\$(B_t' - (1-ρ)B_t) / Y_t\$"))


    col = [
        get(ColorSchemes.davos, 0.2, :extrema)
        get(ColorSchemes.lajolla, 0.6, :extrema)
        get(ColorSchemes.cork, 0.9, :extrema)
    ]

    data = [
        scatter(x = df[df.GEO.=="Spain", :].TIME, y = df[df.GEO.=="Spain", :].net_iss, name = "Observed", legendgroup = 1, yaxis = "y1", marker_color = col[1])
        scatter(x = df[df.GEO.=="Spain", :].TIME, y = df[df.GEO.=="Spain", :].g_spend, name = "Observed", legendgroup = 1, yaxis = "y2", showlegend = false, marker_color = col[1])
        scatter(x = df[df.GEO.=="Spain", :].TIME, y = iss_hat, name = "Fitted", line_dash = "dashdot", legendgroup = 2, yaxis = "y1", marker_color = col[2])
        scatter(x = df[df.GEO.=="Spain", :].TIME, y = g_hat, name = "Fitted", line_dash = "dashdot", legendgroup = 2, yaxis = "y2", showlegend = false, marker_color = col[2])
    ]

    annotations = [
        attr(x = 0.5, xref = "paper", xanchor = "center", y = 1.025, yref = "paper", text = "Government Spending", font_size = 18, showarrow = false)
        attr(x = 0.5, xref = "paper", xanchor = "center", y = 0.45, yref = "paper", text = "Debt Issuances", font_size = 18, showarrow = false)
    ]

    layout = Layout(
		template = template,
		annotations = annotations, 
		# height = style.layout[:height] * yh,
        yaxis1 = attr(domain = [0, 0.4], title = "% of GDP"),
        yaxis2 = attr(domain = [0.55, 0.95], title = "% of GDP"),
    )

    pg = plot([
        scatter(x = df[df.GEO.=="Spain", :].TIME, y = df[df.GEO.=="Spain", :].net_iss, name = "Observed")
        scatter(x = df[df.GEO.=="Spain", :].TIME, y = iss_hat, name = "Fitted")
    ]
    )
    pb = plot([
        scatter(x = df[df.GEO.=="Spain", :].TIME, y = df[df.GEO.=="Spain", :].g_spend, name = "Observed")
        scatter(x = df[df.GEO.=="Spain", :].TIME, y = g_hat, name = "Fitted")
    ]
    )

    plot(data, layout)

    # return coef(reg4G), coef(reg4B)
    # return coef(reg3G), coef(reg3B)
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

function get_AR1(y::Vector; trend::Bool=false)
	y_lag = y[1:end-1]
	y = y[2:end]

	t = 1:length(y)

	data = DataFrame(yt = y, ylag = y_lag, t=t)
	if trend
		form = @formula(yt ~ ylag + t)
	else
		form = @formula(yt ~ ylag)
	end
	OLS = lm(form, data)
	
	ρ = coef(OLS)[2]
	# println(OLS); 
	# println("R2 = $(r2(OLS))")
	ϵ = y - predict(OLS)

	σ2 = sum(ϵ.^2)/(length(ϵ))
	σϵ = sqrt(σ2)

	σy = sqrt(σ2 / (1-ρ^2))
	return ρ, σϵ
end

function load_GDP_SPA()
	data_raw 		= readxl(loaddir * "/Eurostat aggregates/Spain_gdp.xls", "Data2!B12:I103");
	# data_raw = XLSX.readdata(loaddir * "/Eurostat aggregates/Spain_gdp.xls", "Data2", "B12:I103")

	d = readxl(loaddir * "/Eurostat aggregates/Spain_gdp.xls", "Data2!A12:A103");
	# d = XLSX.readdata(loaddir * "/Eurostat aggregates/Spain_gdp.xls", "Data2", "A12:A103")

	dates_GDP 		= Date.(["$(parse(Int, d[jt][1:4]))-$(parse(Int,d[jt][end])*3-2)" for jt in eachindex(d)], "yyyy-mm")

	varnames = readxl(loaddir * "/Eurostat aggregates/Spain_gdp.xls", "Data2!B11:I11")[:]
	# varnames = XLSX.readdata(loaddir * "/Eurostat aggregates/Spain_gdp.xls", "Data2", "B11:I11")[:]

	varlabels = ["GDP", "C", "G", "C_hh", "C_npish", "Kform", "X", "M"]

	vardict = Dict(varnames[jj] => varlabels[jj] for jj in eachindex(varnames))

	df = DataFrame(convert(Matrix{Float64},data_raw), :auto)
	rename!(df, [jj => vardict[vn] for (jj, vn) in enumerate(varnames)])

	df.date = dates_GDP

	for yy in varlabels
		df[!,"$(yy)_hp"] = hp_detrend(log.(df[!,yy]))
	end

	df.CoverY = df.C_hh ./ df.GDP * 100

	return df
end

function load_SPA(; loaddir = "../Data/")
	df = load_GDP_SPA()

	data_rates = readxl(loaddir * "/Eurostat aggregates/Spain_rates.xls", "Data!B10:D101")
	# data_rates = XLSX.readdata(loaddir * "/Eurostat aggregates/Spain_rates.xls", "Data", "B10:D101")
	d2 = readxl(loaddir * "/Eurostat aggregates/Spain_rates.xls", "Data!A10:A101")
	# d2 = XLSX.readdata(loaddir * "/Eurostat aggregates/Spain_rates.xls", "Data", "A10:A101")

	dates_rates = Date.(["$(parse(Int, d2[jt][1:4]))-$(parse(Int,d2[jt][end])*3-2)" for jt in eachindex(d2)], "yyyy-mm")
	spread = 100 * convert(Vector{Float64}, data_rates[:,2] - data_rates[:,1]);

	SPR = DataFrame("date" => dates_rates, "spread" => spread)

	df = innerjoin(df, SPR, on=["date"])

	unemp_base = CSV.read(loaddir * "/Eurostat aggregates/Unemployment/une_rt_q_2_Data.csv", DataFrame, datarow=6)
	unemp_base.TIME = Date.(["$(parse(Int, unemp_base.TIME[jt][1:4]))-$(parse(Int,unemp_base.TIME[jt][end])*3-2)" for jt in eachindex(unemp_base.TIME)], "yyyy-mm")
	unemp_base = DataFrame(date = unemp_base.TIME, unemp = unemp_base.Value)
	df = innerjoin(df, unemp_base, on=["date"])

	debt_base = CSV.read(loaddir * "/spain_debt/data_spain.csv", DataFrame, datarow=2)[1:end-3,:]

	debt_base.Column1 = Date.(debt_base.Column1, "u-yy") .+ Year(2000) .- Month(2)
	rename!(debt_base, "Column1" => "date")

	debt_base = DataFrame([y=>debt_base[!,y] for y in ["date", "Total Foreign", "Total domestic", "Debt to GDP"]])

	df = innerjoin(df, debt_base, on=["date"])
	return df
end

function SPA_targets()
	df = load_all()
	df = df[df.TIME .>= Date("2000-01-01"), :]
	df = df[df.GEO .== "Spain", :]
	rename!(df, "TIME" => "date")

	df_nw = load_nw("Spain")

	df = innerjoin(df, df_nw, on = ["date"])

	ρy, σy = get_AR1(log.(df.gdp))
	ρc, σc = get_AR1(log.(df.cons))
	ρq, σq = get_AR1(df.spread)

	B_avg = mean(df.debt)
	B_std = sqrt(var(df.debt))
	u_avg = mean(df.unemp)
	u_std = sqrt(var(df.unemp))

	w_avg = mean(df.net_worth)

	df_spa = load_SPA()
	median_dom = 100 * QuantEcon.quantile(df_spa[!,"Total domestic"] ./ ( df_spa[!,"Total domestic"] + df_spa[!,"Total Foreign"] ), 0.5)

	SPA_gini = CSV.read(loaddir * "/Gini/icw_sr_05_1_Data.csv", DataFrame)

	Gini = SPA_gini.Value[end]

	[ρy, 100*σy, ρc, 100*σc, ρq, σq, B_avg, B_std, u_avg, u_std, median_dom, w_avg, Gini]
end

function save_SPA_targets()
	df = DataFrame(:x => SPA_targets())

	CSV.write("SPA_targets.csv", df);
end

function SPA_comp(; datadir = "../Data/", init = Date("2010-01-01"))
	df = load_all(datadir)
	df = df[df.GEO.=="Spain",:]
	rename!(df, "gdp" => "GDP", "cons" => "C_hh", "TIME" => "date")

	df = df[df.date .>= Date("2008-01-01"),:]
	df = df[df.date .<= Date("2013-01-01"),:]

	GDP08 = df.GDP[df.date .== init]
	cons = df.C_hh
	cons08 = cons[df.date .== init]

	Yn = 100 * df.GDP ./ GDP08
	Cn = 100 * cons ./ cons08
	
	dfn = DataFrame(:date => df.date, :Y => Yn, :C => Cn, :spread => df.spread, :BoY => df.debt, :unemp => df.unemp, :GoY => df.g_spend)

	dfn = dfn[dfn.date .>= init,:]
	dfn
end

function SPA_CvY(country::String="Spain"; loaddir = "../Data/", slides = true, dark = slides, template::Template=qtemplate(slides=slides, dark=dark), yh = 1, sh=true)

	# df = load_SPA()
	# if country != "Spain"
	df = load_all(loaddir)
	df = df[df.GEO.==country,:]
	rename!(df, "gdp" => "GDP", "cons" => "C_hh", "TIME" => "date")
	# end

	df = df[df.date.>= Date("2000-01-01"),:]

	jT = findfirst(df.date .>= Date("2008-01-01"))

	df.Output = 100 * df.GDP ./ df.GDP[jT]
	df.Consumption = 100 * df.C_hh ./ df.C_hh[jT]

	col = [	
		get(ColorSchemes.davos, 0.2, :extrema)
		get(ColorSchemes.lajolla, 0.5, :extrema)
		get(ColorSchemes.cork, 0.9, :extrema)
	]

	lineattr = Dict("Output" => attr(width = 2.5, color=col[1]), "Consumption" => attr(width = 2.5, color=col[2], dash="dashdot"), "Spread" => attr(width = 2.5, color=col[3], dash="dot"))

	date1Y = Date("2008-01-01")
	date1C = Date("2008-02-01")
	date2Y = Date("2012-10-01")
	date2C = Date("2012-10-01")
	datemid = Date("2010-10-01")

	maxY, maxC = 100, 100
	midY, midC = first(df.Output[df.date .>= datemid]), first(df.Consumption[df.date .>= datemid])
	minY, minC = df.Output[df.date .>= date2Y][1], df.Consumption[df.date .>= date2C][1]

	shapes = []
	annotations = []

	if sh
		shapes = [
				hline([minY], [Date(date1Y)], [Date(date2Y)], line = attr(width = 1.5, dash="dot", color=col[1]))
				hline([minC], [Date(date1C)], [Date(date2C)], line = attr(width = 1.5, dash="dot", color=col[2]))
				vline(date1Y, minY, 100, line = attr(width=1.75, dash="dot", color=col[1]))
				vline(date1C, minC, 100, line = attr(width=1.75, dash="dot", color=col[2]))
				vline(datemid-Dates.Month(1), minY, midY, line = attr(width=1.75, dash="dot", color=col[1]))
				vline(datemid, minC, midC, line = attr(width=1.75, dash="dot", color=col[2]))
			]
		annotations = [
				attr(x = Date(date1Y), y = (maxY+minY)/2, showarrow=false, xanchor="right", text="$(round(Int,100*(minY-maxY)/maxY))%", font_color=col[1])
				attr(x = Date(date1C), y = (maxC+minC)/2, showarrow=false, xanchor="left", text="$(round(Int,100*(minC-maxC)/maxC))%", font_color=col[2])
				attr(x = Date(datemid)-Dates.Month(1), y = minY, yanchor="bottom", showarrow=false, xanchor="right", text="$(round(Int,100*(minY-midY)/midY))%", font_color=col[1])
				attr(x = Date(datemid), y = minC, yanchor="bottom", showarrow=false, xanchor="left", text="$(round(Int,100*(minC-midC)/midC))%", font_color=col[2])
			]
	end


	layout = Layout(
		template = template,
		shapes = shapes, annotations = annotations,
		title = "Output and Consumption in " * country,
		yaxis1 = attr(title="%"),
		yaxis2 = attr(titlefont_size=18, title="bps", overlaying = "y", side="right", zeroline=false),
		legend = attr(x=0.5, xanchor = "center"), hovermode="x",
		)

	plot([
		[scatter(x=df.date, y=df[!,y], line = lineattr[y], name=y) for y in ["Output", "Consumption"]]
		scatter(x=df.date, y=df.spread, name="Spread (right axis)", line = lineattr["Spread"], yaxis="y2")
		], layout)
end

function load_nw(country::String="Spain", loaddir ="../Data/")
	function load_nw(k::String, loaddir)
		# data_raw = readxl(loaddir * "/Eurostat aggregates/wealth_statistics_nasq_10_f_bs.xls", "Data$(k)!A12:F89");
		data_raw = XLSX.readdata(loaddir * "/Eurostat aggregates/wealth_statistics_nasq_10_f_bs.xlsx", "Data$(k)", "A12:F89")

		GEO = data_raw[1,:]
		GEO[1] = "date"

		data = DataFrame(data_raw[2:end, :], :auto)
		rename!(data, [jj => geo for (jj, geo) in enumerate(GEO)])
		data.date = Date.(["$(parse(Int, data.date[jt][1:4]))-$(parse(Int,data.date[jt][end])*3-2)" for jt in 1:length(data.date)], "yyyy-mm")

		df = DataFrame(date = data.date, assets = convert(Vector{Float64},data[!,country]))
	end

	df1 = load_nw("", loaddir)
	df2 = load_nw("2", loaddir)
	rename!(df2, "assets" => "liabilities")

	df = innerjoin(df1, df2, on = "date")
	df.net_worth = df.assets - df.liabilities
	return df
end

function make_nw_levels(; levels=false, slides = true, dark = false, template = qtemplate(;slides, dark))

	df1 = load_nw("Spain")
	dfr = load_all("../Data/")
	dfr = dfr[dfr.GEO.=="Spain",:]
	rename!(dfr, "gdp" => "GDP", "TIME" => "date")

	df2 = DataFrame(date = dfr.date, GDP = dfr.GDP)

	df = innerjoin(df1, df2, on="date")

	df = df[df.date .>= Date("2000-01-01"),:]

	ytitle = "% of GDP"
	if levels
		df.assets .*= df.GDP / 100
		df.liabilities .*= df.GDP / 100
		df.net_worth .*= df.GDP / 100
		ytitle = "Million of euros"
	end

	shapes = [
		rect([Date("2008-01-01"), Date("2010-10-01")], [Date("2009-10-01"), Date("2012-10-01")], 0, 1; fillcolor="gray", opacity=0.15, line_width=0.5, line_dash="dash", line_color="black", xref="x", yref="paper")

	]

	# plot([scatter(x=df.date, y=df[!,yvar] .* df[!,:GDP], name=yvar) for yvar in ("assets", "liabilities", "net_worth")], oth_args...)

	plot([
		bar(x=df.date, y=df.assets, name="Assets", opacity=0.5, marker_color=get(ColorSchemes.oslo, 0.5, :extrema), legendgroup=1)
		scatter(x=df.date, y=df.assets, name="Assets", marker_color=get(ColorSchemes.oslo, 0.5, :extrema), legendgroup=1, line_width = 1.5, hoverinfo="skip", showlegend=false, line_shape="spline",)
		bar(x=df.date, y=-df.liabilities, name="Liabilities", opacity=0.5, marker_color=get(ColorSchemes.lajolla, 0.6, :extrema), legendgroup=2)
		scatter(x=df.date, y=-df.liabilities, name="Assets", marker_color=get(ColorSchemes.lajolla, 0.6, :extrema), legendgroup=2, line_width = 1.5, hoverinfo="skip", showlegend=false, line_shape="spline",)
		scatter(x=df.date, y=df.net_worth, name="Net Worth", marker_color="black", line_width=3)
	], 
	Layout(template = template, barmode="overlay", title="Household sector balance sheet", yaxis_title=ytitle,
	shapes = shapes, xaxis_domain = [0.02, 1],
	legend = attr(orientation="h", x=0.05), hovermode = "x",
	)
	)
end

function make_nw(country::String = "Spain"; loaddir = "../Data/", with_annot = true, slides = true, dark = slides, template::Template=qtemplate(slides=slides, dark=dark))
    df = load_nw(country, loaddir)

    col = [
        get(ColorSchemes.davos, 0.15, :extrema)
        get(ColorSchemes.lajolla, 0.6, :extrema)
        get(ColorSchemes.cork, 0.9, :extrema)
    ]
    wths = [2, 2, 2.5]
	lsh = ["solid", "solid", "dashdot"]

    annot = []
    shapes = []

    if with_annot
        annot = [
            attr(x = Date("2009-02"), y = mean(df.assets) - 1, ax = -40, ay = 40, xanchor = "right", yanchor = "top", text = "Mean $(round(Int, mean(df.assets)))%")
            attr(x = Date("2004-01"), y = mean(df.liabilities) + 1, ax = 40, ay = -40, xanchor = "left", yanchor = "bottom", text = "Mean $(round(Int, mean(df.liabilities)))%")
        ]
        shapes = [
            hline(mean(df.assets), line = attr(width = 1, dash = "dot", color = col[1]))
            hline(mean(df.liabilities), line = attr(width = 1, dash = "dot", color = col[2]))
        ]
    end

    layout = Layout(
		template = template,
		annotations = annot, shapes = shapes,
        yaxis = attr(title = "% of GDP"),
    )

    linenames = filter(x -> x != "date", names(df))
    pNW = plot([
            scatter(; x = df.date, y = df[!, y], line_color = col[jj], line_width = wths[jj], line_dash=lsh[jj], name = uppercasefirst(replace(y, "_" => " "))) for (jj, y) in enumerate(linenames)
        ], layout)
end


function WEO_spark(loaddir = "../Data/"; own, slides = true, dark = slides, template::Template=qtemplate(slides=slides, dark=dark))

	df = CSV.read(loaddir * "/WEO/debts_spark.csv", DataFrame)

	df = df[df.year.>=2003,:]


	if own
		ytitle = "<i>% of own GDP"
	else
		ytitle = "<i>% of World GDP"
	end

	colvec = ["#0098e9", "#f97760", "#5aa800"]

	colvec[1] = "#55779A"

	sc = [
		[bar(x=df[df.ifscode.==k,:].year, y = df[df.ifscode.==k,:].debt_usd ./ df[df.ifscode.==ifelse(own, k, 1),:].ngdpd, marker_color=colvec[jk], name=first(unique(df[df.ifscode.==k,:].country))) for (jk,k) in enumerate([110, 201, 1201])]
		# scatter(x=df[df.ifscode.==1,:].year, y = df[df.ifscode.==1,:].debt_usd ./ df[df.ifscode.==1,:].ngdpd)
		]

	layout = Layout(
		template=template,
		barmode=ifelse(own, "", "stack"), yaxis_title=ytitle, title="Government debts globally (WEO Oct 2020)",
		)

	plot(sc, layout)
end

function make_FLP(loaddir = "../Data/"; slides = true, dark = slides, template::Template=qtemplate(slides=slides, dark=dark))
	vnraw = XLSX.readdata(loaddir * "/Eurostat aggregates/ei_bsin_q_r2.xlsx", "Sheet 1!B9:G9")
	dtraw = XLSX.readdata(loaddir * "/Eurostat aggregates/ei_bsin_q_r2.xlsx", "Sheet 1!A71:G179")

	varnames = ["date"]
	for vn in vnraw
		push!(varnames, vn[35:end])
	end

	df = DataFrame(dtraw, varnames)

	df.date = Date.(["$(parse(Int, df.date[jt][1:4]))-$(parse(Int,df.date[jt][end])*3-2)" for jt in 1:length(df.date)], "yyyy-mm")

	df = df[df.date .>= Date("2002-01-01"), :]
	df = df[df.date .< Date("2020-01-01"), :]

	scats = [
		scatter(x=df.date, y=df[!,key], name = key) for key in names(df) if key != "date"
	]
	layout = Layout(
		template = template,
		legend = attr(x=0.5, xanchor="center", xref="paper"),
		yaxis = attr(title="<i>%")
	)

	plot(scats, layout)
end