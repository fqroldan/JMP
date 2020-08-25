using CSV, DataFrames, Dates, GLM, PlotlyJS, ORCA, FixedEffectModels, RegressionTables

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

	data.cons = data.consnom ./ data.cpi
	data.lcons = log.(data.cons)
	data.lgdp = log.(data.gdp)

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

	df = data
	
	# fs = reg(df, @formula(spread ~ b0 + fe(TIME)), save=true)
	# df.spr_hat = coef(fs)[1] * df.b0 + fe(fs).fe_TIME

	regY1 = reg(df, @formula(lgdp ~ spread + fe(GEO) + fe(TIME)))
	regC1 = reg(df, @formula(lcons ~ spread + fe(GEO) + fe(TIME)))

	regY2 = reg(df, @formula(lgdp ~ spread + debt + fe(GEO) + fe(TIME)))
	regC2 = reg(df, @formula(lcons ~ spread + debt + fe(GEO) + fe(TIME)))

	res1 = reg(df, @formula(lcons ~ spread + lgdp + fe(GEO) + fe(TIME)))
	res2 = reg(df, @formula(lcons ~ spread + debt + lgdp + fe(GEO) + fe(TIME)))

	regtable(regY1, regY2, regC1, regC2, res1, res2, renderSettings = latexOutput())
end