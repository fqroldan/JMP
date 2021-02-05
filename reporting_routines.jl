function print_save(s::String, dir::String = pwd()*"/../Output/")
	print(s)
	output = read(dir * "output.txt", String)
	write(dir * "output.txt", output * s)

	nothing
end

show_float(x; digits=10) = @sprintf("%.3g", round(x, digits=digits))

function time_print(tfloat::Float64)
	t = floor(Int, tfloat)
	if t < 1
		t_print = "no time"
	elseif t < 60
		t_print = "$(t) second"
		t == 1 ? nothing : t_print = t_print * "s"
	elseif t < 3600
		t_print = "$(floor(Int,t/60)) minute"
		floor(Int,t/60) == 1 ? nothing : t_print = t_print * "s"
		if t % 60 != 0
			t_print = t_print * " and $(t%60) second"
			t % 60 == 1 ? nothing : t_print = t_print * "s"
		end
	else
		t_print = "$(floor(Int,t/3600)) hour"
		floor(Int, t/3600) == 1 ? nothing : t_print = t_print * "s"
		t = t % 3600
		t_print = t_print * " and $(floor(Int,t/60)) minute"
		floor(Int,t/60) == 1 ? nothing : t_print = t_print * "s"
		if t % 60 != 0
			t_print = t_print * " and $(t%60) second"
			t%60 == 1 ? nothing : t_print = t_print * "s"
		end
	end
	return t_print
end

function make_params_table(sd::SOEdef)
	table = ""

	td = Dict(
		:β => ["Discount rate of HHs", 			"\$ 1 / \\beta - 1 \$", 		(sd.pars[:β]^(-4)-1), true],
		:γ => ["Risk aversion", 				"\$ \\gamma \$", 				sd.pars[:γ], false],
		:τ => ["Progressivity of tax schedule", "\$ \\tau \$", 					sd.pars[:τ], true],
		:w => ["Wage minimum", 					"\$ \\bar{w} \$", 				sd.pars[:wbar], false],
		:z => ["TFP process", 					"\$ \\rho_z, \\sigma_z \$", 	[sd.pars[:ρz], sd.pars[:σz]], false],
		:ξ0=> ["Mean risk premium", 			"\$\\bar{\\xi}\$", 				sd.pars[:meanξ], true],
		:ξ => ["Risk premium AR(1)", 			"\$\\rho_\\xi, \\sigma_\\xi\$", [sd.pars[:ρξ], sd.pars[:σξ]], false],
		:μ_gov => ["Mean utility cost of default", "\$\\mu_g\$", sd.pars[:μ_gov], false],
		# :σ_gov => ["Std utility cost of default", "\$\\sigma_g\$", sd.pars[:σ_gov], false]
		)

	for key in [:β, :γ, :τ, :w, :z, :ξ0, :ξ, :μ_gov]
		table *= "\n" * rpad(td[key][1],30," ") * "&" * rpad(td[key][2],30," ") * "&"
		if key == :z || key == :ξ
			table *= "($(@sprintf("%0.4g",td[key][3][1])), $(@sprintf("%0.4g",(td[key][3][2]))))"
		elseif td[key][4]
			table *= "$(@sprintf("%0.4g",(100*td[key][3])))" * "\\%"
			if key == :β
				table *= " ann."
			end
		else
			table *= "$(@sprintf("%0.4g",td[key][3]))"
		end

		table *= " & Moments in Table \\ref{tab:calibration} \\\\"
	end
	table *= "\n		\\bottomrule"
	# write(pwd()*"/../../params_table.txt", last_table)
	# nothing
	return table
end

function make_calib_table(v_m)
	table = "\\begin{tabular*}{.5\\textwidth}{@{\\extracolsep{20pt}}lcc@{}} \\toprule \n"

	colnames = ["\\textbf{Target}", "\\textbf{Model}", "\\textbf{Data}"]

	rownames = ["AR(1) autocorr.~coef \$\\log(Y_t)\$"; "AR(1) std coef \$\\log(Y_t)\$"; "AR(1) autocorr.~coef \$\\log(C_t)\$"; "AR(1) std coef \$\\log(C_t)\$"; "AR(1) autocorr.~coef spread"; "AR(1) std coef spread"; "Avg Debt-to-GDP	"; "Std Debt-to-GDP	"; "Avg unemployment"; "Std unemployment"; "Median dom holdings"; "Avg wealth-to-GDP"; "Avg wealth Gini"]

	# data_stats = [ 0.96580506; 0.01294576; 0.96172496; 0.01663608; 0.96656486; 0.10252351; 64.57638889; 23.48323041; 15.94722222; 6.08732167; 56.49; 94.48]

	data_stats = load_SPA_targets()

	list_perc = ones(size(data_stats))
	list_perc[[1,3,5,6]] .= 0

	vs = [@sprintf("%0.3g",vv)*ifelse(list_perc[jv]==1,"\\%", "") for (jv, vv) in enumerate(v_m)]
	ds = [@sprintf("%0.3g",dv)*ifelse(list_perc[jv]==1,"\\%", "") for (jv, dv) in enumerate(data_stats)]

	pad_n = max(length(colnames[1]), maximum(length.(rownames))) + 1
	pad_v = max(length(colnames[2]), maximum(length.( vs ))) + 3
	pad_d = max(length(colnames[3]), maximum(length.( ds ))) + 3

	table *= rpad(colnames[1],pad_n, " ") *  "& " * rpad(colnames[2], pad_v, " ") * "& " *	rpad(colnames[3], pad_d, " ") * "\\\\ \\midrule \n"
	for jj in 1:length(data_stats)
		table *= rpad(rownames[jj],pad_n," ") * "& " * rpad(vs[jj], pad_v, " ") * "& "* rpad(ds[jj], pad_d, " ") * "\\\\\n"
	end
	table *= "\\bottomrule\n"

	# table *= " \\multicolumn{3}{@{}p{.5\\textwidth}@{}}{\\footnotesize All data from Eurostat 2000Q1:2017Q4, except private consumption from OECD 2000Q1:2017Q4, domestic holdings from Banco de España, 2004Q1:2017Q4}\n"

	table *= "\\end{tabular*}"
	return table
end

function make_calib_table_comp(v_m, v_m_nodef, v_m_noΔ=[], v_m_nob=[])
	k = sum(length(v) > 0 for v in [ v_m, v_m_nodef, v_m_noΔ, v_m_nob ])
	table = "\\begin{tabular*}{.65\\textwidth}{@{\\extracolsep{\\fill}}l*{$k}c@{}} \\toprule \n"

	colnames = ["\\textbf{Moment}", "\\textbf{Benchmark}"]
	rownames = ["AR(1) autocorr.~coef \$\\log(Y_t)\$"; "AR(1) std coef \$\\log(Y_t)\$"; "AR(1) autocorr.~coef \$\\log(C_t)\$"; "AR(1) std coef \$\\log(C_t)\$"; "AR(1) autocorr.~coef spread"; "AR(1) std coef spread"; "Avg Debt-to-GDP	"; "Std Debt-to-GDP	"; "Avg unemployment"; "Std unemployment"; "Median dom holdings"; "Avg wealth-to-GDP"; "Avg wealth Gini"; "Default frequency"; "Welfare in repayment"]
	pad_n = maximum(length.(rownames)) + 1

	if length(v_m_noΔ) > 0
		push!(colnames, "\\textbf{\$\\Delta=0\$}")
	end
	if length(v_m_nob) > 0
		push!(colnames, "\\textbf{No dom.~holdings}")
	end
	push!(colnames, "\\textbf{No default}")

	list_perc = ones(length(v_m))
	list_perc[[1,3,5,6]] .= 0.0
	list_perc[end] = 0.0
	vs = [@sprintf("%0.3g",vv)*ifelse(list_perc[jv]==1,"\\%", "") for (jv, vv) in enumerate(v_m)]
	# pad_v = max(length(colnames[2]), maximum(length.( vs ))) + 3
	pad_v = maximum(length.( vs )) + 3

	pad_n = max(length(colnames[1]), maximum(length.(rownames))) + 1

	table *= rpad(colnames[1],pad_n, " ") * "& " * rpad(colnames[2], pad_v, " ")

	jj = 3
	if length(v_m_noΔ) > 0
		v_noΔ = [@sprintf("%0.3g",round(vv,digits=8))*ifelse(list_perc[jv]==1,"\\%", "") for (jv, vv) in enumerate(v_m_noΔ)]
		# pad_Δ = max(length(colnames[jj]), maximum(length.( v_noΔ ))) + 3
		pad_Δ = maximum(length.( v_noΔ )) + 3
		table *= "& " * rpad(colnames[jj], pad_Δ, " ")
		jj += 1
	end
	if length(v_m_nob) > 0
		v_nob = [@sprintf("%0.3g",round(vv,digits=8))*ifelse(list_perc[jv]==1,"\\%", "") for (jv, vv) in enumerate(v_m_nob)]
		# pad_b = max(length(colnames[jj]), maximum(length.( v_nob ))) + 3
		pad_b = maximum(length.( v_nob )) + 3
		table *= "& " * rpad(colnames[jj], pad_b, " ")
		jj += 1
	end
	v_nodef = [@sprintf("%0.3g",round(vv,digits=8))*ifelse(list_perc[jv]==1,"\\%", "") for (jv, vv) in enumerate(v_m_nodef)]
	# pad_def = max(length(colnames[jj]), maximum(length.( v_nodef ))) + 3
	pad_def = maximum(length.( v_nodef )) + 3
	table *= "& " * rpad(colnames[jj], pad_def, " ")
	
	table *= "\\\\ \\midrule \n"


	for jj in 1:length(rownames)
		table *= rpad(rownames[jj],pad_n," ") * "& " * rpad(vs[jj], pad_v, " ")
		if length(v_m_noΔ) > 0
			table *= "& " * rpad(v_noΔ[jj], pad_Δ, " ")
		end
		if length(v_m_nob) > 0
			table *= "& " * rpad(v_nob[jj], pad_b, " ")
		end
		table *= "& " * rpad(v_nodef[jj], pad_def, " ")
		table *= "\\\\\n"

	end

	table *= "\\bottomrule \n\\end{tabular*}"
	return table
end


function make_RH_table(v_m, v_rep)
	k = 2
	table = "\\begin{tabular*}{.8\\textwidth}{@{\\extracolsep{\\fill}}l*{$k}c@{}} \\toprule \n"
	colnames = ["\\textbf{Moment}", "\\textbf{Benchmark}", "\\textbf{Rep.~Agent}"]

	rownames = ["AR(1) coef \$\\log(Y_t)\$"; "Std coef \$\\log(Y_t)\$"; "AR(1) coef \$\\log(C_t)\$"; "Std coef \$\\log(C_t)\$"; "AR(1) coef spread"; "Std coef spread"; "Avg Debt-to-GDP	"; "Std Debt-to-GDP	"; "Avg unemployment"; "Std unemployment"; "Median dom holdings"; "Avg wealth-to-GDP"; "Avg Gini"; "Default frequency"; "Value in repayment"]

	pad_n = maximum(length.(rownames)) + 1

	list_perc = ones(length(v_m))
	list_perc[1:6] .= 0.0
	list_perc[end] .= 0.0
	vs = [@sprintf("%0.3g",vv)*ifelse(list_perc[jv]==1,"\\%", "") for (jv, vv) in enumerate(v_m)]
	# pad_v = max(length(colnames[2]), maximum(length.( vs ))) + 3
	pad_v = maximum(length.( vs )) + 3

	pad_n = max(length(colnames[1]), maximum(length.(rownames))) + 1

	table *= rpad(colnames[1],pad_n, " ") * "& " * rpad(colnames[2], pad_v, " ")

	jj = 3
	vreps = [@sprintf("%0.3g",round(vv,digits=8))*ifelse(list_perc[jv]==1,"\\%", "") for (jv, vv) in enumerate(v_rep)]
	# pad_def = max(length(colnames[jj]), maximum(length.( vreps ))) + 3
	pad_def = maximum(length.( vreps )) + 3
	table *= "& " * rpad(colnames[jj], pad_def, " ")
	
	table *= "\\\\ \\midrule \n"


	for jj in 1:length(rownames)
		table *= rpad(rownames[jj],pad_n," ") * "& " * rpad(vs[jj], pad_v, " ")
		table *= "& " * rpad(vreps[jj], pad_def, " ")
		table *= "\\\\\n"

	end

	table *= "\\bottomrule \n\\end{tabular*}"
	return table
end


function make_welfare_table(p_bench, p_nodef)

	table = "\\begin{tabular*}{.55\\textwidth}{@{\\extracolsep{\\fill}}l*{3}c@{}} \\toprule \n"
	colnames = ["\\textbf{Moment}", "\\textbf{Benchmark}", "\\textbf{No default}", "\\textbf{Gains}"]
	rownames = ["\$\\,p10\$", "\$\\,p25\$", "\$\\,p50\$", "\$\\,p75\$", "\$\\,p90\$", "Average"]
	varnames = [:Wr10, :Wr25, :Wr50, :Wr75, :Wr90, :Wr]

	pad_n = maximum(length.(rownames))
	pad_n = max(pad_n, length(colnames[1])) + 2

	v_bench = [ mean([mean(series(p, key)) for p in p_bench]) for key in varnames ]
	pad_bench = maximum(length.(v_bench))
	pad_bench = max(pad_bench, length(colnames[2])) + 2
	v_nodef = [ mean([mean(series(p, key)) for p in p_nodef]) for key in varnames ]
	pad_nodef = maximum(length.(v_nodef))
	pad_nodef = max(pad_bench, length(colnames[3])) + 2
	v_rel   = [ rel_var(p_bench, p_nodef, key) for key in varnames ]
	pad_rel = maximum(length.(v_rel))
	pad_rel = max(pad_bench, length(colnames[4])) + 2


	pads = [pad_n, pad_bench, pad_nodef, pad_rel]
	seps = ["& ", "& ", "& ", "\\\\ "]
	for jj in eachindex(colnames)
		table *= rpad(colnames[jj], pads[jj], " ") * seps[jj]
	end

	table *= "\\midrule\n"

	for (jj, key) in enumerate(rownames)
		table *= rpad(rownames[jj], pads[1], " ") * seps[1]
		table *= rpad(@sprintf("%0.3g", v_bench[jj]), pads[2], " ") * seps[2]
		table *= rpad(@sprintf("%0.3g", v_nodef[jj]), pads[3], " ") * seps[3]
		table *= rpad(@sprintf("%0.3g", 100*v_rel[jj])*"\\%", pads[4], " ") * seps[4] * "\n"
	end

	table *= "\\bottomrule \n\\end{tabular*}"
	return table
end

function get_def_freq(pvec)
	DF = 0.0
	for pp in pvec
		ζ_vec = series(pp,:ζ)
		Ndefs = sum([ (ζ_vec[jj] == 0) & (ζ_vec[jj-1] == 1) for jj in 2:length(ζ_vec)])
		years = floor(Int64,periods(pp)*0.25)
		def_freq = Ndefs / years

		DF += def_freq / length(pvec)
	end
	DF
end

function save_all_tables(savedir = "../HPC_output/current_best/")

	sd = FileIO.load(savedir * "SOEdef.jld2", "sd")
	p_bench, W_bench = FileIO.load(savedir * "SOEdef.jld2", "pp", "Wr")
	# p_nodef, W_nodef = File.load(savedir * "SOEdef_nodef.jld", "pp", "Wr")
	# p_noΔ, W_noΔ = load(savedir * "SOEdef_nodelta.jld", "pp", "Wr")
	# p_nob, W_nob = load(savedir * "SOEdef_nob.jld", "pp", "Wr")

	write(savedir * "params_table.txt", make_params_table(sd))

	v_bench = simul_stats(p_bench)
	freq_bench = get_def_freq(p_bench)
	v_nodef = simul_stats(p_nodef)
	freq_nodef = get_def_freq(p_nodef)
	# v_noΔ = simul_stats(p_noΔ)
	# freq_noΔ = get_def_freq(p_noΔ)
	# v_nob = simul_stats(p_nob)
	# freq_nob = get_def_freq(p_nob)

	write(savedir * "calib_table.txt", make_calib_table(v_bench))

	write(savedir * "calib_table_comp.txt", 
		make_calib_table_comp(
			[v_bench; 100*freq_bench; W_bench],
			[v_nodef; 100*freq_nodef; W_nodef],
			# [v_noΔ; 100*freq_noΔ; W_noΔ],
			# [v_nob; 100*freq_nob; W_nob],
			))

	write(savedir * "welfare_table.txt", make_welfare_table(p_bench, p_nodef))
end

function save_all_plots(savedir = "../HPC_output/current_best/")

	sd_bench, p_bench, W_bench = FileIO.load(savedir * "SOEdef.jld2", "sd", "pp", "Wr")
	sd_nodef, p_nodef, W_nodef = FileIO.load(savedir * "SOEdef_nodef.jld2", "sd", "pp", "Wr")

	p1 = make_panels(sd, "Earnings_Default", style=paper, leg=false)
	savefig(p1, savedir * "wages_paper.pdf", width=floor(Int, p1.plot.layout[:width]), height=floor(Int,p1.plot.layout[:height]))

	p1 = make_twisted(sd, style=paper, leg=true)	
	savefig(p1, savedir * "twistedprobs_paper.pdf", width=floor(Int, p1.plot.layout[:width]), height=floor(Int,p1.plot.layout[:height]))

	p1 = make_panels(sd, "Wr-Wd", style=paper, leg=false)
	savefig(p1, savedir * "WrWd_paper.pdf", width=floor(Int, p1.plot.layout[:width]), height=floor(Int,p1.plot.layout[:height]))

	p1 = make_panels(sd, "T", style=paper, leg=false)
	savefig(p1, savedir * "transfers_paper.pdf", width=floor(Int, p1.plot.layout[:width]), height=floor(Int,p1.plot.layout[:height]))

	p1 = make_debtprice(sd, style=paper, leg=false)
	savefig(p1, savedir * "debtprice_paper.pdf", width=floor(Int, p1.plot.layout[:width]), height=floor(Int,p1.plot.layout[:height]))

	p1 = make_unemp(sd, style=paper, leg=false)
	savefig(p1, savedir * "unemp_paper.pdf", width=floor(Int, p1.plot.layout[:width]), height=floor(Int,p1.plot.layout[:height]))

	p1 = panels_crises(p_bench, 400, :spread, k=1, thres_back = 350, k_back = 11, style=paper)
	savefig(p1, savedir * "highspreads_paper.pdf", width=floor(Int, p1.plot.layout[:width]), height=floor(Int,p1.plot.layout[:height]))

	p1 = panels_defaults(p_bench, style=paper, indiv=false)
	savefig(p1, savedir * "defaults_paper.pdf", width=floor(Int, p1.plot.layout[:width]), height=floor(Int,p1.plot.layout[:height]))

	p1 = panels_comp(p_bench, p_nodef, 400, :spread, thres_back = 350, k=1, k_back=11, style=paper, yh=0.55, relative=true)
	savefig(p1, savedir * "crises_comp_paper.pdf", width=floor(Int, p1.plot.layout[:width]), height=floor(Int,p1.plot.layout[:height]))

	p1 = distribution_comp(p_bench, p_nodef, 400, :spread, thres_back = 350, k=1, k_back=11, style=paper, yh=0.5, transf=false, response = "W")
	savefig(p1, savedir * "distrib_crises_paper.pdf", width=floor(Int, p1.plot.layout[:width]), height=floor(Int,p1.plot.layout[:height]))

	nothing
end
