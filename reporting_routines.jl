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
		)

	for key in [:β, :γ, :τ, :w, :z, :ξ0, :ξ]
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

	table = "\\begin{tabular}{@{\\extracolsep{20pt}}lcc@{}} \\toprule \n"

	colnames = ["\\textbf{Target}", "\\textbf{Model}", "\\textbf{Data}"]

	rownames = ["AR(1) coef \$\\log(Y_t)\$"; "Std coef \$\\log(Y_t)\$"; "AR(1) coef \$\\log(C_t)\$"; "Std coef \$\\log(C_t)\$"; "AR(1) coef spread"; "Std coef spread	"; "Avg Debt-to-GDP	"; "Std Debt-to-GDP	"; "Avg unemployment"; "Std unemployment"; "Median dom holdings"; "Avg wealth-to-GDP"]

	data_stats = [ 0.96580506; 0.01294576; 0.96172496; 0.01663608; 0.96656486; 0.10252351; 64.57638889; 23.48323041; 15.94722222; 6.08732167; 56.49; 94.48]

	list_perc = ones(size(data_stats))
	list_perc[1:6] = zeros(6)

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

	table *= "\\end{tabular}"
	return table
end

function make_calib_table_comp(v_m, v_m_nodef, v_m_noΔ=[], v_m_nob=[])
	k = sum(length(v) > 0 for v in [ v_m, v_m_nodef, v_m_noΔ, v_m_nob ])
	table = "\\begin{tabular*}{.85\\textwidth}{@{\\extracolsep{\\fill}}l*{$k}c@{}} \\toprule \n"

	colnames = ["\\textbf{Moment}", "\\textbf{Benchmark}"]
	rownames = ["AR(1) coef \$\\log(Y_t)\$"; "Std coef \$\\log(Y_t)\$"; "AR(1) coef \$\\log(C_t)\$"; "Std coef \$\\log(C_t)\$"; "AR(1) coef spread"; "Std coef spread	"; "Avg Debt-to-GDP	"; "Std Debt-to-GDP	"; "Avg unemployment"; "Std unemployment"; "Median dom holdings"; "Avg wealth-to-GDP"; "Avg wealth Gini"; "Default frequency"; "Welfare in repayment"]
	pad_n = maximum(length.(rownames)) + 1

	if length(v_m_noΔ) > 0
		push!(colnames, "\\textbf{\$\\Delta=0\$}")
	end
	if length(v_m_nob) > 0
		push!(colnames, "\\textbf{No dom.~holdings}")
	end
	push!(colnames, "\\textbf{No default}")

	list_perc = ones(length(v_m))
	list_perc[1:6] .= 0.0
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

	rownames = ["AR(1) coef \$\\log(Y_t)\$"; "Std coef \$\\log(Y_t)\$"; "AR(1) coef \$\\log(C_t)\$"; "Std coef \$\\log(C_t)\$"; "AR(1) coef spread"; "Std coef spread	"; "Avg Debt-to-GDP	"; "Std Debt-to-GDP	"; "Avg unemployment"; "Std unemployment"; "Median dom holdings"; "Avg wealth-to-GDP"; "Avg Gini"; "Default frequency"; "Value in repayment"]

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

	# [v_m; 100*def_freq; Gini; Wr], [v_rep, freq_rep, Gini_rep, W_rep]