using Printf
function print_save(s::String, dir::String = pwd()*"/../Output/")
	print(s)
	output = read(dir * "output.txt", String)
	write(dir * "output.txt", output * s)

	nothing
end

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

function make_params_table(params)
	table = ""

	corr = ones(size(params))
	corr[1] = 100
	corr[3] = 100
	corr[7] = 100

	rownames = ["Discount rate of HHs			"; "Risk aversion					"; "Progressivity of tax schedule	"; "Wage minimum					"; "TFP process						"; "Mean risk premium				"; "Risk premium AR(1)				"]
	descrips = ["\$ 1 / \\beta - 1 \$"; "\$ \\gamma \$		"; "\$ \\tau \$			"; "\$ \\bar{w} \$		"; "\$ \\rho_z, \\sigma_z \$"; "\$\\bar{\\xi}\$"; "\$\\rho_\\xi, \\sigma_\\xi\$"]
	single = zeros(size(rownames))

	jpar = 1
	for jj in 1:length(rownames)
		if jj != 5 && jj != 7
			fill_param = "$(@sprintf("%0.3g",params[jpar]*corr[jpar]))"
			if jj == 1
				fill_param *= "\\% ann."
			end
			if jj == 3 || jj == 6
				fill_param *= "\\%"
			end
			jpar += 1
		else
			fill_param = "($(@sprintf("%0.3g",params[jpar])), $(@sprintf("%0.3g",params[jpar+1])))"
			jpar += 2
		end
		table *= "\n		" * rownames[jj] * " & " * descrips[jj] * "	& " * fill_param * " & Moments in Table \\ref{tab:calibration} \\\\"
	end
	table *= "\n		\\bottomrule"
	# write(pwd()*"/../../params_table.txt", last_table)
	# nothing
	print(table)
end

function make_calib_table(v_m)
	table = ""

	rownames = ["AR(1) coef \$\\log(Y_t)\$"; "Std coef \$\\log(Y_t)\$"; "AR(1) coef \$\\log(C_t)\$"; "Std coef \$\\log(C_t)\$"; "AR(1) coef spread	"; "Std coef spread		"; "Avg Debt-to-GDP		"; "Std Debt-to-GDP		"; "Avg unemployment	"; "Std unemployment	"; "Median dom holdings	"; "Avg wealth-to-GDP	"]

	data_stats = [ 0.96580506; 0.01294576; 0.96172496; 0.01663608; 0.96656486; 0.10252351; 64.57638889; 23.48323041; 15.94722222; 6.08732167; 56.49; 94.48]

	list_perc = ones(size(data_stats))
	list_perc[1:6] = zeros(6)

	for jj in 1:length(data_stats)
		list_perc[jj] == 1 ? perc = "\\%" : perc = ""
		table *= "\n		" * rownames[jj] * "	& 	$(@sprintf("%0.3g",v_m[jj]))" * perc * "	& 	$(@sprintf("%0.3g",data_stats[jj]))" * perc * "\\\\"
	end
	table *= "\n		\\bottomrule"

	print(table)
end

function make_calib_table_comp(v_m, v_m_nodef, v_m_nodelta=[])
	table = ""

	rownames = ["AR(1) coef \$\\log(Y_t)\$"; "Std coef \$\\log(Y_t)\$"; "AR(1) coef \$\\log(C_t)\$"; "Std coef \$\\log(C_t)\$"; "AR(1) coef spread	"; "Std coef spread		"; "Avg Debt-to-GDP		"; "Std Debt-to-GDP		"; "Avg unemployment	"; "Std unemployment	"; "Median dom holdings	"; "Avg wealth-to-GDP	"]

	data_stats = [ 0.96580506; 0.01294576; 0.96172496; 0.01663608; 0.96656486; 0.10252351; 64.57638889; 23.48323041; 15.94722222; 6.08732167; 56.49; 94.48]

	list_perc = ones(size(data_stats))
	list_perc[1:6] = zeros(6)

	for jj in 1:length(data_stats)
		list_perc[jj] == 1 ? perc = "\\%" : perc = ""
		table *= "\n		" * rownames[jj] * "	& 	$(@sprintf("%0.3g",v_m[jj]))" * perc
		if length(v_m_nodelta) > 0
			table *=  "& 	$(@sprintf("%0.3g",v_m_nodelta[jj]))" * perc
		end
		table *= "	& 	$(@sprintf("%0.3g",v_m_nodef[jj]))" * perc

		table *= "\\\\"
	end
	table *= "\n		\\bottomrule"

	print(table)
end