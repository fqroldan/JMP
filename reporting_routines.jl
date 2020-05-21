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
		:β => ["Discount rate of HHs			", "\$ 1 / \\beta - 1 \$		", (sd.pars[:β]^(-4)-1), true],
		:γ => ["Risk aversion					", "\$ \\gamma \$				", sd.pars[:γ], false],
		:τ => ["Progressivity of tax schedule 	", "\$ \\tau \$					", sd.pars[:τ], true],
		:w => ["Wage minimum 					", "\$ \\bar{w} \$				", sd.pars[:wbar], false],
		:z => ["TFP process						", "\$ \\rho_z, \\sigma_z \$	", [sd.pars[:ρz], sd.pars[:σz]], false],
		:ξ0=> ["Mean risk premium 				", "\$\\bar{\\xi}\$				", sd.pars[:meanξ], true],
		:ξ => ["Risk premium AR(1)				", "\$\\rho_\\xi, \\sigma_\\xi\$", [sd.pars[:ρξ], sd.pars[:σξ]], false],
		)

	for key in [:β, :γ, :τ, :w, :z, :ξ0, :ξ]
		table *= "\n" * td[key][1] * "&" * td[key][2] * "&"
		if key == :z || key == :ξ
			table *= "($(@sprintf("%0.3g",td[key][3][1])), $(@sprintf("%0.3g",(td[key][3][2]))))"
		elseif td[key][4]
			table *= "$(@sprintf("%0.3g",(100*td[key][3])))" * "\\%"
			if key == :β
				table *= "ann."
			end
		else
			table *= "$(@sprintf("%0.3g",td[key][3]))"
		end

		table *= " & Moments in Table \\ref{tab:calibration} \\\\"
	end
	table *= "\n		\\bottomrule"
	# write(pwd()*"/../../params_table.txt", last_table)
	# nothing
	return table
end

function make_calib_table(v_m)
	table = ""

	rownames = ["AR(1) coef \$\\log(Y_t)\$"; "Std coef \$\\log(Y_t)\$"; "AR(1) coef \$\\log(C_t)\$"; "Std coef \$\\log(C_t)\$"; "AR(1) coef spread"; "Std coef spread	"; "Avg Debt-to-GDP	"; "Std Debt-to-GDP	"; "Avg unemployment"; "Std unemployment"; "Median dom holdings"; "Avg wealth-to-GDP"]

	data_stats = [ 0.96580506; 0.01294576; 0.96172496; 0.01663608; 0.96656486; 0.10252351; 64.57638889; 23.48323041; 15.94722222; 6.08732167; 56.49; 94.48]

	list_perc = ones(size(data_stats))
	list_perc[1:6] = zeros(6)

	for jj in 1:length(data_stats)
		list_perc[jj] == 1 ? perc = "\\%" : perc = ""
		table *= "\n		" * rownames[jj] * "	& 	$(@sprintf("%0.3g",v_m[jj]))" * perc * "	& 	$(@sprintf("%0.3g",data_stats[jj]))" * perc * "\\\\"
	end
	table *= "\n		\\bottomrule"

	return table
end

function make_calib_table_comp(v_m, v_m_nodef, v_m_noΔ=[], v_m_nob=[])
	table = "Moment 			& Benchmark "
	if length(v_m_noΔ) > 0
		table *=  " 	& \$\\Delta=0\$"
	end
	if length(v_m_nob) > 0
		table *=  " 	& No dom. holdings \$"
	end
	table *= "	& No default "
	table *= "\\\\ \\midrule"

	rownames = ["AR(1) coef \$\\log(Y_t)\$"; "Std coef \$\\log(Y_t)\$"; "AR(1) coef \$\\log(C_t)\$"; "Std coef \$\\log(C_t)\$"; "AR(1) coef spread"; "Std coef spread	"; "Avg Debt-to-GDP	"; "Std Debt-to-GDP	"; "Avg unemployment"; "Std unemployment"; "Median dom holdings"; "Avg wealth-to-GDP"; "Default frequency"]

	data_stats = [ 0.96580506; 0.01294576; 0.96172496; 0.01663608; 0.96656486; 0.10252351; 64.57638889; 23.48323041; 15.94722222; 6.08732167; 56.49; 94.48]

	list_perc = ones(length(v_m))
	list_perc[1:6] = zeros(6)

	for jj in 1:length(v_m)
		list_perc[jj] == 1 ? perc = "\\%" : perc = ""
		table *= "\n" * rownames[jj] * "	& $(show_float(v_m[jj]))" * perc
		if length(v_m_noΔ) > 0
			table *=  " 	& $(show_float(v_m_noΔ[jj]))" * perc
		end
		if length(v_m_nob) > 0
			table *=  " 	& $(show_float(v_m_nob[jj]))" * perc
		end
		table *= "	& $(show_float(v_m_nodef[jj]))" * perc

		table *= "\\\\"
	end
	table *= "\n		\\bottomrule"

	return table
end