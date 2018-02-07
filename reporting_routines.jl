function print_save(s::String)
	print(s)
	output = readstring(pwd()*"/../../output.txt")
	write(pwd()*"/../../output.txt", output * s)

	Void
end

function time_print(tfloat::Float64)
	t = floor(Int, tfloat)
	if t < 60
		t_print = "$(t) second"
		t == 1 ? Void : t_print = t_print * "s"
	elseif t < 3600
		t_print = "$(floor(Int,t/60)) minute"
		floor(Int,t/60) == 1 ? Void : t_print = t_print * "s"
		if t % 60 != 0
			t_print = t_print * " and $(t%60) second"
			t % 60 == 1 ? Void : t_print = t_print * "s"
		end
	else
		t_print = "$(floor(Int,t/3600)) hour"
		floor(Int, t/3600) == 1? Void: t_print = t_print * "s"
		t = t % 3600
		t_print = t_print * " and $(floor(Int,t/60)) minute"
		floor(Int,t/60) == 1 ? Void : t_print = t_print * "s"
		if t % 60 != 0
			t_print = t_print * " and $(t%60) second"
			t%60 == 1 ? Void : t_print = t_print * "s"
		end
	end
	return t_print
end
