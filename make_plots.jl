using PlotlyJS, JLD, LaTeXStrings, Rsvg, Interpolations

p_names = ["conv" "hh" "hh_z" "hh_b" "hh_def" "statefuncs" "LoMs" "nontradables_B" "nontradables_z" "labordemand" "simul" "simul_full" "objfunc" "reactions" "aggcons" "varcons" "outconv"]

for name in p_names
	print("$name")
	try
		if name == "varcons"
			p = load(pwd() * "/../HPC_Output/p_" * name * ".jld", "p2")
		else
			p = load(pwd() * "/../HPC_Output/p_" * name * ".jld", "p")
		end
		savefig(p, pwd() * "/../Graphs/" * name * ".png")
		print(": ✓\n")
	catch
		print("\n")
	end
end
