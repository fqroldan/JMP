using PlotlyJS, JLD, LaTeXStrings, Rsvg, Interpolations
# using Plots, JLD, LaTeXStrings, Rsvg
# plotlyjs()

# include("type_def.jl")
# h = load(pwd() * "/../HPC_Output/hank.jld", "h")

# include("plotting_routines.jl")
# plot_state_funcs(h, remote=false)
# plot_LoM(h, remote=false)
# plot_labor_demand(h, remote=false)
# plot_hh_policies(h, remote=false)

p_names = ["conv" "hh" "hh_def" "statefuncs" "LoMs" "labordemand" "simul" "objfunc" "reactions"]

for name in p_names
	print("$name")
	try
		p = load(pwd() * "/../HPC_Output/p_" * name * ".jld", "p")
		savefig(p, pwd() * "/../Graphs/" * name * ".png")
		print(": ✓\n")
	catch
		print("\n")
	end
end
