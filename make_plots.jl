using Plots, JLD, LaTeXStrings
gr()

include("type_def.jl")
include("plotting_routines.jl")


h = load(pwd() * "/../HPC_Output/hank.jld", "h")

plot_labor_demand(h)
plot_hh_policies(h)

