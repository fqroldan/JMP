using QuantEcon, BasisMatrices, Interpolations, Optim, NLopt, MINPACK, Plots, LaTeXStrings, Distributions
pyplot()

# Load codes
include("../../Julia/printtime.jl")
# include("hh_col.jl")
include("hh_reiter.jl")

# Run
h = Hank();
showtext = "\nA Theory of Sovereign Risk\n"
write(pwd()*"/../output.txt", showtext)
print(showtext)

vfi!(h, verbose = true)
# @time A_sample, B_sample, L_sample, Z_sample, R_sample, q_sample, qg_sample, π_sample = simul(h, makefile = true)

Void