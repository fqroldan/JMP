# Initialize output file
write(pwd()*"/../../output.txt", "")

using QuantEcon, BasisMatrices, Interpolations, Optim, NLopt, MINPACK, LaTeXStrings, Distributions, JLD

# Load codes
@everywhere include("reporting_routines.jl")
@everywhere include("type_def.jl")
@everywhere include("reiter.jl")
@everywhere include("comp_eqm.jl")
include("plotting_routines.jl")

print_save("\nA Theory of Sovereign Risk\n")

location = "remote"
if pwd() == "/home/q/Dropbox/NYU/AToSR/Codes"
	location = "local"
	using Rsvg
end

print_save("\nStarting $(location) run on $(nprocs()) cores at "*Dates.format(now(), "HH:MM"))

# Initialize type
h = Hank();
try
	# h2 = load("hank.jld", "h")
	# h.ϕa = h2.ϕa
	# h.ϕb = h2.ϕb
	# h.ϕc = h2.ϕc
	# h.vf = h2.vf
end

print_save("\nβ, RRA, IES: $(round(h.β,2)), $(h.γ), $(h.ψ)")
print_save("\nϵ: $(h.ϵgrid)")
print_save("\nz: $(h.zgrid)")
print_save("\nω: $(h.ωgrid)\n")

# Run
vfi!(h, verbose = true, remote = (location=="remote"))

Void