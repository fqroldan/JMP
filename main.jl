using QuantEcon, BasisMatrices, Interpolations, Optim, NLopt, MINPACK, LaTeXStrings, Distributions, JLD

location = "remote"
if pwd() == "/home/q/Dropbox/NYU/AToSR/Codes"
	location = "local"
	using Rsvg
end

remote = (location=="remote")

# Initialize output file
write(pwd()*"/../../output.txt", "")

# Load codes
@everywhere include("reporting_routines.jl")
@everywhere include("type_def.jl")
@everywhere include("interp_atosr.jl")
@everywhere include("reiter.jl")
@everywhere include("comp_eqm.jl")
include("simul.jl")
include("plotting_routines.jl")

print_save("\nAggregate Demand around Debt Crises\n")

print_save("\nStarting $(location) run on $(nprocs()) cores at "*Dates.format(now(), "HH:MM"))

# Initialize type
h = Hank();
# h = load("hank.jld", "h")
try
	h2 = load("hank.jld", "h")
	if h.ψ == h2.ψ && h.γ == h2.γ && h.Ns == h2.Ns
		h.ϕa = h2.ϕa
		h.ϕb = h2.ϕb
		h.ϕc = h2.ϕc
		h.vf = h2.vf
		h.pN = h2.pN
		h.wage = h2.wage
		h.Ld = h2.Ld
	end
end

# h = load("../HPC_Output/hank.jld", "h")

print_save("\nβ, RRA, IES: $(round(h.β,2)), $(h.γ), $(h.ψ)")
print_save("\nϵ: $(h.ϵgrid)")
print_save("\nz: $(h.zgrid)")
print_save("\nω: $(h.ωgrid)\n")

# Run
vfi!(h, verbose = true, remote = remote)
save(pwd() * "/../../hank.jld", "h", h)

p, jz_series = simul(h; simul_length=200, burn_in=100, only_def_end=true)
plot_simul(p; remote=remote)

Void