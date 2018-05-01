using QuantEcon, BasisMatrices, Interpolations, Optim, NLopt, MINPACK, LaTeXStrings, Distributions, JLD

location = "remote"
if pwd() == "/home/q/Dropbox/NYU/AToSR/Codes"
	location = "local"
	using Rsvg
end

remote = (location=="remote")

if remote
	dir_path = "/../../"
else
	dir_path = "/../"
end

# Initialize output file
write(pwd()*dir_path*"output.txt", "")

# Load codes
@everywhere include("reporting_routines.jl")
@everywhere include("type_def.jl")
@everywhere include("interp_atosr.jl")
@everywhere include("reiter.jl")
@everywhere include("comp_eqm.jl")
include("simul.jl")
include("plotting_routines.jl")

print_save("\nAggregate Demand around Debt Crises\n", remote=remote)

print_save("\nStarting $(location) run on $(nprocs()) cores at "*Dates.format(now(),"HH:MM"), remote=remote)

# Set options
local_run = true

# Initialize type
if remote || local_run
	h = Hank();
	try
		h2 = load("hank.jld", "h")
		if h.ψ == h2.ψ && h.γ == h2.γ && h.Ns == h2.Ns
			print_save("Starting from loaded guess", remote=remote)
			h.ϕa = h2.ϕa
			h.ϕb = h2.ϕb
			h.ϕc = h2.ϕc
			h.vf = h2.vf
			h.pN = h2.pN
			h.μ′ = h2.μ′
			h.σ′ = h2.σ′
			h.wage = h2.wage
			h.Ld = h2.Ld
		end
	end
else
	print_save("\nLoading solved model file\n", remote=remote)
	h = load("../HPC_Output/hank.jld", "h")
end


print_save("\nβ, RRA, IES: $(round(h.β,2)), $(h.γ), $(h.ψ)", remote=remote)
print_save("\nϵ: $(h.ϵgrid)", remote=remote)
print_save("\nz: $(h.zgrid)", remote=remote)
print_save("\nω: $(h.ωgrid)\n", remote=remote)

# Run
if remote || local_run 
	vfi!(h, verbose = true, remote = remote)
	save(pwd() * "/hank.jld", "h", h)
end

p, p_full, jz_series, ols = simul(h; simul_length=4*25, burn_in=4*250, only_def_end=true)
save(pwd()*"/ols.jld", "ols", ols)

plot_simul(p_full; remote=remote, name="_full")
plot_simul(p; remote=remote)

Void