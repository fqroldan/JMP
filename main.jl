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
include("gov_pol.jl")
include("simul.jl")
# include("plotting_routines.jl")

print_save("\nAggregate Demand around Debt Crises\n")

print_save("\nStarting $(location) run on $(nprocs()) cores at "*Dates.format(now(),"HH:MM"))

# Set options
local_run = true
# Initialize type
if remote || local_run
	h = Hank();
	# h = load(pwd() * "/../../hank.jld", "h")
	try
		h2 = load(pwd() * "/../../hank.jld", "h")
		remote? h2 = load(pwd() * "/../../hank.jld", "h"): h2 = load("hank.jld", "h")
		print_save("\nFound JLD file")
		if h.Ns == h2.Ns && h.Nω == h2.Nω && h.Nϵ == h2.Nϵ
			print_save(": loading previous results")
			h.ϕa = h2.ϕa
			h.ϕb = h2.ϕb
			h.ϕc = h2.ϕc
			h.vf = h2.vf
			h.pngrid = h2.pngrid
			h.wgrid = h2.wgrid
			h.pN = h2.pN
			h.μ′ = h2.μ′
			h.σ′ = h2.σ′
			h.w′ = h2.w′
			h.output = h2.output
			h.wage = h2.wage
			h.Ld = h2.Ld
			h.repay = h2.repay
			h.welfare = h2.welfare
			print_save(" ✓")
		end
	catch
		print_save("JLD file incompatible")
	end
else
	print_save("\nLoading solved model file\n")
	h = load("../HPC_Output/hank.jld", "h")
end

print_save("\nβ, RRA, IES: $(round(h.β,2)), $(h.γ), $(h.ψ)")
print_save("\nϵ: $(h.ϵgrid)")
print_save("\nz: $(h.zgrid)")
print_save("\nω: $(h.ωgrid)\n")

# Run
if remote || local_run
	# vfi!(h, verbose = true, remote = remote)
	mpe_iter!(h; remote = remote)
end

p, jz_series = simul(h; simul_length=4*(250+25), only_def_end=true)
# plot_simul(p; trim = 0, remote=remote)
# plot_simul(p; trim = 4*250, remote=remote)

ols = simul_regs(p)
save(pwd()*"/ols.jld", "ols", ols)

Void
