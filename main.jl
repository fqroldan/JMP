using QuantEcon, BasisMatrices, Interpolations, Optim, NLopt, LaTeXStrings, Distributions, JLD, Sobol, HCubature, Distributed, Dates

# Load codes
include("reporting_routines.jl")
include("type_def.jl")
include("interp_atosr.jl")
include("reiter.jl")
include("comp_eqm.jl")
include("gov_pol.jl")
include("simul.jl")
include("handle_guesses.jl")
# include("plotting_routines.jl")

#				 r_loc,   tax,    RRA,     τ,    ρz,    σz,    ρξ,    σξ,    wbar
params_center = [0.094; 0.02 ; 12.032; 0.092; 0.875; 0.007; 0.995; 0.002; 1.10825]

# Set options
nodef     	 = false
noΔ 		 = false
rep_agent 	 = false

# Run
function wrapper_run(params, nodef, noΔ, rep_agent, L)
	push!(L, length(L)+1)
	run_number = L[end]
	saving_folder = pwd() * "/../Output/run$(run_number)/"
	run(`mkdir -p $saving_folder`)

	# Initialize output file
	write("../Output/output.txt", "")

	print_save("\nAggregate Demand and Sovereign Debt Crises\n")
	print_save("\nStarting run number $(run_number) on $(nprocs()) cores and $(Threads.nthreads()) threads at $(Dates.format(now(),"HH:MM")) on $(Dates.monthname(now())) $(Dates.day(now()))")

	save(saving_folder * "params.jld", "params", params)

	r_loc, tax, RRA, τ, ρz, σz, ρξ, σξ, wbar = params
	h = make_guess(nodef, noΔ, rep_agent, r_loc, tax, RRA, τ, ρz, σz, ρξ, σξ, wbar, run_number);

	print_save("\nϵ: $(h.ϵgrid)")
	print_save("\nz: $(h.zgrid)")
	print_save("\nξ: $(h.ξgrid)")
	print_save("\nω: $(h.ωgrid)\n")

	mpe_iter!(h; nodef = nodef, noΔ = noΔ, rep_agent = rep_agent, run_number=run_number)
	g = make_simulated_path(h, run_number)

	s = read("../Output/output.txt", String)
	write(saving_folder * "output.txt", s)

	return g
end

L = Vector{Int64}(undef, 0)
wrapper_run(params_center, nodef, noΔ, rep_agent, L)

nothing
