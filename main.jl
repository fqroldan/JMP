using QuantEcon, BasisMatrices, Interpolations, Optim, NLopt, LaTeXStrings, Distributions, JLD, Sobol, HCubature, Distributed, Dates, ORCA

# Load codes
include("reporting_routines.jl")
include("type_def.jl")
include("interp_atosr.jl")
include("reiter.jl")
include("comp_eqm.jl")
include("gov_pol.jl")
include("simul.jl")
include("handle_guesses.jl")
include("plotting_routines.jl")

#				r_loc,   tax, RRA,     τ,    ρz,    σz,    ρξ,    σξ,  wbar
params_center = [0.09; 0.010;  10; 0.092; 0.970; 0.0035; 0.995; 1e-5; 0.88]

# Set options
nodef     	 = false
noΔ 		 = false
rep_agent 	 = false

# Run
function wrapper_run(params, nodef, noΔ, rep_agent, L, gs; do_all::Bool=true)

	time_init = time()
	
	ρξ, σξ = 0.995, 0.002
	if !do_all
		params = [params[1:6]; ρξ; σξ; params[end]]
	end
	push!(L, length(L)+1)
	run_number = L[end]
	savedir = pwd() * "/../Output/run$(run_number)/"

	s = read("../Output/big_output.txt", String)
	write("../Output/big_output.txt", s * "run number : $(run_number). ")

	# Initialize output file
	write("../Output/output.txt", "\nAggregate Demand and Sovereign Debt Crises\n")

	print_save("\nStarting run number $(run_number) on $(nprocs()) cores and $(Threads.nthreads()) threads at $(Dates.format(now(),"HH:MM")) on $(Dates.monthname(now())) $(Dates.day(now()))")


	r_loc, tax, RRA, τ, ρz, σz, ρξ, σξ, wbar = params
	h = make_guess(nodef, noΔ, rep_agent, r_loc, tax, RRA, τ, ρz, σz, ρξ, σξ, wbar, run_number);

	already_done = false
	try
		# h_done = load(pwd() * "/../Output/run$(run_number)/hank.jld", "h")
		params = load(pwd() * "/../Output/run$(run_number)/params.jld", "params")
		print_save("\nFound params file for run $(run_number).")
		if params == pars(h)
			print_save(" Parameters correct. Looking for g value.")
			try
				g_done = load(pwd() * "/../Output/run$(run_number)/g.jld", "g")
				print_save(" Found g.")
				print_save("\ng = $(g_done)")
				s = read("../Output/big_output.txt", String)
				s *= "g = $g_done"
				push!(gs, g_done)
				if g_done == minimum(gs)
					s *= " ✓"
				end
				s *= "\n"
				write("../Output/big_output.txt", s)
				return g_done
			catch
				print_save(" Couldn't find g.")
			end
		else
			print_save(" Found different parameters, rewriting.")
		end
	catch
		print_save("\nNo previous file found.")
	end

	run(`rm $savedir -rf`)
	run(`mkdir -p $savedir`)

	save(savedir * "params.jld", "params", params)
	params_table = make_params_table(params)
	write(savedir * "params_table.txt", params_table)

	print_save("\nϵ: $(h.ϵgrid)")
	print_save("\nz: $(h.zgrid)")
	print_save("\nξ: $(h.ξgrid)")
	print_save("\nω: $(h.ωgrid)\n")

	mpe_iter!(h; nodef = nodef, noΔ = noΔ, rep_agent = rep_agent, run_number=run_number)
	plot_hh_policies(h, run_number=run_number)
	plot_contour_debtprice(h, savedir)
	plot_contour_unemp(h, savedir)
	
	years = 4000
	g, p_bench, πthres, v_m = make_simulated_path(h, run_number, years)
	run(`cp ../Output/hank.jld ../Output/run$(run_number)/hank.jld`)

	s = read("../Output/big_output.txt", String)
	s *= "g = $(@sprintf("%0.3g",g)) in $(time_print(time()-time_init))"
	push!(gs, g)
	if g == minimum(gs)
		print_save("Minimum g for now. Computing no-def comparison")
		s *= " ✓"

		v_noΔ, v_nodef = make_comparison_simul(h, noΔ, rep_agent, run_number, years, p_bench, "onlyspread", πthres, savedir)

		calib_table_comp = make_calib_table_comp(v_m, v_nodef, v_noΔ)
		write(savedir * "calib_table_comp.txt", calib_table_comp)
	else
		print_save("Suboptimal g. Skipping computation of no-def")
	end

	run(`cp ../Output/run$(run_number)/hank.jld ../Output/hank.jld`)

	s *= "\n"
	write("../Output/big_output.txt", s)

	s = read("../Output/output.txt", String)
	write(savedir * "output.txt", s)

	save(savedir * "g.jld", "g", g)
	params = pars(h)
	save(savedir * "params.jld", "params", params)
	
	return g
end

# wrapper_run(params_center, nodef, noΔ, rep_agent, L)

function SMM(params_center, do_all::Bool=true)
	write("../Output/big_output.txt", "")
	#				 r_loc,   tax,    RRA,     τ,    ρz,    σz,    ρξ,    σξ,    wbar
	# params_center = [0.094; 0.02 ; 12.032; 0.092; 0.970; 0.005; 0.995; 0.002; 0.91]
	if do_all
		mins = 	      [0.05 ; 0.001; 5     ; 0.05 ;  0.85; 0.001;  0.99; 1e-8; 0.82]
		maxs = 		  [0.15 ; 0.05 ; 20    ; 0.35 ;  0.99; 0.012; 0.999; 0.003 ; 1.12]
	else
		mins = 		  [0.05 ; 0.001; 5     ; 0.05 ;  0.85; 0.001; 			    0.82]
		maxs = 		  [0.15 ; 0.05 ; 20    ; 0.35 ;  0.99; 0.012; 			    1.12]
		params_center = [params_center[1:6]; params_center[9]]
	end

	L = Vector{Int64}(undef, 0)
	gs = Vector{Float64}(undef, 0)
	# inner_opt = LBFGS(;linesearch=LineSearches.HagerZhang(linesearchmax=200))
	nlprecon = GradientDescent(alphaguess=Optim.LineSearches.InitialStatic(alpha=1e-4,scaled=true),
                           linesearch=Optim.LineSearches.Static())
	oacc10 = OACCEL(nlprecon=nlprecon, wmax=10)
	res = Optim.optimize(
		params -> wrapper_run(params, false, false, false, L, gs, do_all=do_all)
		# , params_center
		, mins, maxs, params_center, Fminbox(oacc10)
		)

	print("$(res)")
	s = read("../Output/big_output.txt", String)
	write("../Output/big_output.txt", s * "$(res)\n")

	nothing
end

SMM(params_center)
