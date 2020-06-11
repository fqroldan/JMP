using QuantEcon, BasisMatrices, Interpolations, Optim, LaTeXStrings, Distributions, JLD, HCubature, Distributed, Dates, ForwardDiff, Printf, Random, LinearAlgebra, DataFrames, GLM, PlotlyJS, ColorSchemes
using ORCA, JSON

include("type_def.jl")
include("handle_guesses.jl")
include("handle_itps.jl")
include("fiscal.jl")
include("hh_pb.jl")
include("comp_eqm.jl")
include("gov_pol.jl")
include("reporting_routines.jl")
include("simul.jl")
include("plotting_routines.jl")

# print("mpe_iter!(sd)")
params_center = Dict{Symbol, Float64}(
	:β		=> 0.9770096825552423,
	:γ		=> 13,
	:τ		=> 0.24396769268442076,
	:wbar	=> 0.897,
	:ρz		=> 0.97,
	:σz		=> 0.0022,
	:meanξ	=> 0.0018,
	:ρξ		=> 0.95,
	:σξ		=> 0.001,
)

function wrapper_run(par_vec, nodef, noΔ, rep_agent, L, gs; do_all::Bool=true)
	params = Dict{Symbol, Float64}(
		:β		=> par_vec[1],
		:γ		=> par_vec[2],
		:τ		=> par_vec[3],
		:wbar	=> par_vec[4],
		:ρz		=> par_vec[5],
		:σz		=> par_vec[6],
		:meanξ	=> par_vec[7],
		:ρξ		=> par_vec[8],
		:σξ		=> par_vec[9],
	)

	time_init = time()
	
	if !do_all
		params[:ρξ] = 0.995
		params[:σξ] = 0.002
		params[:τ]  = 0.092
		params[:ρz] = 0.970
	end
	push!(L, length(L)+1)
	run_number = L[end]
	savedir = pwd() * "/../Output/run$(run_number)/"

	s = read("../Output/big_output.txt", String)
	write("../Output/big_output.txt", s * "run number : $(run_number). ")

	# Initialize output file
	write("../Output/output.txt", "\nAggregate Demand and Sovereign Debt Crises\n")

	print_save("\nStarting run number $(run_number) on $(nprocs()) cores and $(Threads.nthreads()) threads at $(Dates.format(now(),"HH:MM")) on $(Dates.monthname(now())) $(Dates.day(now()))")

	already_done, g = determine_done(run_number, params)
	if already_done
		s = read("../Output/big_output.txt", String)
		s *= "g = $(@sprintf("%0.3g",g)) in $(time_print(time()-time_init))"
		push!(gs, g)
		if g == minimum(gs)
			s *= " ✓"
		end

		s *= "\n"
		write("../Output/big_output.txt", s)
		return g
	end

	sd = make_guess(nodef, noΔ, rep_agent, params, run_number);

	if !already_done
		run(`rm $savedir -rf`)
		run(`mkdir -p $savedir`)

		save(savedir * "params.jld", "params", pars(sd))
		params_table = make_params_table(sd)
		write(savedir * "params_table.txt", params_table)

		mpe_iter!(sd; nodef = nodef, noΔ = noΔ, rep_agent = rep_agent, run_number=run_number, save_copies=true)
		# plot_hh_policies(h, run_number=run_number)
		# plot_contour_debtprice(h, savedir)
		# plot_contour_unemp(h, savedir)
		
		years = 10000
		g, p_bench, πthres, v_m, def_freq = make_simulated_path(sd, savedir, years)
		save(pwd() * "/../Output/g.jld", "g", g)
		run(`cp ../Output/SOEdef.jld ../Output/run$(run_number)/SOEdef.jld`)
	else
		p_bench = load("../Output/run$(run_number)/p_bench.jld", "pp")
		v_m = simul_stats(p_bench)
		ζ_vec = series(p_bench,:ζ)
		Ndefs = sum([ (ζ_vec[jj] == 0) & (ζ_vec[jj-1] == 1) for jj in 2:length(ζ_vec)])
		years = floor(Int64,periods(p_bench)*0.25)
		def_freq = Ndefs / years
		# _, πthres = plot_episodes(p_bench; episode_type="onlyspread", slides=true, πthres=0.95)
	end

	πthres = 1.0

	s = read("../Output/big_output.txt", String)
	s *= "g = $(@sprintf("%0.3g",g)) in $(time_print(time()-time_init))"
	push!(gs, g)
	if g == minimum(gs)
		print_save("Minimum g for now. Computing no-def comparison")
		s *= " ✓"
		write("../Output/big_output.txt", s)
		
		if length(gs) > 1
			current_best = findmin(gs[1:end-1])[2]
		else
			current_best = 1
		end

		v_noΔ, v_nodef, v_nob, freq_noΔ, freq_nodef, freq_nob = make_comparison_simul(sd, noΔ, rep_agent, run_number, current_best, years, p_bench, "onlyspread", πthres, savedir, already_done)

		calib_table_comp = make_calib_table_comp([v_m; 100*def_freq], [v_nodef; 100*freq_nodef], [v_noΔ; 100*freq_noΔ], [v_nob; 100*freq_nob])
		write(savedir * "calib_table_comp.txt", calib_table_comp)
	else
		print_save("Suboptimal g. Skipping computation of no-def")
	end

	if !already_done
		run(`cp ../Output/run$(run_number)/SOEdef.jld ../Output/SOEdef.jld`)
		run(`cp ../Output/output.txt ../Output/run$(run_number)/output.txt`)
	end

	s *= "\n"
	write("../Output/big_output.txt", s)


	# s = read("../Output/output.txt", String)
	# write(savedir * "output.txt", s)

	params = pars(sd)
	save(savedir * "params.jld", "params", params)
	run(`cp ../Output/g.jld ../Output/run$(run_number)/g.jld`)

	return g
end

function SMM(p_dict; do_all::Bool=true)
	write("../Output/big_output.txt", "")
	params_center = [p_dict[:β], p_dict[:γ], p_dict[:τ], p_dict[:wbar], p_dict[:ρz], p_dict[:σz], p_dict[:meanξ], p_dict[:ρξ], p_dict[:σξ]]
	#				 
	# params_center = [0.094; 0.02 ; 12.032; 0.092; 0.970; 0.005; 0.995; 0.002; 0.91]
	if do_all
		mins = 	      [1.15^(-0.25) ; 5   ; 0.05 ; 0.82 ; 0.85 ; 0.0001; 0.0001; 0.92 ; 1e-8 ]
		maxs = 		  [1.05^(-0.25) ; 20  ; 0.35 ; 1.00 ; 0.99 ; 0.008 ; 0.02  ; 0.999; 0.003]
	else
		mins = 	      [1.15^(-0.25) ; 5   ;      ; 0.82 ;      ; 0.0001; 0.0001;      ;      ]
		maxs = 		  [1.05^(-0.25) ; 20  ;      ; 1.00 ;      ; 0.012 ; 0.05  ;      ;      ]
		params_center = [p_dict[:β], p_dict[:γ], p_dict[:wbar], p_dict[:σz], p_dict[:meanξ]]
	end

	L = Vector{Int64}(undef, 0)
	gs = Vector{Float64}(undef, 0)
	# inner_opt = LBFGS(;linesearch=LineSearches.HagerZhang(linesearchmax=200))
	nlprecon = GradientDescent(alphaguess=Optim.LineSearches.InitialStatic(alpha=1e-4,scaled=true),
                           linesearch=Optim.LineSearches.Static())
	inner_opt = OACCEL(nlprecon=nlprecon, wmax=10)

	inner_opt = NelderMead()
	res = Optim.optimize(
		params -> wrapper_run(params, false, false, false, L, gs, do_all=do_all)
		# , params_center
		, mins, maxs, params_center, Fminbox(inner_opt)
		)

	print("$(res)")
	s = read("../Output/big_output.txt", String)
	write("../Output/big_output.txt", s * "$(res)\n")

	nothing
end
