make_guess(nodef, noΔ, rep_agent, p_dict, run_number) = make_guess(nodef, noΔ, rep_agent, p_dict[:β], p_dict[:meanξ], p_dict[:γ], p_dict[:τ], p_dict[:ρz], p_dict[:σz], p_dict[:ρξ], p_dict[:σξ], p_dict[:wbar], run_number)

function reinterp_b(sd::SOEdef, y, new_bgridzz; agg::Bool=false)
	gr = sd.gr
	knots = (gr[:ω], gr[:ϵ], gr[:b], gr[:μ], gr[:σ], gr[:ξ], gr[:ζ], gr[:z])
	if agg
		knots = (gr[:b], gr[:μ], gr[:σ], gr[:ξ], gr[:ζ], gr[:z])
		y = reshape(y, N(sd,:b), N(sd,:μ), N(sd,:σ), N(sd,:ξ), N(sd,:ζ), N(sd,:z))
	end

	itp_obj_y = interpolate(knots, y, Gridded(Linear()))
	itp_y = extrapolate(itp_obj_y, Line())

	if agg
		y_new = itp_y(new_bgrid, gr[:μ], gr[:σ], gr[:ξ], gr[:ζ], gr[:z])
		return reshape(y_new, length(y_new))
	else
		y_new = itp_y(gr[:ω], gr[:ϵ], new_bgrid, gr[:μ], gr[:σ], gr[:ξ], gr[:ζ], gr[:z])
		return y_new
	end
end

function reinterp_b!(sd::SOEdef, sd_old::SOEdef, nodef::Bool)
	for (key, val) in sd_old.ϕ
		sd.ϕ[key] = reinterp_b(sd_old, val, sd.gr[:b])
	end
	for (key, val) in sd_old.v
		sd.v[key] = reinterp_b(sd_old, val, sd.gr[:b])
	end
	for (key, val) in sd_old.eq
		sd.eq[key] = reinterp_b(sd_old, val, sd.gr[:b], agg=true)
	end

	μ′_new = similar(sd.LoM[:μ])
	σ′_new = similar(sd.LoM[:σ])

	for (jξp, ξpv) in enumerate(sd.gr[:ξ]), (jzp, zpv) in enumerate(sd.gr[:z])
		mvec = [reinterp_b(sd_old, [sd_old.LoM[:μ][js,jξp,jzp][jζp] for js in 1:size(sd_old.LoM[:μ],1)], sd.gr[:b], agg=true) for jζp in 1:2]
		svec = [reinterp_b(sd_old, [sd_old.LoM[:σ][js,jξp,jzp][jζp] for js in 1:size(sd_old.LoM[:σ],1)], sd.gr[:b], agg=true) for jζp in 1:2]
		for js in 1:size(μ′_new,1)
			μ′_new[js,jξp,jzp] = [max(min(mvec[jζp][js], maximum(sd.gr[:μ])), minimum(sd.gr[:μ])) for jζp in 1:2]
			σ′_new[js,jξp,jzp] = [max(min(svec[jζp][js], maximum(sd.gr[:σ])), minimum(sd.gr[:σ])) for jζp in 1:2]
		end
	end
	sd.LoM[:μ] = μ′_new
	sd.LoM[:σ] = σ′_new

	if nodef
		print_save(" Not loading default policy.")
	else
		knots 		= (sd_old.gr[:b], sd_old.gr[:μ], sd_old.gr[:σ], sd_old.gr[:ξ], sd_old.gr[:ζ], sd_old.gr[:z], sd_old.gr[:ξ], sd_old.gr[:z])
		repay_mat 	= reshape(sd_old.gov[:repay], N(sd_old,:b), N(sd_old,:μ), N(sd_old,:σ), N(sd_old,:ξ), N(sd_old,:ζ), N(sd_old,:z), N(sd_old,:ξ), N(sd_old,:z))
		itp_repay 	= extrapolate(interpolate(knots, repay_mat, Gridded(Linear())), Line())
		rep_new 	= itp_repay(sd.gr[:b], sd.gr[:μ], sd.gr[:σ], sd.gr[:ξ], sd.gr[:ζ], sd.gr[:z], sd.gr[:ξ], sd.gr[:z])
		rep_new 	= max.(0, min.(1, rep_new))
		sd.gov[:repay] 	= reshape(rep_new, length(rep_new))
	end
end

function make_guess(nodef, noΔ, rep_agent, β, tax, RRA, τ, ρz, σz, ρξ, σξ, wbar, run_number)

	print_save("\nRun with β, RRA, τ, wbar, ρz, σz, tax, ρξ, σξ = $(round(100*(β^-4-1),digits=4))%, $(round(RRA,digits=4)), $(round(τ,digits=4)), $(round(wbar,digits=4)), $(round(ρz,digits=4)), $(round(σz,digits=4)), $(round(tax,digits=4)), $(round(ρξ,digits=4)), $(round(σξ,digits=4))")
	sd = SOEdef(; β=β, tax = tax, RRA=RRA, τ=τ, nodef = nodef, noΔ = noΔ, rep_agent = rep_agent, ρz=ρz, σz=σz, ρξ=ρξ, σξ=σξ, wbar=wbar
		# , Nω=2,Nϵ=3,Nb=2,Nμ=2,Nσ=2,Nξ=2,Nz=3
		);
	try
		# h2 = load(pwd() * "/../Output/hank_backup.jld", "h")
		sd_old = load(pwd() * "/../Output/SOEdef.jld", "sd")
		print_save("\nFound generic JLD file")
		try
			sd_old = load(pwd() * "/../Output/run$(max(1,run_number-1))/SOEdef.jld", "sd") # set max(1, ...) to use run1 in first go
			print_save("\nFound JLD file from last run")
		catch
		end
		if N(sd,:ω) == N(sd_old,:ω) && N(sd,:ϵ) == N(sd_old,:ϵ)
			print_save(": loading previous results")
			if N(sd,:μ) == N(sd_old,:μ) && N(sd,:σ)==N(sd,:σ)
				sd.gr[:μ] = sd_old.gr[:μ]
				sd.gr[:σ] = sd_old.gr[:σ]
				sd.gr[:pN] = sd_old.gr[:pN]
				print_save("...")
			end


			if N(sd_old, :b) == N(sd, :b)
				for (key, val) in sd_old.ϕ
					sd.ϕ[key] = val
				end
				for (key, val) in sd_old.v
					sd.v[key] = val
				end

				for (key, val) in sd_old.eq
					sd.eq[key] = val
				end
				for (key, val) in sd_old.LoM
					sd.LoM[key] = val
				end
				
				if nodef
					print_save(" Not loading default policy.")
				else
					for (key, val) in sd_old.gov
						sd.gov[key] = val
					end
				end
			else
				reinterp_b!(sd, sd_old, nodef)
			end
			print_save(" ✓")
		end
	catch
		print_save(": JLD file incompatible")
	end
	return sd
end

function eval_GMM(v_o, target_o = vec([0.96580506  0.01294576  0.96172496  0.01663608  0.96656486  0.10252351 64.57638889 23.48323041 15.94722222  6.08732167  56.4851069778397  94.479167]); show_res::Bool=true)
	v_m = copy(v_o)
	targets = copy(target_o)

	targets[4] = targets[4] / targets[2]
	v_m[4] = v_m[4] / v_m[2]
	print("\nTargets ✓")

	W = diagm([1.0/targets[jj] for jj in 1:length(targets)])
	W[2,2] *= 100
	W[4,4] *= 50

	g = (v_m - targets)'*W*(v_m-targets)

	if show_res
		targetnames = ["AR(1) Output"; "σ(Output)"; "AR(1) Cons"; "σ(Cons) / σ(Output)"; "AR(1) Spreads"; "σ(spreads)"; "mean B/Y"; "std B/Y"; "mean unemp"; "std unemp"; "median Dom Holdings"; "mean wealth/Y" ]
		res = [targetnames v_m targets (targets-v_m)./targets]
		for jj in 1:size(res,1)
		    print_save("\n")
		    print_save(rpad(res[jj,1],20," "))
		    for ii in 2:size(res,2)
		    	print_save(rpad("$(@sprintf("%0.3g",res[jj,ii]))",9," "))
		    end
		end
		print_save("\ng = $(@sprintf("%0.3g", g))")
	end
	print_save("\n")
	return g
end

function make_simulated_path(sd::SOEdef, savedir, years=100)
	pp, Ndefs = parsimul(sd; simul_length=4*years, burn_in=1+4*100)
	Tyears = floor(Int64,periods(pp)*0.25)
	def_freq = Ndefs/Tyears
	print_save("\n$Ndefs defaults in $Tyears years: default freq = $(round(1000*def_freq)/10)%")
	print_save("\nAverage Gini coefficient: $(@sprintf("%0.3g",100*mean([mean(series(path,:Gini)) for path in pp])))")
	save(savedir*"p_bench.jld", "pp", pp, "Ndefs", Ndefs)
	
	# pl = plot_simul(path)

	# pl, πthres = plot_episodes(path; episode_type="onlyspread", slides=true, πthres=0.95) # Also here need to make it so there are at least 10 episodes
	πthres = 1.0
	# make_IRF_plots(path; slides = true, create_plots = true, response = resp, savedir=savedir) # for resp = Y, C

	v_m = simul_stats(pp)
	targets = vec([0.96580506  0.01294576  0.96172496  0.01663608  0.96656486  0.10252351 64.57638889 23.48323041 15.94722222  6.08732167  56.4851069778397  94.479167])
	
	g = eval_GMM(v_m, targets)
	calib_table = make_calib_table(v_m)
	write(savedir * "calib_table.txt", calib_table)

	return g, pp, πthres, v_m, def_freq
end
function pass_CE!(sd::SOEdef, sdg::SOEdef)
	sd.gr[:μ] = sdg.gr[:μ]
	sd.gr[:σ] = sdg.gr[:σ]
	sd.gr[:pN] = sdg.gr[:pN]

	for (key, val) in sdg.ϕ
		sd.ϕ[key] = val
	end
	for (key, val) in sdg.v
		sd.v[key] = val
	end

	for (key, val) in sdg.eq
		sd.eq[key] = val
	end
	for (key, val) in sdg.LoM
		sd.LoM[key] = val
	end
end

function try_load(run_number, current_best, sim_name)
	sd = load("../Output/run$(run_number)/SOEdef.jld", "sd")
	try
		print_save("\nLoading results from run $(current_best): ")
		sdg = load("../Output/run$(current_best)/SOEdef_$(sim_name).jld", "sd")
		print_save("found old $(sim_name) results")
		pass_CE!(sd, sdg)
		print_save(" ✓")
	catch
		try
			print_save("\nLoading generic results: ")
			sdg = load("../Output/SOEdef_$(sim_name).jld", "sd")
			print_save("found old $(sim_name) results")
			pass_CE!(sd, sdg)
			print_save(" ✓")
		catch
		end
	end
	return sd
end

function simul_done(run_number, sim_name)
	dir = "Output/run$(run_number)"
	pp, Ndefs = load("../$(dir)/p_$(sim_name).jld", "pp", "Ndefs")
	print_save("\nFound $(sim_name) simul")
	return pp, Ndefs
end

function try_simul(run_number, current_best, sim_name, nodef, nodelta, nob, rep_agent, years, already_done, pars, gov)

	if already_done
		try
			pp, Ndefs = simul_done(run_number, sim_name)
			return pp, Ndefs
		catch
		end
	end

	sd = try_load(run_number, current_best, sim_name)
	for (key, val) in pars
		sd.pars[key] = val
	end
	for (key, val) in gov
		sd.gov[key] = val
	end
	update_probs!(sd)
	if sim_name == "nodelta"
		sd.pars[:Δ] = 0
	end
	print_save("\nSolving $(sim_name) version")
	mpe_iter!(sd; nodef = nodef, noΔ = nodelta, nob = nob, rep_agent = rep_agent, run_number=run_number, save_copies=false)
	pp, Ndefs = parsimul(sd; simul_length=4*years, burn_in=1+4*100)
	save("../Output/run$(run_number)/SOEdef_$(sim_name).jld", "sd", sd)
	save("../Output/SOEdef_$(sim_name).jld", "sd", sd)
	save("../Output/run$(run_number)/p_$(sim_name).jld", "pp", pp, "Ndefs", Ndefs)
	print_save(" ✓")
	return pp, Ndefs
end

make_comparison_simul(sd::SOEdef, noΔ, rep_agent, run_number, current_best, years, p_bench::Path, πthres, savedir, already_done) = make_comparison_simul(sd, noΔ, rep_agent, run_number, current_best, years, [p_bench], πthres, savedir, already_done)
function make_comparison_simul(sd::SOEdef, noΔ, rep_agent, run_number, current_best, years, p_bench::Vector{T}, πthres, savedir, already_done) where T <: AbstractPath

	old_Δ = copy(sd.pars[:Δ])
	sim_mat = [ji==jj for ji in 1:3, jj in 1:3]
	sim_names = ["nodelta", "nodef", "nob"]

	freq = [0.0 for jj in 1:3]
	Wr = [0.0 for jj in 1:3]
	v = [zeros(12) for jj in 1:3]
	for (js, sim_name) in enumerate(sim_names)
		nodelta, nodef, nob = sim_mat[js, :]
		pars, gov = sd.pars, sd.gov
		pp, Ndefs = try_simul(run_number, current_best, sim_name, nodef, nodelta, nob, rep_agent, years, already_done, pars, gov)
		Tyears = floor(Int64,periods(pp)*0.25)
		freq[js] = Ndefs/Tyears
		vsim = simul_stats(pp)
		v[js] = vsim

		Wr[js] = mean([mean(series(p, :Wr)) for p in pp])

		eval_GMM(vsim)
		sd.pars[:Δ] = old_Δ
	end

	freq_nodelta, freq_nodef, freq_nob = freq
	v_nodelta, v_nodef, v_nob = v
	W_nodelta, W_nodef, W_nob = Wr

	# for (jj, slides) in enumerate([true; false])
	# 	pcomp = plot_comparison_episodes(p_bench, p_nodef; episode_type = episode_type, slides = slides, πthres=πthres)
	# 	savejson(pcomp, savedir * "comparison_crisis_nodef$(jj).json")
	# end

	return v_nodelta, v_nodef, v_nob, freq_nodelta, freq_nodef, freq_nob, W_nodelta, W_nodef, W_nob
end

function determine_done(run_number, params)

	already_done = false
	g = zeros(12)
	try
		pars_new = load(pwd() * "/../Output/run$(run_number)/params.jld", "params")
		print_save("\nFound params file for run $(run_number).")
		if pars_new == params
			print_save(" Parameters correct. Looking for g value.")
			try
				g = load(pwd() * "/../Output/run$(run_number)/g.jld", "g")
				print_save(" Found g.")
				print_save("\ng = $(g)")
				print_save("\nLooking for path")
				path = load("../Output/run$(run_number)/p_bench.jld", "pp")
				print_save(": ✓")
				already_done = true
			catch
				print_save(" Couldn't find g.")
			end
		else
			print_save(" Found different parameters, rewriting.")
		end
	catch
		print_save("\nNo previous file found.")
	end

	return already_done, g
end

# 	pp, Ndefs = load("../Output/run$(run_number)/p_nodelta.jld", "pp", "Ndefs")

# 	h.Δ = old_Δ
# 	try 
# 		pp, Ndefs = load("../Output/run$(run_number)/p_nodef.jld", "pp", "Ndefs")
# 		print_save("\nFound nodef simul")
# 	catch
# 		print_save("\nSolving nodef version")
# 		h = load(savedir * "hank.jld", "h")
# 		mpe_iter!(h; nodef = true, noΔ = false, rep_agent = rep_agent, run_number=run_number, maxiter = 21, save_copies=false)
# 		pp, _, Ndefs = simul(h; simul_length=4*(years+25), only_def_end=false)
# 		save("../Output/run$(run_number)/p_nodef.jld", "pp", pp, "Ndefs", Ndefs)
# 		print_save(" ✓")
# 	end
# 	pp, Ndefs = load("../Output/run$(run_number)/p_nodef.jld", "pp", "Ndefs")
# 	freq_nodef = Ndefs/Tyears

# 	v_nodef = simul_stats(p_nodef)

# 	h.Δ = old_Δ
# 	try
# 		p_nob, Ndefs = load("../Output/run$(run_number)/p_nob.jld", "p_nob", "Ndefs")
# 		print_save("\nFound nob simul")
# 	catch
# 		print_save("\nSolving nob version")
# 		h = load(savedir * "hank.jld", "h")
# 		mpe_iter!(h; nodef = false, noΔ = false, rep_agent = rep_agent, run_number=run_number, maxiter = 21, save_copies=false, nob=true)
# 		p_nob, _, Ndefs = simul(h; simul_length=4*(years+25), only_def_end=false)
# 		save("../Output/run$(run_number)/p_nob.jld", "p_nob", p_nob, "Ndefs", Ndefs)
# 		print_save(" ✓")
# 	end
# 	p_nob, Ndefs = load("../Output/run$(run_number)/p_nob.jld", "p_nob", "Ndefs")
# 	freq_nob = Ndefs/Tyears

# 	v_nob = simul_stats(p_nob)

# end


# function make_something_else()
# 	savedir = pwd() * "/../Output/run$(run_number)/"
# 	save(savedir * "path.jld", "pp", pp)
	
# 	pl = plot_simul(path)
# 	savejson(pl, savedir * "path.json")	
# 	try
# 		savefig(savedir * "path.pdf")
# 	catch
# 	end

# 	pl, πthres = plot_episodes(path; episode_type="onlyspread", slides=true, πthres=0.95)
# 	savejson(pl, savedir * "onlyspread_slides.json")
# 	try
# 		savefig(savedir * "onlyspread_slides.pdf")
# 	catch
# 	end

# 	# pl = plot_episodes(path; episode_type="onlyspread", slides=false, πthres=0.975)
# 	# savejson(pl, savedir * "onlyspread_paper.json")

# 	for resp in ["C"; "Y"]
# 		try
# 			make_IRF_plots(path; slides = true, create_plots = true, response = resp, savedir=savedir)
# 		catch
# 		end
# 	end

# 	params = pars(h)

# 	try
# 		v_m = simul_stats(path)
# 		targetnames = ["AR(1) Output"; "σ(Output)"; "AR(1) Cons"; "σ(Cons) / σ(Output)"; "AR(1) Spreads"; "σ(spreads)"; "mean B/Y"; "std B/Y"; "mean unemp"; "std unemp"; "median Dom Holdings"; "mean wealth/Y" ]
# 		targets[4] = targets[4] / targets[2]
# 		v_m[4] = v_m[4] / v_m[2]
# 		print("\nTargets ✓")

# 		W = zeros(length(v_m),length(v_m))
# 		[W[jj,jj] = 1.0/targets[jj] for jj in 1:length(targets)]
# 		W[2,2] *= 100
# 		W[4,4] *= 50

# 		print("\nW ✓")
# 		g = (v_m - targets)'*W*(v_m-targets)
# 		print_save("\nObjective function = $(@sprintf("%0.3g",g))")
# 		print_save("\n")
# 		for jj in 1:length(targets)
# 			print_save("$(@sprintf("%0.3g",v_m[jj]))  ")
# 		end
# 		print_save("\nParams: ")
# 		for jj in 1:length(params)
# 			print_save("$(@sprintf("%0.3g",params[jj]))  ")
# 		end
# 		res = [targetnames v_m targets (targets-v_m)./targets]
# 		for jj in 1:size(res,1)
# 		    print_save("\n")
# 		    print_save(res[jj,1])
# 		    for ii in 2:size(res,2)
# 		    	print_save("     ")
# 		    	print_save("$(@sprintf("%0.3g",res[jj,ii]))")
# 		    end
# 		end
# 	catch
# 		print_save("\nWARNING: Found problems computing simulation statistics")
# 	end

# 	v_m[4] = v_m[4] * v_m[2]
# 	# save(savedir * "v_m.jld", "v_m", v_m)
# 	# write(savedir * "stats.txt", "$(v_m)")
	
# 	calib_table = make_calib_table(v_m)
# 	write(savedir * "calib_table.txt", calib_table)

# 	return g, path, πthres, v_m, def_freq
# end

# pars(sd::SOEdef) = [(1/sd.pars[:β])^4-1; sd.pars[:γ]; sd.pars[:τ]; sd.pars[:wbar]; sd.pars[:ρz]; sd.pars[:σz]; sd.pars[:meanξ]; sd.pars[:ρξ]; sd.pars[:σξ]]

pars(sd::SOEdef) = Dict(
	:β 		=> sd.pars[:β],
	:γ		=> sd.pars[:γ],
	:τ		=> sd.pars[:τ],
	:wbar	=> sd.pars[:wbar],
	:ρz		=> sd.pars[:ρz],
	:σz		=> sd.pars[:σz],
	:meanξ	=> sd.pars[:meanξ],
	:ρξ		=> sd.pars[:ρξ],
	:σξ		=> sd.pars[:σξ],
	)