function make_guess(nodef, noΔ, rep_agent, r_loc, tax, RRA, τ, ρz, σz, ρξ, σξ, wbar, run_number)

	print_save("\nRun with r_loc, RRA, τ, wbar, ρz, σz, tax, ρξ, σξ = $(round(r_loc,digits=3)), $(round(RRA,digits=3)), $(round(τ,digits=3)), $(round(wbar,digits=3)), $(round(ρz,digits=3)), $(round(σz,digits=3)), $(round(tax,digits=3)), $(round(ρξ,digits=3)), $(round(σξ,digits=3))")
	h = Hank(; β=(1.0/(1.0+r_loc))^0.25, tax = tax, RRA=RRA, τ=τ, nodef = nodef, noΔ = noΔ, rep_agent = rep_agent, ρz=ρz, σz=σz, ρξ=ρξ, σξ=σξ, wbar=wbar
		# , Nω=2,Nϵ=3,Nb=2,Nμ=2,Nσ=2,Nξ=2,Nz=3
		);
	try
		h2 = load(pwd() * "/../Output/hank.jld", "h")
		try
			h2 = load(pwd() * "/../Output/run$(run_number-1)/hank.jld", "h")
			print_save("\nFound JLD file from last run")
		catch
			print_save("\nFound generic JLD file")
		end
		if h.Nω == h2.Nω && h.Nϵ == h2.Nϵ
			print_save(": loading previous results")
			h.μgrid = h2.μgrid
			h.σgrid = h2.σgrid
			update_grids!(h2, new_μgrid = h.μgrid, new_σgrid = h.σgrid, new_zgrid = h.zgrid)
			print_save("...")
			h.ϕa = h2.ϕa
			h.ϕb = h2.ϕb
			h.ϕa_ext = h2.ϕa_ext
			h.ϕb_ext = h2.ϕb_ext
			h.ϕc = h2.ϕc
			h.vf = h2.vf
			h.pngrid = h2.pngrid
			h.pN = h2.pN
			h.μ′ = h2.μ′
			h.σ′ = h2.σ′
			h.output = h2.output
			h.wage = h2.wage
			h.Ld = h2.Ld
			if !nodef
				h.repay = h2.repay
				h.welfare = h2.welfare
			else
				print_save(" Not loading default policy.")
			end
			print_save(" ✓")
		end
	catch
		print_save(": JLD file incompatible")
	end
	return h
end

function make_simulated_path(h::Hank, run_number, years=100)
	path, jz_series, Ndefs = simul(h; simul_length=4*(years+25), only_def_end=false)
	Tyears = floor(Int64,size(path.data, 1)*0.25)
	print_save("\n$Ndefs defaults in $Tyears years: default freq = $(floor(100*Ndefs/Tyears))%")
	trim_path!(path, 4*25)
	v_m = 0
	g = 0

	savedir = pwd() * "/../Output/run$(run_number)/"
	save(savedir * "path.jld", "path", path)
	
	pl = plot_simul(path)
	savejson(pl, savedir * "path.json")	
	try
		savefig(savedir * "path.pdf")
	catch
	end

	pl, πthres = plot_episodes(path; episode_type="onlyspread", slides=true, πthres=0.975)
	savejson(pl, savedir * "onlyspread_slides.json")
	try
		savefig(savedir * "onlyspread_slides.pdf")
	catch
	end

	# pl = plot_episodes(path; episode_type="onlyspread", slides=false, πthres=0.975)
	# savejson(pl, savedir * "onlyspread_paper.json")

	for resp in ["C"; "Y"]
		try
			make_IRF_plots(path; slides = true, create_plots = true, response = resp, savedir=savedir)
		catch
		end
	end

	params = pars(h)

	try
		v_m = simul_stats(path)
		targets = vec([0.96580506  0.01294576  0.96172496  0.01663608  0.96656486  0.10252351 64.57638889 23.48323041 15.94722222  6.08732167  56.4851069778397  94.479167])
		targetnames = ["AR(1) Output"; "σ(Output)"; "AR(1) Cons"; "σ(Cons) / σ(Output)"; "AR(1) Spreads"; "σ(spreads)"; "mean B/Y"; "std B/Y"; "mean unemp"; "std unemp"; "median Dom Holdings"; "mean wealth/Y" ]
		targets[4] = targets[4] / targets[2]
		v_m[4] = v_m[4] / v_m[2]
		print("\nTargets ✓")

		W = zeros(length(v_m),length(v_m))
		[W[jj,jj] = 1.0/targets[jj] for jj in 1:length(targets)]
		W[2,2] *= 100
		W[4,4] *= 50

		print("\nW ✓")
		g = (v_m - targets)'*W*(v_m-targets)
		print_save("\nObjective function = $(@sprintf("%0.3g",g))")
		print_save("\n")
		for jj in 1:length(targets)
			print_save("$(@sprintf("%0.3g",v_m[jj]))  ")
		end
		print_save("\nParams: ")
		for jj in 1:length(params)
			print_save("$(@sprintf("%0.3g",params[jj]))  ")
		end
		res = [targetnames v_m targets (targets-v_m)./targets]
		for jj in 1:size(res,1)
		    print_save("\n")
		    print_save(res[jj,1])
		    for ii in 2:size(res,2)
		    	print_save("     ")
		    	print_save("$(@sprintf("%0.3g",res[jj,ii]))")
		    end
		end
	catch
		print_save("\nWARNING: Found problems computing simulation statistics")
	end

	v_m[4] = v_m[4] * v_m[2]
	# save(savedir * "v_m.jld", "v_m", v_m)
	# write(savedir * "stats.txt", "$(v_m)")
	
	calib_table = make_calib_table(v_m)
	write(savedir * "calib_table.txt", calib_table)

	return g, path, πthres, v_m
end


function make_comparison_simul(h::Hank, noΔ, rep_agent, run_number, years, p_bench::Path, episode_type, πthres, savedir)

	mpe_iter!(h; nodef = true, noΔ = noΔ, rep_agent = rep_agent, run_number=run_number, maxiter = 21, save_copies=false)
	p_nodef, _, _ = simul(h; simul_length=4*(years+25), only_def_end=false)

	for (jj, slides) in enumerate([true; false])
		pcomp = plot_comparison_episodes(p_bench, p_nodef; episode_type = episode_type, slides = slides, πthres=πthres)
		savejson(pcomp, savedir * "comparison_crisis_nodef$(jj).json")
	end
	v_nodef = simul_stats(p_nodef)

	return v_nodef
end

pars(h::Hank) = [(1/h.β)^4-1; h.γ; h.τ; h.wbar; h.ρz; h.σz; h.tax; h.ρξ; h.σξ]



