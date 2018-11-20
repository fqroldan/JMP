using QuantEcon, BasisMatrices, Interpolations, Optim, NLopt, MINPACK, LaTeXStrings, Distributions, JLD, Sobol

function establish_run()
	location = "remote"
	run_number = 1
	try
		if pwd() == "/home/q/Dropbox/NYU/AToSR/Codes"
			location = "local"
		else
			for jj in 0:3
				if pwd()[end-jj-4] == '/' && pwd()[end-jj-3:end-jj-1] == "run"
					run_number = parse(Int64, pwd()[end-jj:end])
					break
				end
			end
		end
	catch
		print_save("\nWARNING: not able to assign run_number, reverting to default")
		run_number = 1
		location = "remote"
	end

	return location, run_number, (location=="remote")
end

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


location, run_number, remote = establish_run()
if location == "local"
	using Rsvg
end

print_save("\nAggregate Demand and Sovereign Debt Crises\n")

print_save("\nStarting $(location) run on $(nprocs()) cores at "*Dates.format(now(),"HH:MM"))
print_save("\nRun number: $run_number")

# Set options
local_run 	 = true
update_start = true
nodef     	 = false
rep_agent 	 = false

# Initialize type
function set_params(run_number, xcenter, xdist)
	N = length(xcenter)
	s = SobolSeq(N)
	x = zeros(N)
	for j in 1:run_number
		x = next!(s, xcenter-xdist, xcenter+xdist)
	end
	return x
end
#				 r_loc,   tax, RRA,    τ,  ρz,      σz,   ρξ,     σξ, w_bar
params_center = [0.055; 0.025; 7.5; 0.15; 0.9;   0.025; 0.95; 0.0025; 1.175]
xdist = 		[0.015; 0.01;  2.5; 0.05; 0.01; 0.0025; 0.01;  0.001; 0.025]

function find_new_cube(targets::Vector, W::Matrix; K::Int64=25)
	old_center = load(pwd() * "/../../../params_center.jld", "params_center")
	old_dist = load(pwd() * "/../../../xdist.jld", "xdist")

	k = 0
	curr_min = 1e10
	best_run = 0
	for jj in 1:K
		try
			v_m = load(pwd() * "/../../../v_m$(jj).jld", "v_m")

			gGMM = (v_m - targets)' * W * (v_m-targets)
			print_save("\ng = $gGMM on job $(jj)")
			if gGMM < curr_min
				curr_min = gGMM
				best_run = jj
			end
			k += 1
		end
	end
	
	if best_run == 0
		print_save("\nNo guesses available, keeping original parameters")
		new_dist = old_dist
		new_center = old_center
	else
		print_save("\nBest run from $k recovered trials = $best_run")
		if best_run == 1
			print_save(". Shrinking the cube.")
			new_dist = old_dist * 0.9
			new_center = old_center
		else
			print_save(". Computing new starting point.")
			s = SobolSeq(length(old_center))
			x = zeros(size(old_center))
			for j in 1:best_run
				x = next!(s, old_center-old_dist, old_center+old_dist)
			end
			η = 0.5
			new_center = x * η + old_center * (1.0 - η)
			new_dist = old_dist
		end
	end

	new_center += 0.05 * new_dist .* rand(size(new_center))

	return new_center, new_dist
end

if update_start
	targets = vec([0.96580506  0.01294576  0.96172496  0.01663608  0.96656486  0.10252351 64.57638889 23.48323041 15.94722222  6.08732167 60  1.42985986])

	W = zeros(length(targets),length(targets))
	[W[jj,jj] = 1.0/targets[jj] for jj in 1:length(targets)]

	params_center, xdist = find_new_cube(targets, W)
end

params = set_params(run_number, params_center, xdist)
r_loc, tax, RRA, τ, ρz, σz, ρξ, σξ = params

function make_guess(remote, local_run, nodef, rep_agent, r_loc, tax, RRA, τ, ρz, σz, ρξ, σξ)
	if remote || local_run
		h = Hank(; β=(1.0/(1.0+r_loc))^0.25, tax = tax, RRA=RRA, τ=τ, nodef = nodef, rep_agent = rep_agent, ρz=ρz, σz=σz, ρξ=ρξ, σξ=σξ
			# , Nω=2,Nϵ=3,Nb=2,Nμ=2,Nσ=2,Nξ=2,Nz=3
			);
		print_save("\nRun with r_loc, tax, RRA, τ, ρz, σz, ρξ, σξ = $(round(r_loc,3)), $(round(tax,3)), $(round(RRA,3)), $(round(τ,3)), $(round(ρz,3)), $(round(σz,3)), $(round(ρξ,3)), $(round(σξ,3))")
		# h = load(pwd() * "/../../hank.jld", "h")
		try
			h2 = load(pwd() * "/../../hank.jld", "h")
			remote ? h2 = load(pwd() * "/../../hank.jld", "h") : h2 = load("hank.jld", "h")
			print_save("\nFound JLD file")
			if h.Ns == h2.Ns && h.Nω == h2.Nω && h.Nϵ == h2.Nϵ
				print_save(": loading previous results")
				h2.μgrid = h.μgrid
				h2.σgrid = h.σgrid
				update_grids!(h2, new_μgrid = h.μgrid, new_σgrid = h.σgrid)
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
				h.repay = h2.repay
				h.welfare = h2.welfare
				print_save(" ✓")
			end
		catch
			print_save(": JLD file incompatible")
		end
	else
		print_save("\nLoading solved model file\n")
		h = load("../HPC_Output/hank.jld", "h")
	end
	return h
end
h = make_guess(remote, local_run, nodef, rep_agent, r_loc, tax, RRA, τ, ρz, σz, ρξ, σξ);

print_save("\nϵ: $(h.ϵgrid)")
print_save("\nz: $(h.zgrid)")
print_save("\nω: $(h.ωgrid)\n")

function make_simulated_path(h::Hank, run_number)
	path, jz_series, Ndefs = simul(h; simul_length=4*(1000+25), only_def_end=true)
	T = size(path.data, 1)
	print_save("\n$Ndefs defaults in $T years")
	trim_path!(path, 4*25)
	save(pwd() * "/../../path.jld", "path", path)
	v_m = 0
	try
		v_m = simul_stats(path)
		targets = vec([0.96580506  0.01294576  0.96172496  0.01663608  0.96656486  0.10252351 64.57638889 23.48323041 15.94722222  6.08732167 60  1.42985986])
		targetnames = ["AR(1) Output"; "σ(Output)"; "AR(1) Cons"; "σ(Cons)"; "AR(1) Spreads"; "σ(spreads)"; "mean B/Y"; "std B/Y"; "mean unemp"; "std unemp"; "mean Dom Holdings"; "std G/Y" ]

		W = zeros(length(v_m),length(v_m))
		[W[jj,jj] = 1.0/targets[jj] for jj in 1:length(targets)]
		if v_m[end-1] > targets[end-1]
			W[end-1,end-1] += 0.5
		end
		W[2,2] *= 10
		W[4,4] *= 10

		print_save("\nObjective function = $(@sprintf("%0.3g",(v_m - targets)'*W*(v_m-targets)))")
		res = [targetnames v_m targets (targets-v_m)./targets]
		for jj in 1:size(res,1)
		    print_save("\n")
		    for ii in 1:size(res,2)
		    	print_save("$(@sprintf("%0.3g",h.res[jj,ii]))")
		    	print_save("     ")
		    end
		end
	catch
		print_save("\nWARNING: Found problems computing simulation statistics")
	end
	save(pwd() * "/../../../v_m$(run_number).jld", "v_m", v_m)
	write(pwd()*"/stats.txt", "$(v_m)")
	nothing
end

# Run
if run_number == 1
	save(pwd() * "/../../../params_center.jld", "params_center", params_center)
	save(pwd() * "/../../../xdist.jld", "xdist", xdist)
end
if remote || local_run
	# vfi!(h, verbose = true, remote = remote)
	mpe_iter!(h; remote = remote, nodef = nodef, rep_agent = rep_agent, run_number=run_number)
end

make_simulated_path(h)

nothing