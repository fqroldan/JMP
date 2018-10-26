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

print_save("\nAggregate Demand around Debt Crises\n")

print_save("\nStarting $(location) run on $(nprocs()) cores at "*Dates.format(now(),"HH:MM"))
print_save("\nRun number: $run_number")

# Set options
local_run = true
nodef     = false
rep_agent = false

# Initialize type
function set_params(run_number)
	# 		β 	   tax	  RRA   τ
	xmin = [0.041; 0.02; 5.0;  0.1]
	xmax = [0.044; 0.03; 10.0; 0.2]
	N = length(xmin)
	s = SobolSeq(N)
	x = zeros(N)
	for j in 1:run_number
		x = next!(s, xmin, xmax)
	end
	return x
end
r_loc, tax, RRA, τ = set_params(run_number)

function make_guess(remote, local_run, nodef, rep_agent, r_loc, tax, RRA, τ)
	if remote || local_run
		h = Hank(; β=(1.0/(1.0+r_loc))^0.25, tax = tax, RRA=RRA, τ=τ, nodef = nodef, rep_agent = rep_agent
			# , Nω=2,Nϵ=3,Nb=2,Nμ=2,Nσ=2,Nξ=2,Nz=3
			);
		print_save("\nRun with r, tax, RRA, τ = $(round(r_loc, 3)), $(round(tax, 3)), $(round(RRA, 2)), $(round(τ, 2))")
		# h = load(pwd() * "/../../hank.jld", "h")
		try
			h2 = load(pwd() * "/../../hank.jld", "h")
			remote ? h2 = load(pwd() * "/../../hank.jld", "h") : h2 = load("hank.jld", "h")
			print_save("\nFound JLD file")
			if h.Ns == h2.Ns && h.Nω == h2.Nω && h.Nϵ == h2.Nϵ
				print_save(": loading previous results")
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
				# h.repay = h2.repay
				# h.welfare = h2.welfare
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
h = make_guess(remote, local_run, nodef, rep_agent, r_loc, tax, RRA, τ);

print_save("\nϵ: $(h.ϵgrid)")
print_save("\nz: $(h.zgrid)")
print_save("\nω: $(h.ωgrid)\n")

function make_simulated_path(h::Hank)
	path, jz_series = simul(h; simul_length=4*(1000+25), only_def_end=true)
	trim_path!(path, 4*25)
	save(pwd() * "/../../path.jld", "path", path)
	v_m = 0
	try
		v_m = simul_stats(path)
	catch
		print_save("\nWARNING: Found problems computing simulation statistics")
		v_m = 0
	end
	write(pwd()*"/stats.txt", "$(v_m)")
	Void
end

# Run
if remote || local_run
	# vfi!(h, verbose = true, remote = remote)
	mpe_iter!(h; remote = remote, nodef = nodef, rep_agent = rep_agent)
end

make_simulated_path(h)
Void
