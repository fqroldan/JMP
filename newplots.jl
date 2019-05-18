using PlotlyJS, Interpolations, JLD, Distributions, HCubature, ORCA, KernelDensity, Printf

include("type_def.jl")
# include("interp_atosr.jl")
include("plotting_routines.jl")
include("simul.jl")
include("reporting_routines.jl")

function maketable_bench(; fecha::String="mar20")
	p = load("../HPC_Output/" * fecha * "/path_bench.jld", "path")
	v_m = simul_stats(p)
	targets = vec([0.96580506  0.01294576  0.96172496  0.01663608  0.96656486  0.10252351 64.57638889 23.48323041 15.94722222  6.08732167  56.4851069778397  94.479167])
	targetnames = ["AR(1) Output"; "σ(Output)"; "AR(1) Cons"; "σ(C)/σ(Y)"; "AR(1) Spreads"; "σ(spreads)"; "mean B/Y"; "std B/Y"; "mean unemp"; "std unemp"; "mean Dom Holdings"; "mean wealth/Y" ]
	targets[4] = targets[4] / targets[2]

	W = zeros(length(v_m),length(v_m))
	[W[jj,jj] = 1.0/targets[jj] for jj in 1:length(targets)]
	W[2,2] *= 100
	W[4,4] *= 50
	v_m[4] = v_m[4] / v_m[2]


	print("\nObjective function = $(@sprintf("%0.3g",(v_m - targets)'*W*(v_m-targets)))")
	print("\n")
	for jj in 1:length(targets)
		print("$(@sprintf("%0.3g",v_m[jj]))  ")
	end
	res = [targetnames v_m targets (targets-v_m)./targets]
	for jj in 1:size(res,1)
		print("\n")
		print(res[jj,1])
		for ii in 2:size(res,2)
			print("     ")
			print("$(@sprintf("%0.3g",res[jj,ii]))")
		end
	end

	ζvec = series(p, :ζ)

	def = [(ζvec[jj] == 2) .& (ζvec[jj-1] == 1) for jj in 2:length(ζvec)]

	print("\n"); print(sum(def)/4000)

	v_m[4] *= v_m[2]
	make_calib_table(v_m)

	nothing
end

function maketable_params(; fecha::String="mar20")
	h = load(pwd() * "/../HPC_Output/" * fecha * "/hank_bench.jld", "h")

	params = [(1/h.β)^4-1; h.γ; h.τ; h.wbar; h.ρz; h.σz; h.tax; h.ρξ; h.σξ]
	make_params_table(collect(params))

	nothing
end

function makeplot_episodes(; fecha::String="mar20", ext::Vector=[], slides::Bool=false, πthres::Float64=0.5, episode_type::String="onlyspread")

	path = load("../HPC_Output/" * fecha * "/path_bench.jld", "path")

	p1 = plot_episodes(path; episode_type = episode_type, slides = slides, πthres=πthres)
	relayout!(p1, height = 800)
	if slides
		name = "_slides"
		relayout!(p1, width=1450, font_family="Fira Sans Light", font_size=14)
	else
		name = "_paper"
		relayout!(p1, font_family = "STIX Two Text", font_size=14, width = 1200, height = 700)
	end
	
	if length(ext) > 0
		for (jj, jext) in enumerate(ext)
			savefig(p1, pwd()*"/../Graphs/" * episode_type * name * "." * jext)
		end
	end

	return p1
end

function maketable_comparison(p_bench, p_nodef, p_nodelta=p_nodef)
	v_m = simul_stats(p_bench, nodef=false)
	ζ_vec = series(p_bench, :ζ)
	v_m_nodef = simul_stats(p_nodef; ζ_vec = ζ_vec, nodef=false)
	if p_nodelta == p_nodef
		v_m_nodelta = []
		print("\nOnly No-def")
	else
		v_m_nodelta = simul_stats(p_nodelta; ζ_vec = ζ_vec, nodef=false)
		print("\nBench, no def, and no delta")
	end
	targets = vec([0.96580506  0.01294576  0.96172496  0.01663608  0.96656486  0.10252351 64.57638889 23.48323041 15.94722222  6.08732167  56.4851069778397  94.479167])
	targetnames = ["AR(1) Output"; "σ(Output)"; "AR(1) Cons"; "σ(C)/σ(Y)"; "AR(1) Spreads"; "σ(spreads)"; "mean B/Y"; "std B/Y"; "mean unemp"; "std unemp"; "mean Dom Holdings"; "mean wealth/Y" ]

	print("\n\nComparison table\n")
	make_calib_table_comp(v_m, v_m_nodef, v_m_nodelta)
end

function makeplot_comparison(; fecha::String="mar20", ext::Vector=[], slides::Bool=false, πthres::Float64=0.5, levels::Bool=true, no_nodelta=true)
	p_bench = load("../HPC_Output/" * fecha * "/path_bench.jld", "path");
	p_nodef = load("../HPC_Output/" * fecha * "/path_nodef.jld", "path");
	if !no_nodelta
		p_nodelta = load("../HPC_Output/" * fecha * "/path_nodelta.jld", "path");
	else
		p_nodelta = p_nodef
	end

	p1 = plot_comparison_episodes(p_bench, p_nodef, p_nodelta; episode_type = "onlyspread", slides = slides, πthres=πthres, levels=levels)

	if slides
		name = "_slides"
		relayout!(p1, width=1450, height = 700, font_family="Fira Sans Light", font_size=14)
	else
		name = "_paper"
		relayout!(p1, font_family = "STIX Two Text", font_size=14, width = 1200, height = 700)
	end

	if levels
		name *= "_levels"
	else
		name *= "_changes"
	end

	if length(ext) > 0
		for (jj, jext) in enumerate(ext)
			savefig(p1, pwd() * "/../Graphs/comparison_crisis_nodef" * name * "." * jext)
		end
	end

	maketable_comparison(p_bench, p_nodef, p_nodelta)

	p1
end


function make_IRF_plots(; fecha::String="mar20", slides::Bool=true, response::String="Y", impulse::String="z", verbose::Bool=false, create_plots::Bool=false, cond::String="bench")
	
	p = load("../HPC_Output/" * fecha * "/path_"*cond*".jld", "path");

	zvec = series(p, :z)
	ξvec = series(p, :ξ)
	Wvec = series(p, :Wr) - series(p, :Wd)
	Cvec = log.(series(p, :C))
	Yvec = log.(series(p, :Y))
	πvec = series(p, :π)*100
	Bvec = series(p, :B)./series(p, :Y)
	t0 = 11
	T = length(zvec)

	if impulse == "z"
		shock_vec = zvec
		ρ_shock = 0.9
	elseif impulse == "ξ"
		shock_vec = ξvec
		ρ_shock = 0.95
	else
		throw(error("Keyword impulse has to be 'z' or 'ξ'"))
	end
	
	if response == "Y"
		yvec = Yvec
	elseif response == "C"
		yvec = Cvec
	elseif response == "π"
		yvec = πvec
	end

	E_shock = shock_vec * ρ_shock
	y = yvec[t0:T]
	X = Bvec[t0:T]

	eps_z = [shock_vec[jt] - E_shock[jt-1] for jt in t0:T]

	Bmed = quantile(Bvec, 0.75)
	zlow = quantile(eps_z, 0.25)
	# eps_z = min.(eps_z, zlow)
	# eps_z = max.(eps_z, quantile(eps_z, 0.9))

	p1 = plot(
		scatter(;x=1:T, y=eps_z)
		, Layout(;shapes=[hline(minimum(eps_z)); hline(quantile(eps_z,0.01))])
	)

	H = 20
	β = Matrix{Float64}(undef, H+1,3)
	βhB = Matrix{Float64}(undef, H+1,6)
	βlz = Matrix{Float64}(undef, H+1,6)
	βboth = Matrix{Float64}(undef, H+1,3)
	for jh in 1:H+1
		yh = y[jh:end]
		Bh = X[1:end-jh+1]
		e_z = eps_z[1:end-jh+1]

		dummy_highB = (Bh.>Bmed)
		dummy_lowz = (e_z.<zlow)
		data = DataFrame(yh=yh, eps_z=e_z, X = Bh, ind_B = dummy_highB, ind_z = dummy_lowz)

		#data = data[data[:eps_z].<= zlow,:]
		# println(size(data[(data[:eps_z].<=zlow) .* (data[:X].>=Bmed),:]))
		OLS = glm(@formula(yh ~ eps_z), data, Normal(), IdentityLink())
		if verbose
			println(OLS)
		end
		β[jh,1] = coef(OLS)[2]
		β[jh,2] = coef(OLS)[2] + stderror(OLS)[2]*1.96
		β[jh,3] = coef(OLS)[2] - stderror(OLS)[2]*1.96
		#data_cond = data[(data[:X].>=Bmed),:]
		#OLS = glm(@formula(yh ~ eps_z), data_cond, Normal(), IdentityLink())    
		OLS = glm(@formula(yh ~ eps_z + eps_z&ind_B), data, Normal(), IdentityLink())
		if verbose
			println(OLS)
		end
		βhB[jh,1] = coef(OLS)[2]
		βhB[jh,2] = coef(OLS)[2] + stderror(OLS)[2]*1.96
		βhB[jh,3] = coef(OLS)[2] - stderror(OLS)[2]*1.96
		βhB[jh,4] = coef(OLS)[2] + coef(OLS)[3]
		βhB[jh,5] = coef(OLS)[2] + coef(OLS)[3] + stderror(OLS)[3]*1.96 + stderror(OLS)[2]*1.96
		βhB[jh,6] = coef(OLS)[2] + coef(OLS)[3] - stderror(OLS)[3]*1.96 - stderror(OLS)[2]*1.96
		#data_cond = data[(data[:X].<=Bmed),:]
		#OLS = glm(@formula(yh ~ eps_z), data_cond, Normal(), IdentityLink())    
		OLS = glm(@formula(yh ~ eps_z + eps_z&ind_z), data, Normal(), IdentityLink())
		if verbose
			println(OLS)
		end
		βlz[jh,1] = coef(OLS)[2]
		βlz[jh,2] = coef(OLS)[2] + stderror(OLS)[2]*1.96
		βlz[jh,3] = coef(OLS)[2] - stderror(OLS)[2]*1.96
		βlz[jh,4] = coef(OLS)[2] + coef(OLS)[3]
		βlz[jh,5] = coef(OLS)[2] + coef(OLS)[3] + stderror(OLS)[3]*1.96 #+ stderror(OLS)[2]*1.96
		βlz[jh,6] = coef(OLS)[2] + coef(OLS)[3] - stderror(OLS)[3]*1.96 #- stderror(OLS)[2]*1.96
		if verbose
			println(coef(OLS))
			println(βlz[jh,:])
		end
		#data_cond = data[(data[:eps_z].<zlow) .* (data[:X].>=Bmed),:]
		#OLS = glm(@formula(yh ~ eps_z), data_cond, Normal(), IdentityLink())    
		OLS = glm(@formula(yh ~ eps_z + ind_B + eps_z*ind_B + ind_z + eps_z*ind_z + eps_z*ind_B*ind_z), data, Normal(), IdentityLink())
		if verbose
			println(OLS)
		end
		βboth[jh,1] = coef(OLS)[2]
		βboth[jh,2] = coef(OLS)[2] + stderror(OLS)[2]*1.96
		βboth[jh,3] = coef(OLS)[2] - stderror(OLS)[2]*1.96
	end
	#β_h *= sqrt(var(zvec))
	
	yaxistitle = "∂log " * response * " / ∂log z"
	
	pYz = plot([
			scatter(;x=0:H, y=β[:,1], line_color=col[1], name="βₕ")        
			scatter(;x=0:H, y=β[:,3], line_width=0, showlegend=false, name="βl")
			scatter(;x=0:H, y=β[:,2], fill="tonexty", line_color="#bbd6e8", line_width=0, showlegend=false, name="βh")
		], Layout(;xaxis_zeroline=false, yaxis_zeroline=false, xaxis_title="Quarters", yaxis_title=yaxistitle,
			legend_orientation="h", legend_x = 0.1, width = 600, height = 250, font_family = "STIX Two Text",
			shapes=[hline(0, line_dash="dot", line_width=1)]))

	pYzz = plot([
			scatter(;x=0:H, y=βlz[:,3], line_width=0, showlegend=false, name="βl")
			scatter(;x=0:H, y=βlz[:,2], fill="tonexty", line_color="#bbd6e8", line_width=0, showlegend=false, name="βh")
			scatter(;x=0:H, y=βlz[:,6], line_width=0, showlegend=false, name="βl")
			scatter(;x=0:H, y=βlz[:,5], fill="tonexty", line_color="#bfe2bf", line_width=0, showlegend=false, name="βh")
			scatter(;x=0:H, y=βlz[:,4], line_color=col[3], name="βₕ (low z)", line_dash="dashdot")        
			scatter(;x=0:H, y=βlz[:,1], line_color=col[1], name="βₕ")        
		], Layout(;xaxis_zeroline=false, yaxis_zeroline=false, xaxis_title="Quarters", yaxis_title=yaxistitle,
			legend=attr(;orientation="h", x = 0.1, traceorder="reversed"), width = 600, height = 400, font_family = "STIX Two Text",
			shapes=[hline(0, line_dash="dot", line_width=1)]))

	pYzB = plot([
			scatter(;x=0:H, y=β[:,3], line_width=0, showlegend=false, name="βl")
			scatter(;x=0:H, y=β[:,2], fill="tonexty", line_color="#bbd6e8", line_width=0, showlegend=false, name="βh")
			scatter(;x=0:H, y=βhB[:,6], line_width=0, showlegend=false, name="βl")
			scatter(;x=0:H, y=βhB[:,5], fill="tonexty", line_color="#bfe2bf", line_width=0, showlegend=false, name="βh")
			scatter(;x=0:H, y=βhB[:,4], line_color=col[3], name="βₕ (high B)", line_dash="dashdot")        
			scatter(;x=0:H, y=β[:,1], line_color=col[1], name="βₕ")        
		], Layout(;xaxis_zeroline=false, yaxis_zeroline=false, xaxis_title="Quarters", yaxis_title=yaxistitle,
			legend=attr(;orientation="h", x = 0.1, traceorder="reversed"), width = 600, height = 400, font_family = "STIX Two Text",
			shapes=[hline(0, line_dash="dot", line_width=1)]))
	
	if slides
		relayout!(pYz, font_size = 16, plot_bgcolor = "rgb(250,250,250)", paper_bgcolor = "rgb(250,250,250)", font_family = "Fira Sans Light")
		relayout!(pYzz, font_size = 16, plot_bgcolor = "rgb(250,250,250)", paper_bgcolor = "rgb(250,250,250)", font_family = "Fira Sans Light")
		relayout!(pYzB, font_size = 16, plot_bgcolor = "rgb(250,250,250)", paper_bgcolor = "rgb(250,250,250)", font_family = "Fira Sans Light")
	end

	relayout!(pYzz, yaxis_range = [-0.1; 1.5])
	relayout!(pYzB, yaxis_range = [-0.1; 1.5])

	slides ? name = "_slides" : name = ""

	if create_plots
		savefig(pYz, pwd()*"/../Graphs/elast_" * response * "z" * name * "_" * cond * ".pdf")
		savefig(pYzz, pwd()*"/../Graphs/elast_" * response * "zz" * name * "_" * cond * ".pdf")
		savefig(pYzB, pwd()*"/../Graphs/elast_" * response * "zB" * name * "_" * cond * ".pdf")
	end
	[pYzz pYzB]
end


function makeplot_simulpath(; fecha="mar20", slides::Bool=true)
	cond="bench"
	path_entry = load("../HPC_Output/" * fecha * "/path_"*cond*".jld", "path");

	Tmin = 0*size(path_entry.data, 1) #- 800
	p1 = plot_simul(path_entry; remote=false, trim=Tmin)

	# B_vec = series(path_entry,:NX)
	# p1 = plot(B_vec)

	return p1
end


print("maketable_bench, maketable_params, makeplot_episodes, maketable_comparison, makeplot_comparison, make_IRF_plots, makeplot_simulpath")

print("\nDates: mar20, apr28, apr29, may4, may6, may10, may12, may15")