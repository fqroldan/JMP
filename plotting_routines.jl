using Interpolations

function plot_hh_policies(h::Hank; remote::Bool=false)
	leg = Array{LaTeXStrings.LaTeXString}(1, h.Nϵ)
	for jϵ in 1:h.Nϵ
		leg[jϵ] = latexstring("\\epsilon = $(round(h.ϵgrid[jϵ],2))")
	end

	knots = (h.ωgrid, h.ϵgrid, h.bgrid, h.μgrid, h.σgrid, h.wgrid, h.ζgrid, h.zgrid)
	itp_ϕa  = interpolate(knots, h.ϕa, Gridded(Linear()))
	itp_ϕb  = interpolate(knots, h.ϕb, Gridded(Linear()))
	itp_ϕc  = interpolate(knots, h.ϕc, Gridded(Linear()))
	itp_vf  = interpolate(knots, h.vf, Gridded(Linear()))

	show_b, show_μ, show_σ, show_w, show_ζ, show_z = mean(h.bgrid), mean(h.μgrid), mean(h.σgrid), mean(h.wgrid), h.ζgrid[1], h.zgrid[end]

	ϕc_mat = itp_ϕc[h.ωgrid_fine, h.ϵgrid, show_b, show_μ, show_σ, show_w, show_ζ, show_z]
	ϕa_mat = itp_ϕa[h.ωgrid_fine, h.ϵgrid, show_b, show_μ, show_σ, show_w, show_ζ, show_z]
	ϕb_mat = itp_ϕb[h.ωgrid_fine, h.ϵgrid, show_b, show_μ, show_σ, show_w, show_ζ, show_z]
	vf_mat = itp_vf[h.ωgrid_fine, h.ϵgrid, show_b, show_μ, show_σ, show_w, show_ζ, show_z]

	pc = plot(h.ωgrid_fine, ϕc_mat[:,:,1,1,1,1,1,1], title = "Consumption", label = "")
	pa = plot(h.ωgrid_fine, ϕa_mat[:,:,1,1,1,1,1,1], title = "Private Savings", label = "")
	pb = plot(h.ωgrid_fine, ϕb_mat[:,:,1,1,1,1,1,1], title = "Debt Purchases", label = leg)
	pv = plot(h.ωgrid_fine, vf_mat[:,:,1,1,1,1,1,1], title = "Value Function", label = "")

	l = @layout([a b; c d])

	plot(pc, pv, pa, pb, layout=l, lw = 2, xlabel = L"\omega_t", size = (600,700))
	plot!(h.ωgrid_fine, ones(h.ωgrid_fine), lw = 1, lc = "black", ls = :dashdot, label = "")
	plot!(bg_outside = RGBA(0.99,0.99,0.99, 0.75))
	plot!(titlefont=font(11), guidefont=font(8), tickfont=font(7), titlefont=font(12))
	# plot!(titlefont=font(11,"Palatino"), guidefont=font(8,"Palatino"), tickfont=font(7,"Palatino"), titlefont=font(12,"Palatino"))
	if remote
		path = pwd() * "/../../Graphs/"
	else
		path = pwd() * "/../Graphs/"
	end
	# savefig(path * "hh2.pdf")
	if remote
		path = pwd() * "/../../Graphs/"
	else
		path = pwd() * "/../Graphs/"
	end
	savefig(path * "hh2.png")

	show_b, show_μ, show_σ, show_w, show_ζ, show_z = mean(h.bgrid), mean(h.μgrid), mean(h.σgrid), mean(h.wgrid), h.ζgrid[2], h.zgrid[1]

	ϕc_mat = itp_ϕc[h.ωgrid_fine, h.ϵgrid, show_b, show_μ, show_σ, show_w, show_ζ, show_z]
	ϕa_mat = itp_ϕa[h.ωgrid_fine, h.ϵgrid, show_b, show_μ, show_σ, show_w, show_ζ, show_z]
	ϕb_mat = itp_ϕb[h.ωgrid_fine, h.ϵgrid, show_b, show_μ, show_σ, show_w, show_ζ, show_z]
	vf_mat = itp_vf[h.ωgrid_fine, h.ϵgrid, show_b, show_μ, show_σ, show_w, show_ζ, show_z]

	pc = plot(h.ωgrid_fine, ϕc_mat[:,:,1,1,1,1,1,1], title = "Consumption", label = "")
	pa = plot(h.ωgrid_fine, ϕa_mat[:,:,1,1,1,1,1,1], title = "Private Savings", label = "")
	pb = plot(h.ωgrid_fine, ϕb_mat[:,:,1,1,1,1,1,1], title = "Debt Purchases", label = leg)
	pv = plot(h.ωgrid_fine, vf_mat[:,:,1,1,1,1,1,1], title = "Value Function", label = "")

	l = @layout([a b; c d])

	plot(pc, pv, pa, pb, layout=l, lw = 2, xlabel = L"\omega_t", size = (600,700))
	plot!(h.ωgrid_fine, ones(h.ωgrid_fine), lw = 1, lc = "black", ls = :dashdot, label = "")
	plot!(bg_outside = RGBA(0.99,0.99,0.99, 0.85))
	plot!(titlefont=font(11,"Palatino"), guidefont=font(8,"Palatino"), tickfont=font(7,"Palatino"), titlefont=font(12,"Palatino"))
	if remote
		path = pwd() * "/../../Graphs/"
	else
		path = pwd() * "/../Graphs/"
	end
	savefig(path * "hh3.png")


#=	leg = Array{LaTeXStrings.LaTeXString}(1, h.Nω)
	for jω in 1:h.Nω
		leg[jω] = latexstring("ω = $(round(h.ωgrid[jω],2))")
	end=#

	leg = latexstring("\\omega = $(round(mean(h.ωgrid),2))")

	ϕc_mat = itp_ϕc[mean(h.ωgrid), mean(h.ϵgrid), show_b, show_μ, show_σ, show_w, show_ζ, h.zgrid]
	ϕa_mat = itp_ϕa[mean(h.ωgrid), mean(h.ϵgrid), show_b, show_μ, show_σ, show_w, show_ζ, h.zgrid]
	ϕb_mat = itp_ϕb[mean(h.ωgrid), mean(h.ϵgrid), show_b, show_μ, show_σ, show_w, show_ζ, h.zgrid]
	vf_mat = itp_vf[mean(h.ωgrid), mean(h.ϵgrid), show_b, show_μ, show_σ, show_w, show_ζ, h.zgrid]

	pc = plot(h.zgrid, ϕc_mat[:,1,1,1,1,1,1,:]', title = "Consumption", label = "")
	pa = plot(h.zgrid, ϕa_mat[:,1,1,1,1,1,1,:]', title = "Private Savings", label = "")
	pb = plot(h.zgrid, ϕb_mat[:,1,1,1,1,1,1,:]', title = "Debt Purchases", label = leg)
	pv = plot(h.zgrid, vf_mat[:,1,1,1,1,1,1,:]', title = "Value Function", label = "")

	l = @layout([a b; c d])

	plot(pc, pv, pa, pb, layout=l, lw = 2, xlabel = L"z_t", size = (600,700))
	plot!(bg_outside = RGBA(0.99,0.99,0.99, 0.85))
	plot!(titlefont=font(11,"Palatino"), guidefont=font(8,"Palatino"), tickfont=font(7,"Palatino"), titlefont=font(12,"Palatino"))
	if remote
		path = pwd() * "/../../Graphs/"
	else
		path = pwd() * "/../Graphs/"
	end
	savefig(path * "hh_z.png")



	return Void
end

function plot_state_funcs(h::Hank; remote::Bool=false)
	R, T, ℓ, Rep, qˢ, qᵇ, Π = _unpackstatefs(h)

	A⁺_mat	= reshape(h.A⁺, h.Nb, h.Nμ, h.Nσ, h.Nz)
	A⁻_mat	= reshape(h.A⁻, h.Nb, h.Nμ, h.Nσ, h.Nz)

	P = kron(h.Ps, kron(h.Pϵ, speye(h.Na)))

	EβR = (P * (R ./ Π)) ./ qˢ
	EβR = h.β * reshape(EβR, h.Na, h.Nϵ, h.Nb, h.Nμ, h.Nσ, h.Nz)[1,1,:,:,:,:]

	R	= reshape(R, 		h.Na, h.Nϵ, h.Nb, h.Nμ, h.Nσ, h.Nz)[1,1,:,:,:,:]
	T	= reshape(T, 		h.Na, h.Nϵ, h.Nb, h.Nμ, h.Nσ, h.Nz)[1,1,:,:,:,:]
	qˢ	= reshape(qˢ, 		h.Na, h.Nϵ, h.Nb, h.Nμ, h.Nσ, h.Nz)[1,1,:,:,:,:]
	qᵇ	= reshape(qᵇ, 		h.Na, h.Nϵ, h.Nb, h.Nμ, h.Nσ, h.Nz)[1,1,:,:,:,:]
	Π	= reshape(Π, 		h.Na, h.Nϵ, h.Nb, h.Nμ, h.Nσ, h.Nz)[1,1,:,:,:,:]
	w 	= reshape(h.wage, 	h.Na, h.Nϵ, h.Nb, h.Nμ, h.Nσ, h.Nz)[1,1,:,:,:,:]
	qᵍ 	= reshape(h.qᵍ, 	h.Na, h.Nϵ, h.Nb, h.Nμ, h.Nσ, h.Nz)[1,1,:,:,:,:]

	Πdev = (Π)/h.Πstar

	Z = zeros(h.Nb, h.Nμ, h.Nσ, h.Nz)
	for (jz, zv) in enumerate(h.zgrid)
		Z[:,:,:,jz] = zv
	end	
	L = (w * (1-h.τ)/h.θ * h.Ξ).^(1/h.χ)
	Y = Z .* L

	ψ = Y .* (1 - w./Z - 0.5*h.η*(Π/h.Πstar - 1).^2)

	Π = (Π.^4 - 1)*100

	i = ((1./qᵇ).^4 - 1) * 100 # Annualized percent nominal rate
	iˢ = ((1./qˢ).^4 - 1) * 100
	
	j = 3

	# l = @layout([a b; c d; e f; g h i])

	pEβR = plot(h.bgrid, vec(EβR[:,j,j,j]), title=L"\beta E[R^S]", label = "")
	pqᵍ	 = plot(h.bgrid, vec(qᵍ[:,j,j,j]), 	title=L"q^g", label = "")
	pΠ	 = plot(h.bgrid, vec(Π[:,j,j,j]), 	title=L"Π", label = "")
	pAp	 = plot(h.bgrid, vec(A⁺_mat[:,j,j,j]), title=L"A^+", label = "")
	pAm	 = plot(h.bgrid, vec(A⁻_mat[:,j,j,j]), title=L"A^-", label = "")
	p_i	 = plot(h.bgrid, [vec(iˢ[:,j,j,j]) vec(i[:,j,j,j])], title=L"i", label = [L"i^s" L"i^b"])
	pT	 = plot(h.bgrid, vec(T[:,j,j,j]), 	title=L"T", label = "")
	pL	 = plot(h.bgrid, [vec(L[:,j,j,j]) vec(Πdev[:,j,j,j])], title=L"L", label = [L"L" L"\tilde{\pi}"])
	pw 	 = plot(h.bgrid, vec(w[:,j,j,j]), title=L"w", label = "")
	pψ = plot(h.bgrid, vec(ψ[:,j,j,j]), title=L"\psi", label = "")
	
	plot(pEβR, pqᵍ, pΠ, pAp, pAm, p_i, pT, pL, pψ, xlabel = L"B_t", layout = (3,3), lw = 1.5)
	if remote
		path = pwd() * "/../../Graphs/"
	else
		path = pwd() * "/../Graphs/"
	end
	savefig(path * "fs_b.png")

	# l = @layout([a b; c d; e f g])
	pEβR = plot(h.μgrid, vec(EβR[j,:,j,j]), title=L"\beta E[R^S]", label = "")
	pqᵍ	 = plot(h.μgrid, vec(qᵍ[j,:,j,j]), title=L"q^g", label = "")
	pΠ	 = plot(h.μgrid, vec(Π[j,:,j,j]), title=L"Π", label = "")
	pAp	 = plot(h.μgrid, vec(A⁺_mat[j,:,j,j]), title=L"A^+", label = "")
	pAm	 = plot(h.μgrid, vec(A⁻_mat[j,:,j,j]), title=L"A^-", label = "")
	p_i	 = plot(h.μgrid, [vec(iˢ[j,:,j,j]) vec(i[j,:,j,j])], title=L"i", label = [L"i^s" L"i^b"])
	pT	 = plot(h.μgrid, vec(T[j,:,j,j]), title=L"T", label = "")
	pL	 = plot(h.μgrid, [vec(L[j,:,j,j]) vec(Πdev[j,:,j,j])], title=L"L", label = [L"L" L"\tilde{\pi}"])
	pw = plot(h.μgrid, vec(w[j,:,j,j]), title=L"w", label = "")
	pψ = plot(h.μgrid, vec(ψ[j,:,j,j]), title=L"\psi", label = "")
	
	plot(pEβR, pqᵍ, pΠ, pAp, pAm, p_i, pT, pL, pψ, xlabel = L"\mu_t", layout = (3,3), lw = 1.5)
	if remote
		path = pwd() * "/../../Graphs/"
	else
		path = pwd() * "/../Graphs/"
	end
	savefig(path * "fs_mu.png")

	# l = @layout([a b; c d; e f g])
	pEβR = plot(h.σgrid, vec(EβR[j,j,:,j]), title=L"\beta E[R^S]", label = "")
	pqᵍ	 = plot(h.σgrid, vec(qᵍ[j,j,:,j]), title=L"q^g", label = "")
	pΠ	 = plot(h.σgrid, vec(Π[j,j,:,j]), title=L"Π", label = "")
	pAp	 = plot(h.σgrid, vec(A⁺_mat[j,j,:,j]), title=L"A^+", label = "")
	pAm	 = plot(h.σgrid, vec(A⁻_mat[j,j,:,j]), title=L"A^-", label = "")
	p_i	 = plot(h.σgrid, [vec(iˢ[j,j,:,j]) vec(i[j,j,:,j])], title=L"i", label = [L"i^s" L"i^b"])
	pT	 = plot(h.σgrid, vec(T[j,j,:,j]), title=L"T", label = "")
	pL	 = plot(h.σgrid, [vec(L[j,j,:,j]) vec(Πdev[j,j,:,j])], title=L"L", label = [L"L" L"\tilde{\pi}"])
	pw = plot(h.σgrid, vec(w[j,j,:,j]), title=L"w", label = "")
	pψ = plot(h.σgrid, vec(ψ[j,j,:,j]), title=L"\psi", label = "")
	
	plot(pEβR, pqᵍ, pΠ, pAp, pAm, p_i, pT, pL, pψ, xlabel = L"\sigma_t", layout = (3,3), lw = 1.5)
	if remote
		path = pwd() * "/../../Graphs/"
	else
		path = pwd() * "/../Graphs/"
	end
	savefig(path * "fs_sigma.png")

	# l = @layout([a b; c d; e f g])
	pEβR = plot(h.zgrid, vec(EβR[j,j,j,:]), title=L"\beta E[R^S]", label = "")
	pqᵍ	 = plot(h.zgrid, vec(qᵍ[j,j,j,:]), title=L"q^g", label = "")
	pΠ	 = plot(h.zgrid, vec(Π[j,j,j,:]), title=L"Π", label = "")
	pAp	 = plot(h.zgrid, vec(A⁺_mat[j,j,j,:]), title=L"A^+", label = "")
	pAm	 = plot(h.zgrid, vec(A⁻_mat[j,j,j,:]), title=L"A^-", label = "")
	p_i	 = plot(h.zgrid, [vec(iˢ[j,j,j,:]) vec(i[j,j,j,:])], title=L"i", label = [L"i^s" L"i^b"])
	pT	 = plot(h.zgrid, vec(T[j,j,j,:]), title=L"T", label = "")
	pL	 = plot(h.zgrid, [vec(L[j,j,j,:]) vec(Πdev[j,j,j,:])], title=L"L", label = [L"L" L"\tilde{\pi}"])
	pw = plot(h.zgrid, vec(w[j,j,j,:]), title=L"w", label = "")
	pψ = plot(h.zgrid, vec(ψ[j,j,j,:]), title=L"\psi", label = "")
	
	plot(pEβR, pqᵍ, pΠ, pAp, pAm, p_i, pT, pL, pψ, xlabel = L"z_t", layout = (3,3), lw = 1.5)
	if remote
		path = pwd() * "/../../Graphs/"
	else
		path = pwd() * "/../Graphs/"
	end
	savefig(path * "fs_z.png")

	Void
end

function plot_LoM(h::Hank, μ′, σ′; remote::Bool=false)

	μ′ = reshape(μ′, h.Nb, h.Nμ, h.Nσ, h.Nz)
	σ′ = reshape(σ′, h.Nb, h.Nμ, h.Nσ, h.Nz)

	j 	= ceil(Int, length(h.bgrid)/2)
	pμb = plot(h.bgrid, vec(μ′[:,j,j,j]), title=L"\mu", xlabel = L"B_t", label = "")
	pσb = plot(h.bgrid, vec(σ′[:,j,j,j]), title=L"\sigma", xlabel = L"B_t", label = "")

	j 	= ceil(Int, length(h.μgrid)/2)
	pμμ = plot(h.μgrid, vec(μ′[j,:,j,j]), title=L"\mu", xlabel = L"\mu_t", label = "")
	pσμ = plot(h.μgrid, vec(σ′[j,:,j,j]), title=L"\sigma", xlabel = L"\mu_t", label = "")

	j 	= ceil(Int, length(h.σgrid)/2)
	pμσ = plot(h.σgrid, vec(μ′[j,j,:,j]), title=L"\mu", xlabel = L"\sigma_t", label = "")
	pσσ = plot(h.σgrid, vec(σ′[j,j,:,j]), title=L"\sigma", xlabel = L"\sigma_t", label = "")

	j 	= ceil(Int, length(h.zgrid)/2)
	pμz = plot(h.zgrid, vec(μ′[j,j,j,:]), title=L"\mu", xlabel = L"z_t", label = "")
	pσz = plot(h.zgrid, vec(σ′[j,j,j,:]), title=L"\sigma", xlabel = L"z_t", label = "")
	
	plot(pμb, pσb, pμμ, pσμ, pμσ, pσσ, pμz, pσz, layout = (4,2), lw = 1.5, size = (600,800))
	if remote
		path = pwd() * "/../../Graphs/"
	else
		path = pwd() * "/../Graphs/"
	end
	savefig(path * "LoMs.png")
	Void
end

function labor_demand(h::Hank, w, z, pN; get_both::Bool = false)

	Ld_nontradables = (h.α_N * pN ./ w).^(1.0/(1.0-h.α_N))
	Ld_tradables    = (h.α_T * z  ./ w).^(1.0/(1.0-h.α_T))

	if get_both
		return Ld_nontradables, Ld_tradables
	else
		return Ld_nontradables + Ld_tradables
	end
end

function plot_labor_demand(h::Hank; remote::Bool=false)
	jz = floor(Int, h.Nz/2)
	z_show = h.zgrid[jz]

	plot(title = "Labor demand", xlabel = L"w_t")
	plot!(h.wgrid, ones(h.wgrid), lw = 1, ls=:dashdot, lc=:black, label = "")
	for (jpN, pNv) in enumerate(h.pngrid)
		Ld = labor_demand(h, h.wgrid, z_show, pNv)

		label = "p_N = $(round(pNv,2))"
		# label = latexstring("p_N = $(round(pNv,2))")
		plot!(h.wgrid, Ld, label = label)
	end
	plot!(titlefont=font(11,"Palatino"), guidefont=font(8,"Palatino"), tickfont=font(7,"Palatino"), titlefont=font(12,"Palatino"))

	if remote
		path = pwd() * "/../../Graphs/"
	else
		path = pwd() * "/../Graphs/"
	end
	savefig(path * "labordemand.png")
	Void
end

function plot_nontradables_demand(h::Hank; remote::Bool=false)
	jb = floor(Int, h.Nb/2)
	jμ = floor(Int, h.Nμ/2)
	jσ = floor(Int, h.Nσ/2)
	jw = floor(Int, h.Nw/2)
	jζ = floor(Int, h.Nζ/2)
	jz = floor(Int, h.Nz/2)
	
	bv, μv, σv, wv, ζv, zv = h.bgrid[jb], h.μgrid[jμ], h.σgrid[jσ], h.wgrid[jw], h.ζgrid[jζ], h.zgrid[jz]

	if remote
		path = pwd() * "/../../Graphs/"
	else
		path = pwd() * "/../Graphs/"
	end
	savefig(path * "N_demand.png")
	Void
end