using QuantEcon, BasisMatrices, Interpolations, Optim, LaTeXStrings, Distributions, JLD, HCubature, Distributed, Dates, ORCA, ForwardDiff, Printf, Random, LinearAlgebra, DataFrames, GLM

include("type_def.jl")
include("handle_itps.jl")
include("fiscal.jl")
include("hh_pb.jl")
include("comp_eqm.jl")
include("gov_pol.jl")
include("reporting_routines.jl")
include("simul.jl")



function solve_opt_value(sd::SOEdef)
	qᵍ_mat = reshape_long(sd, sd.eq[:qᵍ])
	ϕ = similar(qᵍ_mat)

	obj_f(x) = first(x)^2 - 2*first(x) + 1

	res = Optim.optimize(obj_f, -1, 1, GoldenSection())

	return size(ϕ)
end

