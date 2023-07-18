include("main_serial.jl")


function lower_defcost!(sd::SOEdef)

    @assert sd.pars[:Δ] == 0.1

    sd.pars[:Δ] = 0.0

    comp_eqm!(sd, tol = 5e-4)
end

function more_progressive!(sd::SOEdef)

    @assert sd.pars[:τ] == 0.31

    sd.pars[:τ] = 0.32

    comp_eqm!(sd, tol = 5e-4)
end

function more_risk!(sd::SOEdef)
    
    sd.pars[:σϵ] = 0.133

    update_prob_ϵ!(sd)

    comp_eqm!(sd, tol = 5e-4)
end