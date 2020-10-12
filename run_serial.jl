include("main_serial.jl")

sd = load("../Output/SOEdef.jld", "sd")

mpe_iter!(sd, run_number = 1)

save("../Output/SOEdef.jld", "sd", sd)