include("main_serial.jl")

write("../Output/output.txt", "")
sd = load("../Output/SOEdef.jld", "sd")

mpe_iter!(sd, run_number = 1)

save("../Output/SOEdef.jld", "sd", sd)