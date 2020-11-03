include("main_serial.jl")

write("../Output/output.txt", "")
sd = load("../Output/SOEdef.jld", "sd")

mpe_iter!(sd, run_number = 1)

g, p_bench, _, v_m, def_freq = make_simulated_path(sd, "../Output/run1/", 30000);

for (key, val) in pars(sd)
	print_save("$(rpad(key, 6, " ")) => $val\n")
end

save("../Output/SOEdef.jld", "sd", sd)