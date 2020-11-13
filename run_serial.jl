include("main_serial.jl")

write("../Output/output.txt", "")
sd = load("../Output/SOEdef.jld", "sd")
update_probs!(sd)

mpe_iter!(sd, run_number = 1)

g, p_bench, _, v_m, def_freq = make_simulated_path(sd, "../Output/run1/", 30000);
Wr = mean([mean(series(p, :Wr)) for p in p_bench])

for (key, val) in pars(sd)
	print_save("$(rpad(key, 6, " ")) => $val\n")
end

save("../Output/SOEdef.jld", "sd", sd, "g", g, "pp", p_bench, "pars", pars(sd), "Wr", Wr)
print_save("$(sd.gr[:z])\n")

sd2 = load("../Output/SOEdef_nodef.jld", "sd")
for (key, val) in pars(sd)
	sd2.pars[key] = val
end
update_probs!(sd2)
mpe_iter!(sd2, run_number = 2, nodef = true)
g, p_nodef, _, v_m, def_freq = make_simulated_path(sd2, "../Output/run1/", 30000);
Wr = mean([mean(series(p, :Wr)) for p in p_nodef])

print_save("$(sd2.gr[:z])\n")

save("../Output/SOEdef_nodef.jld", "sd", sd2, "g", g, "pp", p_bench, "pars", pars(sd2), "Wr", Wr)
