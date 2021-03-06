include("main_serial.jl")

write("../Output/output.txt", "")
sd = FileIO.load("../Output/SOEdef.jld2", "sd")
update_probs!(sd)

mpe_iter!(sd, run_number = 1)

g, p_bench, _, v_m, def_freq = make_simulated_path(sd, "../Output/run1/", 30000);
Wr = mean([mean(series(p, :Wr)) for p in p_bench])

for (key, val) in pars(sd)
	print_save("$(rpad(key, 6, " ")) => $val\n")
end

FileIO.save("../Output/SOEdef.jld2", "sd", sd, "g", g, "pp", p_bench, "pars", pars(sd), "Wr", Wr)
print_save("$(sd.gr[:z])\n")

sd2 = FileIO.load("../Output/SOEdef_nodef.jld2", "sd")
for (key, val) in pars(sd)
	sd2.pars[key] = val
end
update_probs!(sd2)
mpe_iter!(sd2, run_number = 2, nodef = true)
g, p_nodef, _, v_m, def_freq = make_simulated_path(sd2, "../Output/run2/", 30000);
Wr = mean([mean(series(p, :Wr)) for p in p_nodef])

print_save("$(sd2.gr[:z])\n")
FileIO.save("../Output/SOEdef_nodef.jld2", "sd", sd2, "g", g, "pp", p_nodef, "pars", pars(sd2), "Wr", Wr)



# g, p_benchnodef, _, _, _ = make_simulated_path(sd, "../Output/run5/", 30000, ϕ = sd2.ϕ);
# Wr = mean([mean(series(p, :Wr)) for p in p_benchnodef])
# FileIO.save("../Output/SOEdef_benchnodef.jld2", "g", g, "pp", p_benchnodef, "Wr", Wr)


# sd3 = FileIO.load("../Output/SOEdef.jld2", "sd")
# for (key, val) in pars(sd)
# 	sd3.pars[key] = val
# end
# update_probs!(sd3)
# mpe_iter!(sd3, run_number = 3, nodef = false, noΔ = true)
# g, p_noΔ, _, v_m, def_freq = make_simulated_path(sd3, "../Output/run3/", 30000);
# Wr = mean([mean(series(p, :Wr)) for p in p_noΔ])
# FileIO.save("../Output/SOEdef_nodelta.jld2", "sd", sd3, "g", g, "pp", p_noΔ, "pars", pars(sd3), "Wr", Wr)

# sd4 = FileIO.load("../Output/SOEdef.jld2", "sd")
# for (key, val) in pars(sd)
# 	sd4.pars[key] = val
# end
# update_probs!(sd4)
# mpe_iter!(sd4, run_number = 4, nodef = false, noΔ = false, nob = true)
# g, p_nob, _, v_m, def_freq = make_simulated_path(sd4, "../Output/run4/", 30000);
# Wr = mean([mean(series(p, :Wr)) for p in p_nob])
# FileIO.save("../Output/SOEdef_nob.jld2", "sd", sd4, "g", g, "pp", p_nob, "pars", pars(sd4), "Wr", Wr)


