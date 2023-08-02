### Make sure folder "Replication" exists and is empty
print("\n\nReplication package for 'The Aggregate-Demand Doom Loop.'\nLast updated Jan 29, 2023\n\n\n")
print("Make sure folder '../Replication' exists and is empty and '../Output/' contains SOEdef.jld2, SOEdef_nodef.jld2, and IRF.jld2 (otherwise give 'resolve_resimulate()' with appropriate loaddir). Then give 'replicate()'.\n\n")

include("main_serial.jl")
include("data_jmp.jl")
include("data_rates.jl")
# include("minimal_oneagent.jl")

function resolve_resimulate(folder = "../Replication/"; loaddir = "../Output/", datadir = "../Data/")
    sd_bench = load(loaddir*"SOEdef.jld2", "sd");
    mpe_iter!(sd_bench, run_number = 1)
    K = 500
    T = 400*K
    g, p_bench, _, _, _ = make_simulated_path(sd_bench, loaddir*"run1/", T, K = K, datadir = datadir);
    Wr_bench = mean([mean(series(p, :Wr)) for p in p_bench])
    save(folder * "SOEdef.jld2", "sd", sd_bench, "pp", p_bench, "Wr", Wr_bench)

    sd_nodef = load(loaddir*"SOEdef_nodef.jld2", "sd");
    mpe_iter!(sd_nodef, run_number = 2, nodef = true)
    g, p_nodef, _, _, _ = make_simulated_path(sd_nodef, loaddir*"run2/", T, K = K, datadir = datadir);
    Wr_nodef = mean([mean(series(p, :Wr)) for p in p_nodef])
    save(folder * "SOEdef_nodef.jld2", "sd", sd_nodef, "pp", p_nodef, "Wr", Wr_nodef)

    pIRF_bench, t1, t2, pIRF_nodef, pIRF_samep = IRF_default(sd_bench, sd_nodef, 1, 11, 9, B0 = 4, K = 5000) # used to be B0=5.5, K = 5000
    # save(folder * "IRF.jld2", "pIRF_bench", pIRF_bench, "t1", t1, "t2", t2, "pIRF_nodef", pIRF_nodef, "pIRF_samep", pIRF_samep)
    save("../Rep2/IRF.jld2", "pIRF_bench", pIRF_bench, "t1", t1, "t2", t2, "pIRF_nodef", pIRF_nodef, "pIRF_samep", pIRF_samep)


    sd_hi = load(loaddir*"SOEdef.jld2", "sd");
    sd_hi.pars[:τ] = 0.36
    comp_eqm!(sd_hi, re_q = false, verbose = true)
    sd_lo = load(loaddir*"SOEdef.jld2", "sd");
    sd_lo.pars[:τ] = 0.26
    comp_eqm!(sd_lo, re_q = false, verbose = true)
    save("../Rep2/SOEdef_alt.jld2", "sd_hi", sd_hi, "sd_lo", sd_lo)
    # g_alt, p_alt, _, _, _ = make_simulated_path(sd_alt, loaddir*"run2/", T, K = K, datadir = datadir);
    # Wr_alt = mean([mean(series(p, :Wr)) for p in p_alt])
    # save("../Rep2/SOEdef_alt.jld2", "sd", sd_alt, "pp", p_alt, "Wr", Wr_alt)

    sd_hi, sd_lo = load("../Rep2/SOEdef_alt.jld2", "sd_hi", "sd_lo");

    pIRF_bench1, _, _, pIRF_hi, pIRF_lo = IRF_default_comp(sd_bench, sd_nodef, sd_hi, sd_lo, 1, 11, 9, B0=4, K = 5000) # Default B = 4
    save("../Rep2/IRF_cs.jld2", "pIRF_bench", pIRF_bench1, "pIRF_hi", pIRF_hi, "pIRF_lo", pIRF_lo)


    # panels_IRF(pIRF_bench1, pIRF_hi, pIRF_lo, cond_Y=0.96, slides=false, name_samep="<i>τ</i> = $(sd_lo.pars[:τ])")

    # pIRF_bench2, pIRF_nodef2, pIRF_alt = load("../Rep2/IRF_nodefcost_noq.jld2", "pIRF_bench", "pIRF_nodef", "pIRF_alt");
    nothing
end

# Need: output, gini, spreads, welfare/pN for bench, hi, and lo
# Need: unemployment out, Gini in fig 12


function replicate(folder = "../Replication/"; loaddir = "../Output/", datadir = "../Data/")
    df_all = load_all(datadir);

    sd_bench, p_bench, W_bench = load(loaddir * "SOEdef.jld2", "sd", "pp", "Wr");
    sd_nodef, p_nodef, W_nodef = load(loaddir * "SOEdef_nodef.jld2", "sd", "pp", "Wr");

    pIRF_bench, t1, t2, pIRF_nodef, pIRF_samep = load(loaddir * "/IRF.jld2", "pIRF_bench", "t1", "t2", "pIRF_nodef", "pIRF_samep");

    pIRF_bench1, pIRF_hi, pIRF_lo = load(loaddir * "/IRF_cs.jld2", "pIRF_bench", "pIRF_hi", "pIRF_lo");

    freq_bench = get_def_freq(p_bench)
    v_nodef = simul_stats(p_nodef)
    freq_nodef = get_def_freq(p_nodef)


    # Figure 1: Spanish Output and Consumption in the 2000s
    fig1 = SPA_CvY(loaddir = datadir, slides = false, sh = false)
    savefig(fig1, folder * "YvsC_paper.pdf", width = 900, height = 300)

    # Figure 2: Net Worth of Spanish Households
    fig2 = make_nw(loaddir = datadir, with_annot = false, slides = false)
    savefig(fig2, folder * "networth_ALN_SPA.pdf", width = 900, height = 300)

    fig23_app = make_nw_levels(;levels = true, slides = false)
    savefig(fig23_app, folder * "networth_levels.pdf", width = 900, height = 350)

    # Table 1: Correlation of Spreads and Macroeconomic Outcomes
    regs_sec2(df_all, savedir = folder)

    # Figure 3: Anticipation of TFP costs of default
    sm = SOEmin()
    fig3 = makeplots_minimal(sm, move = "d", slides = false)
    savefig(fig3, folder * "minimal_tfp_d.pdf", width = 900, height = 300)

    sm_norep = SOEmin(ξ_d = 0)
    fig15_app = makeplots_minimal(sm_norep, move="d", slides = false)
    savefig(fig15_app, folder * "minimal_tfp_d_norep.pdf", width=900, height=300)

    # Figure 4: Anticipation of redistribution in case of default
    fig4 = minimal_twoagents(slides = false)
    savefig(fig4, folder * "minimal_2agents.pdf", width = 900, height = 300)

    # Figure 5: Labor Income and Expected returns
    fig5 = make_panels(sd_bench, "Earnings_Default", slides = false, leg = false)
    savefig(fig5, folder * "wages_paper.pdf", width = 800, height = 300)

    # Table 2: Estimated Fiscal Rules
    fig19 = regs_fiscalrules(df_all, savedir = folder, slides = false)
    savefig(fig19, folder * "fiscalrules_paper.pdf", width = 800, height = 400)

    # Table 3: Parameter Values
    write(folder * "params_table.txt", make_params_table(sd_bench))

    # Table 4: Model Fit
    v_bench = simul_stats(p_bench)
    write(folder * "calib_table.txt", make_calib_table(v_bench, loaddir = datadir))

    # Figure 6: Welfare Functions
    fig6 = make_panels(sd_bench, "Wr-Wd", slides = false, leg = false)
    savefig(fig6, folder * "WrWd_paper.pdf", width = 800, height = 300)

    # Figure 7: Price of Debt
    fig7 = make_debtprice(sd_bench, slides = false, leg = false)
    savefig(fig7, folder * "debtprice_paper.pdf", width = 900, height = 350)

    # Figure 8: Unemployment
    fig8 = make_unemp(sd_bench, slides = false, leg = false)
    savefig(fig8, folder * "unemp_paper.pdf", width = 900, height = 350)

    # Table 5: Models
    write(folder * "calib_table_comp.txt",
        make_calib_table_comp(
            [v_bench; 100 * freq_bench; W_bench],
            [v_nodef; 100 * freq_nodef; W_nodef],
        )
    )

    # Table 6: The Welfare Costs of Sovereign Risk
    write(folder * "welfare_table.txt", make_welfare_table(p_bench, p_nodef))

    # Figure 9: Times of High Spreads
    fig9 = panels_crises_small(p_bench, 400, :spread, k = 1, k_back = 11, thres_back = 350, slides = false)
    savefig(fig9, folder * "panels_crises_paper.pdf", width = 1100, height = 550)

    # Figure 10: Crisis dynamics in model and data
    fig10 = panels_crises_data(p_bench, 400, :spread, k = 1, k_back = 11, thres_back = 350, slides = false)
    savefig(fig10, folder * "panels_wdata_paper.pdf", width=900, height=500)

    # Figure 11: Crises
    fig11 = panels_comp(p_bench, p_nodef, 400, :spread, thres_back = 350, k = 1, k_back = 11, slides = false, yh = 0.55, relative = true)
    savefig(fig11, folder * "panels_comp_paper.pdf", width = 1100, height = 550)

    # Figure 12: Crises across the Distribution
    fig12 = dist_CW(p_bench, 400, :spread, k=1, k_back=11, thres_back=350, slides = false, cw = true)
    savefig(fig12, folder * "dist_CW_paper.pdf", width = 900, height = 350)

    fig12b = dist_CW(p_bench, 400, :spread, k=1, k_back=11, thres_back=350, slides = false, cw = false)
    savefig(fig12b, folder * "dist_AB_paper.pdf", width = 900, height = 350)

    # Figure 13: Default-risk IRF
    fig13 = panels_IRF(pIRF_bench, pIRF_nodef, pIRF_samep, cond_Y = 0.95, slides = false)
    savefig(fig13, folder * "defaultriskIRF_paper.pdf", width = 1100, height = 550)

    fig19_app = panels_IRF_wdata(pIRF_bench, pIRF_nodef, pIRF_samep, cond_Y = 0.95, slides=false)
    savefig(fig19_app, folder * "panelsIRF_wdata_paper.pdf", width=900, height=500)

    # Figure 14: Default-risk IRFs and tax progressivity
    fig14 = panels_IRF_cs(pIRF_bench1, pIRF_hi, pIRF_lo, cond_Y = 0.95, slides=false, τ_hi=sd_hi.pars[:τ], τ_lo=sd_lo.pars[:τ], give_stats=true)
    savefig(fig14, folder * "/panels_IRF_cs.pdf", width = 900, height = 500)

    # Figure 16: Value functions in the crisis
    fig16_app = distribution_IRF(pIRF_bench, pIRF_nodef, pIRF_samep, cond_Y = 0.95, slides = false)
    savefig(fig16_app, folder * "distribIRF_paper.pdf", width = 800, height = 400)

    # Figure 17: Subjective probabilities of default
    fig17 = make_twisted(sd_bench, slides = false, leg = true)
    savefig(fig17, folder * "twistedprobs_paper.pdf", width = 800, height = 350)

    # Figures 18 and 19 are above

    # Figure 20: Market Clearing in the Nontradable Sector
    ## This one is done in TikZ on the paper itself

    # Table 7: Discrepancies in Simulation
    # Get values from make_simulated_path(sd_bench,...)

    # Figure 21: Transfers
    fig21 = make_panels(sd_bench, "T", slides = false, leg = false)
    savefig(fig21, folder * "transfers_paper.pdf", width = 800, height = 250)

    # Figure 22: Factors limiting production
    fig22 = make_FLP(; slides = false)
    savefig(fig22, folder * "factorsLP_paper.pdf", width = 900, height = 350)

    # Figure 24: Estimated Fiscal rules
    ## Done along with table 2

    # Figure 25: Interest rates in Spain
    fig25_bor = plot_borrowing(slides = false)
    savefig(fig25_bor, folder * "borrowingrates_paper.pdf", width = 1000, height = 500)
    fig25_dep = plot_deposit(slides = false)
    savefig(fig25_dep, folder * "depositrates_paper.pdf", width = 1000, height = 500)
end
