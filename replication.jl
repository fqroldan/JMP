### Make sure folder "Replication" exists and is empty
print("\n\nReplication package for 'The Aggregate-Demand Doom Loop.'\nLast updated Feb 9, 2022\n\n\n")
print("Make sure folder '../Replication' exists and is empty and '../Output/' contains SOEdef.jld2, SOEdef_nodef.jld2, and IRF.jld2 (otherwise give 'resolve_resimulate()' with appropriate loaddir). Then give 'replicate()'.\n\n")

include("main_serial.jl")
include("data_jmp.jl")
include("data_rates.jl")
include("minimal_oneagent.jl")

function resolve_resimulate(folder = "../Replication/"; loaddir = "../Output/")
    sd_bench = load(loaddir*"SOEdef.jld2", "sd");
    mpe_iter!(sd_bench, run_number = 1)
    K = 500
    T = 200*K
    g, p_bench, _, _, _ = make_simulated_path(sd_bench, loaddir*"run1/", T, K = K);
    Wr_bench = mean([mean(series(p, :Wr)) for p in p_bench])
    save(folder * "SOEdef.jld2", "sd", sd_bench, "pp", p_bench, "Wr", Wr_bench)

    sd_nodef = load(loaddir*"SOEdef_nodef.jld2", "sd")
    mpe_iter!(sd_nodef, run_number = 2, nodef = true)
    g, p_nodef, _, _, _ = make_simulated_path(sd_nodef, loaddir*"run2/", T, K = K)
    Wr_nodef = mean([mean(series(p, :Wr)) for p in p_nodef])
    save(folder * "SOEdef_nodef.jld2", "sd", sd_nodef, "pp", p_nodef, "Wr", Wr_nodef)

    pIRF_bench, t1, t2, pIRF_nodef, pIRF_samep = IRF_default(sd_bench, sd_nodef, 1, 11, 9, B0 = 5.5, K = 5000)
    save(folder * "IRF.jld2", "pIRF_bench", pIRF_bench, "t1", t1, "t2", t2, "pIRF_nodef", pIRF_nodef, "pIRF_samep", pIRF_samep)

    nothing
end


function replicate(folder = "../Replication/"; loaddir = "../Output/")
    df_all = load_all()

    sd_bench, p_bench, W_bench = load(loaddir * "SOEdef.jld2", "sd", "pp", "Wr")
    sd_nodef, p_nodef, W_nodef = load(loaddir * "SOEdef_nodef.jld2", "sd", "pp", "Wr")

    pIRF_bench, t1, t2, pIRF_nodef, pIRF_samep = load(loaddir * "IRF.jld2", "pIRF_bench", "t1", "t2", "pIRF_nodef", "pIRF_samep")

    freq_bench = get_def_freq(p_bench)
    v_nodef = simul_stats(p_nodef)
    freq_nodef = get_def_freq(p_nodef)


    # Figure 1: Spanish Output and Consumption in the 2000s
    fig1 = SPA_CvY(style = paper, sh = false)
    savefig(fig1, folder * "YvsC_paper.pdf", width = 900, height = 300)

    # Figure 2: Net Worth of Spanish Households
    fig2 = make_nw(with_annot = false, style = paper)
    savefig(fig2, folder * "networth_ALN_SPA.pdf", width = 900, height = 300)

    # Table 1: Correlation of Spreads and Macroeconomic Outcomes
    regs_sec2(df_all, savedir = folder)

    # Figure 3: Anticipation of TFP costs of default
    sm = SOEmin()
    fig3 = makeplots_minimal(sm, move = "d", style = paper)
    savefig(fig3, folder * "minimal_tfp_d.pdf", width = 900, height = 300)

    # Figure 4: Anticipation of redistribution in case of default
    fig4 = minimal_twoagents(style = paper)
    savefig(fig4, folder * "minimal_2agents.pdf", width = 900, height = 300)

    # Figure 5: Labor Income and Expected returns
    fig5 = make_panels(sd_bench, "Earnings_Default", style = paper, leg = false)
    savefig(fig5, folder * "wages_paper.pdf", width = 800, height = 400)

    # Table 2: Estimated Fiscal Rules
    fig19 = regs_fiscalrules(df_all, savedir = folder, style = paper)
    savefig(fig19, folder * "fiscalrules_paper.pdf", width = 800, height = 400)

    # Table 3: Parameter Values
    write(folder * "params_table.txt", make_params_table(sd_bench))

    # Table 4: Model Fit
    v_bench = simul_stats(p_bench)
    write(folder * "calib_table.txt", make_calib_table(v_bench))

    # Figure 6: Welfare Functions
    fig6 = make_panels(sd_bench, "Wr-Wd", style = paper, leg = false)
    savefig(fig6, folder * "WrWd_paper.pdf", width = 800, height = 400)

    # Figure 7: Price of Debt
    fig7 = make_debtprice(sd_bench, style = paper, leg = false)
    savefig(fig7, folder * "debtprice_paper.pdf", width = 900, height = 350)

    # Figure 8: Unemployment
    fig8 = make_unemp(sd_bench, style = paper, leg = false)
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
    fig9 = panels_crises_small(p_bench, 400, :spread, k = 1, k_back = 11, thres_back = 350, style = paper)
    savefig(fig9, folder * "panels_crises_paper.pdf", width = 1100, height = 550)

    # Figure 10: Crises
    fig10 = panels_comp(p_bench, p_nodef, 400, :spread, thres_back = 350, k = 1, k_back = 11, style = paper, yh = 0.55, relative = true)
    savefig(fig10, folder * "panels_comp_paper.pdf", width = 1100, height = 550)

    # Figure 11: Crises across the Distribution
    fig11 = dist_CW(p_bench, 400, :spread, k = 1, k_back = 11, thres_back = 350, style = paper)
    savefig(fig11, folder * "dist_CW_paper.pdf", width = 900, height = 350)

    # Figure 12: Distributional Amplification
    fig12 = distribution_crises_new(p_bench, 400, :spread, k = 1, k_back = 11, thres_back = 350, style = paper)
    savefig(fig12, folder * "distribution_crises_paper.pdf", width = 600, height = 400)

    # Figure 13: Default-risk IRF
    fig13 = panels_IRF(pIRF_bench, pIRF_nodef, pIRF_samep, height = 900 * 0.65, width = 1900 * 0.65, cond_Y = 0.9474, style = paper)
    savefig(fig13, folder * "defaultriskIRF_paper.pdf", width = 1100, height = 550)

    # Figure 14: Value functions in the crisis
    fig14 = distribution_IRF(pIRF_bench, pIRF_nodef, pIRF_samep, height = 900 * 0.65, width = 1900 * 0.65, cond_Y = 0.9474, style = paper)
    savefig(fig14, folder * "distribIRF_paper.pdf", width = 800, height = 400)

    # Figure 15: Subjective probabilities of default
    fig15 = make_twisted(sd_bench, style = paper, leg = true)
    savefig(fig15, folder * "twistedprobs_paper.pdf", width = 800, height = 350)

    # Figure 16: Market Clearing in the Nontradable Sector
    ## This one is done in TikZ on the paper itself

    # Table 7: Discrepancies in Simulation

    # Figure 17: Transfers
    fig17 = make_panels(sd_bench, "T", style = paper, leg = false)
    savefig(fig17, folder * "transfers_paper.pdf", width = 800, height = 400)

    # Figure 18: Factors limiting production
    fig18 = make_FLP(; style = paper)
    savefig(fig18, folder * "factorsLP_paper.pdf", width = 900, height = 350)

    # Figure 19: Estimated Fiscal rules
    ## Done along with table 2

    # Figure 20: Interest rates in Spain
    fig20_bor = plot_borrowing(style = paper)
    savefig(fig20_bor, folder * "borrowingrates_paper.pdf", width = 1000, height = 500)
    fig20_dep = plot_deposit(style = paper)
    savefig(fig20_dep, folder * "depositrates_paper.pdf", width = 1000, height = 500)
end