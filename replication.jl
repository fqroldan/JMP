### Make sure folder "Replication" exists and is empty
print("The Aggregate-Demand Doom Loop. Replication package. Last updated Feb 8, 2022\n")
print("Make sure folder 'Replication' exists and is empty. Then give 'replicate()'.\n")

include("main_serial.jl")
include("data_jmp.jl")
include("data_rates.jl")

function resolve_all(folder = "../Replication/")
    sd_bench = load("../Output/SOEdef.jld", "sd")
    mpe_iter!(sd_bench, run_number = 1)
    g, p_bench, _, _, _ = make_simulated_path(sd_bench, "../Output/run1/", 100_000, K = 500)
    Wr_bench = mean([mean(series(p, :Wr)) for p in p_bench])
    save(folder * "SOEdef.jld2", "sd", sd_bench, "pp", p_bench, "Wr", Wr_bench)

    sd_nodef = load("../Output/SOEdef.jld", "sd")
    mpe_iter!(sd_nodef, run_number = 2, nodef = true)
    g, p_nodef, _, _, _ = make_simulated_path(sd_nodef, "../Output/run2/", 10_000, K = 50)
    Wr_nodef = mean([mean(series(p, :Wr)) for p in p_nodef])
    save(folder * "SOEdef_nodef.jld2", "sd", sd_nodef, "pp", p_nodef, "Wr", Wr_nodef)

    pIRF_bench, t1, t2, pIRF_nodef, pIRF_samep = IRF_default(sd_bench, sd_nodef, 1, 11, 9, B0 = 5.5, K = 5000)
    save(folder * "IRF.jld2", "pIRF_bench", pIRF_bench, "t1", t1, "t2", t2, "pIRF_nodef", pIRF_nodef, "pIRF_samep", pIRF_samep)
end


function replicate(folder = "../Replication/")
    df_all = load_all()

    sd_bench, p_bench, W_bench = load("../Output/SOEdef.jld2", "sd", "pp", "Wr")
    sd_nodef, p_nodef, W_nodef = load("../Output/SOEdef_nodef.jld2", "sd", "pp", "Wr")

    pIRF_bench, t1, t2, pIRF_nodef, pIRF_samep = load("../Output/IRF.jld2", "pIRF_bench", "t1", "t2", "pIRF_nodef", "pIRF_samep")

    freq_bench = get_def_freq(p_bench)
    v_nodef = simul_stats(p_nodef)
    freq_nodef = get_def_freq(p_nodef)


    # Figure 1: Spanish Output and Consumption in the 2000s
    fig1 = SPA_CvY(style = paper, sh = false)

    # Figure 2: Net Worth of Spanish Households
    fig2 = make_nw(with_annot = false, style = paper)

    # Table 1: Correlation of Spreads and Macroeconomic Outcomes
    regs_sec2(df_all, savedir = folder)

    # Figure 3: Anticipation of TFP costs of default

    # Figure 4: Anticipation of redistribution in case of default




    # Figure 5: Labor Income and Expected returns
    fig5 = make_panels(sd_bench, "Earnings_Default", style = paper, leg = false)

    # Table 2: Estimated Fiscal Rules
    fig18 = regs_fiscalrules(df_all, savedir = folder)

    # Table 3: Parameter Values

    # Table 4: Model Fit
    v_bench = simul_stats(p_bench)
    write(folder * "calib_table.txt", make_calib_table(v_bench))

    # Figure 6: Welfare Functions
    fig6 = make_panels(sd_bench, "Wr-Wd", style = paper, leg = false)

    # Figure 7: Price of Debt
    fig7 = make_debtprice(sd_bench, style = paper, leg = false)

    # Figure 8: Unemployment
    fig8 = make_unemp(sd_bench, style = paper, leg = false)

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

    # Figure 10: Crises
    fig10 = panels_comp(p_bench, p_nodef, 400, :spread, thres_back = 350, k = 1, k_back = 11, style = paper, yh = 0.55, relative = true)

    # Figure 11: Distributional Amplification
    fig11 = distribution_crises_new(p_bench, 400, :spread, k = 1, k_back = 11, thres_back = 350, style = paper)

    # Figure 12: Crises across the Distribution
    fig12 = dist_CW(p_bench, 400, :spread, k = 1, k_back = 11, thres_back = 350, style = paper)

    # Figure 13: Default-risk IRF
    fig13 = panels_IRF(pIRF_bench, pIRF_nodef, pIRF_samep, height = 900 * 0.65, width = 1900 * 0.65, cond_Y = 0.96, style = paper)

    # Figure 14: Value functions in the crisis
    fig14 = distribution_IRF(pIRF_bench, pIRF_nodef, pIRF_samep, height = 900 * 0.65, width = 1900 * 0.65, cond_Y = 0.96, style = paper)

    # Figure 15: Subjective probabilities of default
    fig15 = make_twisted(sd_bench, style = paper, leg = true)

    # Figure 16: Market Clearing in the Nontradable Sector
    ## This one is done in TikZ on the paper itself

    # Table 7: Discrepancies in Simulation

    # Figure 17: Transfers
    fig17 = make_panels(sd_bench, "T", style = paper, leg = false)

    # Figure 18: Slack in the Spanish economy

    # Figure 19: Estimated Fiscal rules

    # Figure 20: Interest rates in Spain
    fig20_bor = plot_borrowing(style = paper)
    fig20_dep = plot_deposit(style = paper)
end