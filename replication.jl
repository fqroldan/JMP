### Make sure folder "Replication" exists and is empty
print("The Aggregate-Demand Doom Loop. Replication package. Last updated Feb 8, 2022\n")
print("Make sure folder 'Replication' exists and is empty. Then give 'replicate()'.\n")

include("main_serial.jl")
include("data_jmp.jl")
include("data_rates.jl")

function replicate(folder = "../Replication/")
    df_all = load_all()

    sd_bench, p_bench, W_bench = load("../Output/SOEdef.jld2", "sd", "pp", "Wr")
    sd_nodef, p_nodef, W_nodef = load("../Output/SOEdef_nodef.jld2", "sd", "pp", "Wr")

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

    # Figure 11: Default-risk IRF
    pIRF_bench, t1, t2, pIRF_nodef, pIRF_samep = IRF_default(sd_bench, sd_nodef, 1, 11, 9, B0 = 5.5, K = 5000)

    fig11 = panels_IRF(pIRF_bench, pIRF_nodef, pIRF_samep, height = 900 * 0.65, width = 1900 * 0.65, cond_Y = 0.96, style = paper)

    # Figure 13: Value functions in the crisis
    fig13 = distribution_IRF(pIRF_bench, pIRF_nodef, pIRF_samep, height = 900 * 0.65, width = 1900 * 0.65, cond_Y = 0.96, style = paper)

    # Figure 14: Subjective probabilities of default
    fig14 = make_twisted(sd_bench, style = paper, leg = true)

    # Figure 15: Market Clearing in the Nontradable Sector
    ## This one is done in TikZ on the paper itself

    # Table 7: Discrepancies in Simulation

    # Figure 16: Transfers
    fig16 = make_panels(sd_bench, "T", style = paper, leg = false)

    # Figure 17: Slack in the Spanish economy

    # Figure 18: Estimated Fiscal rules

    # Figure 19: Interest rates in Spain
    fig18_bor = plot_borrowing(style = paper)
    fig18_dep = plot_deposit(style = paper)
end