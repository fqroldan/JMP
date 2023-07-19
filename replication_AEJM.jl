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

    pIRF_bench, t1, t2, pIRF_nodef, pIRF_samep = IRF_default(sd_bench, sd_nodef, 1, 11, 9, B0 = 4, K = 2500) # used to be B0=5.5, K = 5000
    save(folder * "IRF.jld2", "pIRF_bench", pIRF_bench, "t1", t1, "t2", t2, "pIRF_nodef", pIRF_nodef, "pIRF_samep", pIRF_samep)


    sd_alt = load(loaddir*"SOEdef.jld2", "sd");
    # sd_alt.pars[:τ] = 0.23
    sd_alt.pars[:τ] = 0.35
    comp_eqm!(sd_alt, re_q = false, verbose = true)
    # g_alt, p_alt, _, _, _ = make_simulated_path(sd_alt, loaddir*"run2/", T, K = K, datadir = datadir);
    # Wr_alt = mean([mean(series(p, :Wr)) for p in p_alt])
    # save("../Rep2/SOEdef_alt.jld2", "sd", sd_alt, "pp", p_alt, "Wr", Wr_alt)

    pIRF_bench2, t1, t2, pIRF_nodef2, pIRF_samep2, pIRF_alt = IRF_default_comp(sd_bench, sd_nodef, sd_alt, 1, 11, 9, B0 = 4, K = 1000);
    panels_IRF(pIRF_bench2, pIRF_nodef2, pIRF_alt, cond_Y=0.96, slides=false, name_samep="<i>τ</i> = $(sd_alt.pars[:τ])")
    panels_IRF(pIRF_bench2, pIRF_nodef2, pIRF_samep2, cond_Y = 0.96, slides = false)

    save("../Rep2/IRF_nodefcost_noq.jld2", "pIRF_bench", pIRF_bench2, "pIRF_nodef", pIRF_nodef2, "pIRF_alt", pIRF_alt)

    # pIRF_bench2, pIRF_nodef2, pIRF_alt = load("../Rep2/IRF_nodefcost_noq.jld2", "pIRF_bench", "pIRF_nodef", "pIRF_alt");
    nothing
end


function replicate(folder = "../Replication/"; loaddir = "../Output/", datadir = "../Data/")
    df_all = load_all(datadir);

    sd_bench, p_bench, W_bench = load(loaddir * "SOEdef.jld2", "sd", "pp", "Wr");
    sd_nodef, p_nodef, W_nodef = load(loaddir * "SOEdef_nodef.jld2", "sd", "pp", "Wr");

    pIRF_bench, t1, t2, pIRF_nodef, pIRF_samep = load(loaddir * "IRF.jld2", "pIRF_bench", "t1", "t2", "pIRF_nodef", "pIRF_samep");

    freq_bench = get_def_freq(p_bench)
    v_nodef = simul_stats(p_nodef)
    freq_nodef = get_def_freq(p_nodef)


    # Figure 1: Spanish Output and Consumption in the 2000s
    fig1 = SPA_CvY(loaddir = datadir, slides = false, sh = false)
    savefig(fig1, folder * "YvsC_paper.pdf", width = 900, height = 300)

    # Figure 2: Net Worth of Spanish Households
    fig2 = make_nw(loaddir = datadir, with_annot = false, slides = false)
    savefig(fig2, folder * "networth_ALN_SPA.pdf", width = 900, height = 300)

    # Table 1: Correlation of Spreads and Macroeconomic Outcomes
    regs_sec2(df_all, savedir = folder)

    # Figure 3: Anticipation of TFP costs of default
    sm = SOEmin()
    fig3 = makeplots_minimal(sm, move = "d", slides = false)
    savefig(fig3, folder * "minimal_tfp_d.pdf", width = 900, height = 300)

    sm_norep = SOEmin(ξ_d = 0)
    fig13_app = makeplots_minimal(sm_norep, move="d", slides = false)
    savefig(fig13_app, folder * "minimal_tfp_d_norep.pdf", width=900, height=300)

    # Figure 4: Anticipation of redistribution in case of default
    fig4 = minimal_twoagents(slides = false)
    savefig(fig4, folder * "minimal_2agents.pdf", width = 900, height = 300)

    # Figure 5: Labor Income and Expected returns
    fig5 = make_panels(sd_bench, "Earnings_Default", slides = false, leg = false)
    savefig(fig5, folder * "wages_paper.pdf", width = 800, height = 250)

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
    savefig(fig6, folder * "WrWd_paper.pdf", width = 800, height = 250)

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

    # Figure 10: Crises
    fig10 = panels_comp(p_bench, p_nodef, 400, :spread, thres_back = 350, k = 1, k_back = 11, slides = false, yh = 0.55, relative = true)
    savefig(fig10, folder * "panels_comp_paper.pdf", width = 1100, height = 550)

    # Figure 11: Crises across the Distribution
    fig11 = dist_CW(p_bench, 400, :spread, k=1, k_back=11, thres_back=350, slides = false, cw = true)
    savefig(fig11, folder * "dist_CW_paper.pdf", width = 900, height = 350)

    fig11b = dist_CW(p_bench, 400, :spread, k=1, k_back=11, thres_back=350, slides = false, cw = false)
    savefig(fig11b, folder * "dist_AB_paper.pdf", width = 900, height = 350)


    # Figure 12: Distributional Amplification
    # fig12 = distribution_crises_new(p_bench, 400, :spread, k = 1, k_back = 11, thres_back = 350, slides = false)
    # savefig(fig12, folder * "distribution_crises_paper.pdf", width = 600, height = 400)

    # Figure 12: Default-risk IRF
    fig12 = panels_IRF(pIRF_bench, pIRF_nodef, pIRF_samep, height = 900 * 0.65, width = 1900 * 0.65, cond_Y = 0.96, slides = false)
    savefig(fig12, folder * "defaultriskIRF_paper.pdf", width = 1100, height = 550)

    # fig12 = panels_IRF(pIRF_highτ, pIRF_highτ, pIRF_highτ, height = 900 * 0.65, width = 1900 * 0.65, cond_Y = 0.95, slides = false)

    # Figure 14: Value functions in the crisis
    fig14 = distribution_IRF(pIRF_bench, pIRF_nodef, pIRF_samep, height = 900 * 0.65, width = 1900 * 0.65, cond_Y = 0.96, slides = false)
    savefig(fig14, folder * "distribIRF_paper.pdf", width = 800, height = 400)

    # Figure 15: Subjective probabilities of default
    fig15 = make_twisted(sd_bench, slides = false, leg = true)
    savefig(fig15, folder * "twistedprobs_paper.pdf", width = 800, height = 350)

    # Figure 16: Market Clearing in the Nontradable Sector
    ## This one is done in TikZ on the paper itself

    # Table 7: Discrepancies in Simulation

    # Figure 17: Transfers
    fig17 = make_panels(sd_bench, "T", slides = false, leg = false)
    savefig(fig17, folder * "transfers_paper.pdf", width = 800, height = 250)

    # Figure 18: Factors limiting production
    fig18 = make_FLP(; slides = false)
    savefig(fig18, folder * "factorsLP_paper.pdf", width = 900, height = 350)

    # Figure 19: Estimated Fiscal rules
    ## Done along with table 2

    # Figure 20: Interest rates in Spain
    fig20_bor = plot_borrowing(slides = false)
    savefig(fig20_bor, folder * "borrowingrates_paper.pdf", width = 1000, height = 500)
    fig20_dep = plot_deposit(slides = false)
    savefig(fig20_dep, folder * "depositrates_paper.pdf", width = 1000, height = 500)
end



function compstats_inequality(; τr = 0.03, Nτ = 4, loaddir = "../Output/")
    sd_bench = load(loaddir * "SOEdef.jld2", "sd")
    sd_nodef = load(loaddir * "SOEdef_nodef.jld2", "sd")
    sd_alt = load(loaddir * "SOEdef.jld2", "sd")

    τ0 = sd_bench.pars[:τ]
    τvec1 = range(τ0, τ0 - τr, length=Nτ)
    τvec2 = range(τ0, τ0 + τr, length=Nτ)[2:end]

    Gmat = zeros(2*Nτ-1, 3)
    Ymat = zeros(2*Nτ-1, 3)
    Wmat = zeros(2*Nτ-1, 3)
    τvec = zeros(2*Nτ-1)
    
    for jτ in axes(Gmat,1)
        if jτ <= Nτ
            τv = τvec1[jτ]
        else
            if jτ == Nτ + 1
                sd_alt = load(loaddir * "SOEdef.jld2", "sd")
            end
            τv = τvec2[jτ - Nτ]
        end
        print("Iter $jτ\n")
        
        sd_alt.pars[:τ] = τv
        comp_eqm!(sd_alt, re_q = false, verbose = false)

        pIRF_bench, t1, t2, pIRF_nodef, _, pIRF_alt = IRF_default_comp(sd_bench, sd_nodef, sd_alt, 1, 11, 9, B0 = 4, K = 1000);
        Wr1, Wr2, G1, G2, Y1, Y2 = panels_IRF(pIRF_bench, pIRF_nodef, pIRF_alt, cond_Y = 0.96, give_stats = true)

        Wmat[jτ, :] .= Wr1, Wr2, Wr2 / Wr1 - 1
        Gmat[jτ, :] .= G1, G2, G2 / G1 - 1
        Ymat[jτ, :] .= Y1, Y2, Y2 / Y1 - 1
        τvec[jτ] = τv
    end

    return τvec, Gmat, Wmat, Ymat
end