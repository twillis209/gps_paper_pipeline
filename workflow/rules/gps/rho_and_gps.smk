import pandas as pd

rule simulate_correlated_p_values_for_rho:
    output:
        "results/rho_simulations/{rho}/{zmean}_{zsd}/{rep}/pvalues.tsv"
    params:
        rho = lambda w: float(w.rho.replace('_', '.')),
        zmean = lambda w: float(w.zmean),
        zsd = lambda w: float(w.zsd),
        N00 = 4e4,
        N01 = 400,
        N10 = 400
    localrule: True
    conda: "../../envs/gps_paper.yaml"
    script:
        "../../scripts/gps/simulate_correlated_p_values.R"

rule compute_gps_for_simulated_correlated_p_values:
    input:
        "results/rho_simulations/{rho}/{zmean}_{zsd}/{rep}/pvalues.tsv"
    output:
        "results/rho_simulations/{rho}/{zmean}_{zsd}/{rep}/gps_value.tsv"
    params:
        a_colname = 'V1',
        b_colname = 'V2'
    localrule: True
    shell:
        "workflow/scripts/gps_cpp/build/apps/computeGpsCLI -i {input} -a {params.a_colname} -b {params.b_colname} -c {params.a_colname} -d {params.b_colname} -n {threads} -f pp -o {output}"

rule permute_for_simulated_correlated_p_values:
    input:
        "results/rho_simulations/{rho}/{zmean}_{zsd}/{rep}/pvalues.tsv"
    output:
        "results/rho_simulations/{rho}/{zmean}_{zsd}/{rep}/{draws}_draws/permutations.tsv"
    params:
        a_colname = 'V1',
        b_colname = 'V2'
    threads: 2
    resources:
        mem_mb = get_mem_mb,
        runtime = 20
    group: "rho_permutation"
    shell:
        "workflow/scripts/gps_cpp/build/apps/permuteTraitsCLI -i {input} -o {output} -a {params.a_colname} -b {params.b_colname} -c {threads} -n {wildcards.draws}"

rule compute_li_gps_pvalue_for_simulated_correlated_p_values:
    input:
        "results/rho_simulations/{rho}/{zmean}_{zsd}/{rep}/gps_value.tsv"
    output:
        "results/rho_simulations/{rho}/{zmean}_{zsd}/{rep}/li_gps_pvalue.tsv"
    localrule: True
    script:
        "../../scripts/compute_li_gps_pvalue.R"

rule fit_gev_and_compute_gps_pvalue_for_simulated_correlated_p_values:
    input:
        gps_file = "results/rho_simulations/{rho}/{zmean}_{zsd}/{rep}/gps_value.tsv",
        perm_file = "results/rho_simulations/{rho}/{zmean}_{zsd}/{rep}/{draws}_draws/permutations.tsv"
    output:
        "results/rho_simulations/{rho}/{zmean}_{zsd}/{rep}/{draws}_draws/gps_pvalue.tsv"
    params:
        trait_A = 'V1',
        trait_B = 'V2'
    localrule: True
    conda: "../../envs/gps_paper.yaml"
    script:
        "../../scripts/fit_gev_and_compute_gps_pvalue.R"

rule collate_gps_test_statistics_for_varying_rho_simulations:
    input:
        gev = [f"results/rho_simulations/{rho}/{zmean}_{zsd}/{rep}/3000_draws/gps_pvalue.tsv" for rho in ["0_0", "0_05", "0_1", "0_15", "0_2", "0_25"] for zmean in [1, 2, 3] for zsd in [1,2] for rep in range(1, 101)],
        li = [f"results/rho_simulations/{rho}/{zmean}_{zsd}/{rep}/li_gps_pvalue.tsv" for rho in ["0_0", "0_05", "0_1", "0_15", "0_2", "0_25"] for zmean in [1, 2, 3] for zsd in [1,2] for rep in range(1, 101)]
    output:
        "results/rho_simulations/compiled_test_statistics.tsv"
    localrule: True
    run:
        gev = []

        for x in input.gev:
            print(x.split('/'))
            rho, zmean_zsd, rep = x.split('/')[2:5]

            rho = float(rho.replace('_', '.'))
            zmean = int(zmean_zsd.split('_')[0])
            zsd = int(zmean_zsd.split('_')[1])
            rep = int(rep)

            with open(x, 'r') as infile:
                line = infile.readline()
                line = infile.readline()
                gps = line.split()[0]
                pvalue = line.split()[8]

            gev.append(
                {
                    'rho' : rho,
                    'zmean' : zmean,
                    'zsd' : zsd,
                    'rep' : rep,
                    'gps' : gps,
                    'pvalue' : pvalue,
                    'stat' : 'GPS-GEV'
                }
             )

        gev_daf = pd.DataFrame(gev)

        li = []

        for x in input.li:
            rho, zmean_zsd, rep = x.split('/')[2:5]

            rho = float(rho.replace('_', '.'))
            zmean = int(zmean_zsd.split('_')[0])
            zsd = int(zmean_zsd.split('_')[1])
            rep = int(rep)

            with open(x, 'r') as infile:
                line = infile.readline()
                line = infile.readline()
                gps = line.split()[0]
                pvalue = line.split()[1]

                li.append(
                    {
                        'rho' : rho,
                        'zmean' : zmean,
                        'zsd' : zsd,
                        'rep' : rep,
                        'gps' : gps,
                        'pvalue' : pvalue,
                        'stat' : 'GPS-Exp'
                    }
                )

        li_daf = pd.DataFrame(li)

        pd.concat([gev_daf, li_daf]).to_csv(output[0], index = False, sep = '\t')

rule plot_fig_s10:
    input:
        test_statistics = "results/rho_simulations/compiled_test_statistics.tsv",
        pvalue_files = ["results/rho_simulations/0_0/1_1/1/pvalues.tsv",
                      "results/rho_simulations/0_0/2_1/1/pvalues.tsv",
                      "results/rho_simulations/0_0/3_1/1/pvalues.tsv",
                      "results/rho_simulations/0_0/1_2/1/pvalues.tsv",
                      "results/rho_simulations/0_0/2_2/1/pvalues.tsv",
                      "results/rho_simulations/0_0/3_2/1/pvalues.tsv"
                      ]
    output:
        "results/rho_simulations/fig_s10.png"
    threads: 1
    localrule: True
    conda: "../../envs/gps_paper.yaml"
    script:
        "../../scripts/gps/plot_fig_s10.R"
