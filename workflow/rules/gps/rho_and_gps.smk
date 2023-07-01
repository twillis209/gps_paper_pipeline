rule simulate_correlated_p_values_for_rho:
    output:
        temp("results/rho_simulations/{rho}/{zmean}_{zsd}/{rep}/pvalues.tsv")
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
        temp("results/rho_simulations/{rho}/{zmean}_{zsd}/{rep}/gps_value.tsv")
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
        temp("results/rho_simulations/{rho}/{zmean}_{zsd}/{rep}/{draws}_draws/permutations.tsv")
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
        temp("results/rho_simulations/{rho}/{zmean}_{zsd}/{rep}/li_gps_pvalue.tsv")
    localrule: True
    script:
        "../../scripts/compute_li_gps_pvalue.R"

rule fit_gev_and_compute_gps_pvalue_for_simulated_correlated_p_values:
    input:
        gps_file = "results/rho_simulations/{rho}/{zmean}_{zsd}/{rep}/gps_value.tsv",
        perm_file = "results/rho_simulations/{rho}/{zmean}_{zsd}/{rep}/{draws}_draws/permutations.tsv"
    output:
        temp("results/rho_simulations/{rho}/{zmean}_{zsd}/{rep}/{draws}_draws/gps_pvalue.tsv"),
    params:
        trait_A = 'V1',
        trait_B = 'V2'
    localrule: True
    conda: "../../envs/gps_paper.yaml"
    script:
        "../../scripts/fit_gev_and_compute_gps_pvalue.R"

rule collate_gps_results_for_varying_rho_simulations:
    input:
        gev = [f"results/rho_simulations/{rho}/{zmean}_{zsd}/{rep}/3000_draws/gps_pvalue.tsv" for rho in ["0_0", "0_05", "0_1", "0_15", "0_2", "0_25"] for zmean in [1, 2, 3] for zsd in [1,2] for rep in range(1, 101)],
        li = [f"results/rho_simulations/{rho}/{zmean}_{zsd}/{rep}/li_gps_pvalue.tsv" for rho in ["0_0", "0_05", "0_1", "0_15", "0_2", "0_25"] for zmean in [1, 2, 3] for zsd in [1,2] for rep in range(1, 101)]

        # TODO amend to bring in line with rest of pipeline
rule plot_rho_effect_on_gps_type_1_error:
    output:
        "results/rho_simulations/correlated_null_simulations.png"
    threads: 1
    localrule: True
    conda: "../../envs/gps_paper.yaml"
    script:
        "../../scripts/gps/plot_rho_effect_on_gps.R"
