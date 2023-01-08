import os

rule compute_gps_for_sim_pair:
    input:
        sum_stats_file = "results/simgwas/simulated_sum_stats/whole_genome_sum_stats/{no_reps}_reps/randomised/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/window_{window}_step_{step}/seed_{seed}_pruned_sum_stats_tags_{tag_A}-{tag_B}.tsv"
    output:
        temp("results/gps/simgwas/{no_reps}_reps/randomised/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/window_{window}_step_{step}/seed_{seed}_tags_{tag_A,\d+}-{tag_B,\d+}_gps_value.tsv")
    params:
        a_colname = lambda wildcards: f"p.{wildcards.tag_A}",
        b_colname = lambda wildcards: f"p.{wildcards.tag_B}",
        no_of_pert_iterations = 0
    threads: 12
    resources:
        runtime = 5
    priority: 1
    group: "ldsc_hoeffding_sumher_gps_sans_permutation"
    shell:
        "workflow/scripts/gps_cpp/build/apps/computeGpsCLI -i {input.sum_stats_file} -a {params.a_colname} -b {params.b_colname} -c {wildcards.effect_blocks_A} -d {wildcards.effect_blocks_B} -p {params.no_of_pert_iterations} -n {threads} -o {output}"

rule permute_sim_pair:
    input:
        sum_stats_file = "results/simgwas/simulated_sum_stats/whole_genome_sum_stats/{no_reps}_reps/randomised/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/window_{window}_step_{step}/seed_{seed}_pruned_sum_stats_tags_{tag_A}-{tag_B}.tsv"
    output:
        "results/gps/simgwas/{no_reps}_reps/randomised/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/window_{window}_step_{step}/{draws,\d+}_permutations/seed_{seed}_tags_{tag_A,\d+}-{tag_B,\d+}.tsv"
    params:
        a_colname = lambda wildcards: f"p.{wildcards.tag_A}",
        b_colname = lambda wildcards: f"p.{wildcards.tag_B}",
        no_of_pert_iterations = 1
    threads: 12
    resources:
        mem_mb = get_mem_mb,
        runtime = get_permute_time,
    group: "permutation"
    priority: 1
    shell:
        "workflow/scripts/gps_cpp/build/apps/permuteTraitsCLI -i {input.sum_stats_file} -o {output} -a {params.a_colname} -b {params.b_colname} -c {threads} -n {wildcards.draws} -p {params.no_of_pert_iterations}"

rule compute_li_gps_pvalue_for_sim_pair:
    input:
        "results/gps/simgwas/{no_reps}_reps/randomised/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/window_{window}_step_{step}/seed_{seed}_tags_{tag_A}-{tag_B}_gps_value.tsv",
    output:
        "results/gps/simgwas/{no_reps}_reps/randomised/{ncases_A,\d+}_{ncontrols_A,\d+}_{ncases_B,\d+}_{ncontrols_B,\d+}/{effect_blocks_A,[smlh\d-]+}_{effect_blocks_B,[smlh\d-]+}_{shared_effect_blocks,[smlh\d-]+}/window_{window}_step_{step}/{draws,\d+}_permutations/seed_{seed,\d+}_tags_{tag_A,\d+}-{tag_B,\d+}_li_gps_pvalue.tsv"
    group: "permutation"
    script:
        "../../scripts/compute_li_gps_pvalue.R"

rule fit_gev_and_compute_gps_pvalue_for_sim_pair:
    input:
        gps_file = "results/gps/simgwas/{no_reps}_reps/randomised/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/window_{window}_step_{step}/seed_{seed}_tags_{tag_A}-{tag_B}_gps_value.tsv",
        perm_file = "results/gps/simgwas/{no_reps}_reps/randomised/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/window_{window}_step_{step}/{draws}_permutations/seed_{seed}_tags_{tag_A}-{tag_B}.tsv"
    output:
        "results/gps/simgwas/{no_reps}_reps/randomised/{ncases_A,\d+}_{ncontrols_A,\d+}_{ncases_B,\d+}_{ncontrols_B,\d+}/{effect_blocks_A,[smlh\d-]+}_{effect_blocks_B,[smlh\d-]+}_{shared_effect_blocks,[smlh\d-]+}/window_{window}_step_{step}/{draws,\d+}_permutations/seed_{seed,\d+}_tags_{tag_A,\d+}-{tag_B,\d+}_gps_pvalue.tsv"
    resources:
        runtime = 2
    group: "permutation"
    shell:
      "Rscript workflow/scripts/fit_gev_and_compute_gps_pvalue.R -g {input.gps_file} -p {input.perm_file} -a {wildcards.effect_blocks_A} -b {wildcards.effect_blocks_B} -o {output}"

rule compute_hoeffdings_for_sim_pair:
    input:
        sum_stats_file = "results/simgwas/simulated_sum_stats/whole_genome_sum_stats/{no_reps}_reps/randomised/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/window_{window}_step_{step}/seed_{seed}_pruned_sum_stats_tags_{tag_A}-{tag_B}.tsv.gz"
    output:
        "results/hoeffdings/simgwas/{no_reps}_reps/randomised/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/window_{window}_step_{step}/seed_{seed,\d+}_tags_{tag_A,\d+}-{tag_B,\d+}_hoeffdings.tsv"
    params:
        a_colname = lambda wildcards: f"p.{wildcards.tag_A}",
        b_colname = lambda wildcards: f"p.{wildcards.tag_B}"
    priority: 1
    group: "ldsc_hoeffding_sumher_gps_sans_permutation"
    resources:
        runtime = 2
    shell:
        """
        Rscript workflow/scripts/compute_hoeffdings.R -i {input.sum_stats_file} -a {params.a_colname} -b {params.b_colname} -o {output} -nt 1
        sed -i 's/{params.a_colname}/{wildcards.effect_blocks_A}_{wildcards.shared_effect_blocks}/' {output}
        sed -i 's/{params.b_colname}/{wildcards.effect_blocks_B}_{wildcards.shared_effect_blocks}/' {output}
        """
