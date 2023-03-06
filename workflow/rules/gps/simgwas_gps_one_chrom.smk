localrules: compute_li_gps_pvalue_for_chrom_for_pair, fit_gev_and_compute_gps_pvalue_for_chrom_for_pair

rule compute_gps_for_chrom_for_pair:
    input:
        sum_stats_file = "results/simgwas/simulated_sum_stats/per_chrom_sum_stats/{no_reps}_reps/randomised/{chr}/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/window_{window}_step_{step}_r2_{r2}/seed_{seed}_pruned_sum_stats_tags_{tag_A}-{tag_B}.tsv"
    output:
        temp("results/gps/simgwas/{no_reps}_reps/randomised/{chr}/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/window_{window}_step_{step}_r2_{r2}/seed_{seed}_tags_{tag_A}-{tag_B}_gps_value.tsv")
    params:
        a_colname = lambda wildcards: f"p.{wildcards.tag_A}",
        b_colname = lambda wildcards: f"p.{wildcards.tag_B}",
        no_of_pert_iterations = 1
    threads: 1
    priority: 1
    group: "one_chrom_analysis"
    shell:
        "workflow/scripts/gps_cpp/build/apps/computeGpsCLI -i {input.sum_stats_file} -a {params.a_colname} -b {params.b_colname} -c {wildcards.effect_blocks_A} -d {wildcards.effect_blocks_B} -p {params.no_of_pert_iterations} -n {threads} -f pp -o {output}"

rule permute_for_chrom_for_pair:
    input:
        sum_stats_file = "results/simgwas/simulated_sum_stats/per_chrom_sum_stats/{no_reps}_reps/randomised/{chr}/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/window_{window}_step_{step}_r2_{r2}/seed_{seed}_pruned_sum_stats_tags_{tag_A}-{tag_B}.tsv"
    output:
        "results/gps/simgwas/{no_reps}_reps/randomised/{chr}/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/window_{window}_step_{step}_r2_{r2}/{draws}_permutations/seed_{seed}_tags_{tag_A}-{tag_B}.tsv"
    params:
        a_colname = lambda wildcards: f"p.{wildcards.tag_A}",
        b_colname = lambda wildcards: f"p.{wildcards.tag_B}",
        no_of_pert_iterations = 1
    threads: 12
    resources:
        mem_mb = get_mem_mb,
        runtime = 10
    group: "permutation"
    priority: 1
    shell:
        "workflow/scripts/gps_cpp/build/apps/permuteTraitsCLI -i {input.sum_stats_file} -o {output} -a {params.a_colname} -b {params.b_colname} -c {threads} -n {wildcards.draws} -p {params.no_of_pert_iterations}"

rule compute_li_gps_pvalue_for_chrom_for_pair:
    input:
        "results/gps/simgwas/{no_reps}_reps/randomised/{chr}/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/window_{window}_step_{step}_r2_{r2}/seed_{seed}_tags_{tag_A}-{tag_B}_gps_value.tsv",
    output:
        "results/gps/simgwas/{no_reps}_reps/randomised/{chr}/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/window_{window}_step_{step}_r2_{r2}/seed_{seed}_tags_{tag_A}-{tag_B}_li_gps_pvalue.tsv"
    group: "one_chrom_analysis"
    script:
        "../../scripts/compute_li_gps_pvalue.R"

rule fit_gev_and_compute_gps_pvalue_for_chrom_for_pair:
    input:
        gps_file = "results/gps/simgwas/{no_reps}_reps/randomised/{chr}/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/window_{window}_step_{step}_r2_{r2}/seed_{seed}_tags_{tag_A}-{tag_B}_gps_value.tsv",
        perm_file = "results/gps/simgwas/{no_reps}_reps/randomised/{chr}/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/window_{window}_step_{step}_r2_{r2}/{draws}_permutations/seed_{seed}_tags_{tag_A}-{tag_B}.tsv"
    output:
        "results/gps/simgwas/{no_reps}_reps/randomised/{chr}/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/window_{window}_step_{step}_r2_{r2}/{draws}_permutations/seed_{seed}_tags_{tag_A}-{tag_B}_gps_pvalue.tsv"
    resources:
        runtime = 2
    group: "permutation"
    shell:
      "Rscript workflow/scripts/fit_gev_and_compute_gps_pvalue.R -g {input.gps_file} -p {input.perm_file} -a {wildcards.effect_blocks_A} -b {wildcards.effect_blocks_B} -o {output}"

rule compute_hoeffdings_for_chrom_for_pair:
    input:
        sum_stats_file = "results/simgwas/simulated_sum_stats/per_chrom_sum_stats/{no_reps}_reps/randomised/{chr}/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/window_{window}_step_{step}_r2_{r2}/seed_{seed}_pruned_sum_stats_tags_{tag_A}-{tag_B}.tsv.gz"
    output:
        "results/hoeffdings/simgwas/{no_reps}_reps/randomised/{chr}/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/window_{window}_step_{step}_r2_{r2}/seed_{seed}_tags_{tag_A}-{tag_B}_hoeffdings.tsv"
    params:
        a_colname = lambda wildcards: f"p.{wildcards.tag_A}",
        b_colname = lambda wildcards: f"p.{wildcards.tag_B}"
    priority: 1
    group: "one_chrom_analysis"
    resources:
        runtime = 2
    script: "../../scripts/simgwas/compute_hoeffdings.R"
