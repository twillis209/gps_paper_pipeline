rule compute_per_chr_gps_for_trait_pair:
    input:
        "results/simgwas/simulated_sum_stats/whole_genome_sum_stats/{no_reps}_reps/randomised/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/window_{window}_step_{step}/seed_{seed}_pruned_sum_stats_tags_{tag_A}-{tag_B}_per_chr/{chr}.tsv"
    output:
        "results/gps/simgwas/{no_reps}_reps/randomised/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/window_{window}_step_{step}/seed_{seed}_tags_{tag_A,\d+}-{tag_B,\d+}_per_chr/{chr}_{ecdf}_pert_{pert}_gps_value.tsv"
    params:
        epsilon_multiple = 2.0
    threads: lambda wildcards: 12 if wildcards.ecdf == 'naive' else 1
    resources:
        runtime = 20
    group: "gps"
    shell:
        "workflow/scripts/gps_cpp/build/apps/computeGpsCLI -i {input} -a P.A -b P.B -c {wildcards.effect_blocks_A} -d {wildcards.effect_blocks_B} -f {wildcards.ecdf} -n {threads} -p {wildcards.pert} -e {params.epsilon_multiple} -o {output}"

rule compute_per_chr_gps_with_intermediate_values:
    input:
        "results/simgwas/simulated_sum_stats/whole_genome_sum_stats/{no_reps}_reps/randomised/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/window_{window}_step_{step}/seed_{seed}_pruned_sum_stats_tags_{tag_A}-{tag_B}_per_chr/{chr}.tsv"
    output:
        "results/gps/simgwas/{no_reps}_reps/randomised/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/window_{window}_step_{step}/seed_{seed}_tags_{tag_A,\d+}-{tag_B,\d+}_per_chr/{chr}_{ecdf}_pert_{pert}_gps_intermediates.tsv",
    threads: 1
    resources:
        runtime = 30
    group: "gps"
    shell:
        "workflow/scripts/gps_cpp/build/apps/fitAndEvaluateEcdfsCLI -i {input} -a P.A -b P.B -f {wildcards.ecdf} -n {threads} -o {output}"

rule permute_per_chr_gps_for_trait_pair:
    input:
        "results/simgwas/simulated_sum_stats/whole_genome_sum_stats/{no_reps}_reps/randomised/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/window_{window}_step_{step}/seed_{seed}_pruned_sum_stats_tags_{tag_A}-{tag_B}_per_chr/{chr}.tsv"
    output:
        "results/gps/simgwas/{no_reps}_reps/randomised/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/window_{window}_step_{step}/seed_{seed}_tags_{tag_A,\d+}-{tag_B,\d+}_per_chr/{draws}_permutations/{chr}.tsv",
    params:
        a_colname = lambda wildcards: f"p.{wildcards.tag_A}",
        b_colname = lambda wildcards: f"p.{wildcards.tag_B}",
        no_of_pert_iterations = 1
    threads: 12
    resources:
        mem_mb = get_mem_mb,
        runtime = get_permute_time,
    group: "gps"
    shell:
        "workflow/scripts/gps_cpp/build/apps/permuteTraitsCLI -i {input.sum_stats_file} -o {output} -a {params.a_colname} -b {params.b_colname} -c {threads} -n {wildcards.draws} -p {params.no_of_pert_iterations}"

rule test_per_chr_gps:
    input:
        "results/gps/simgwas/400_reps/randomised/10000_10000_10000_10000/s400_s400_s0/window_1000kb_step_50/seed_48801_tags_1-2_per_chr/chr22_naive_pert_0_gps_value.tsv"
#        "results/gps/simgwas/400_reps/randomised/10000_10000_10000_10000/s400_s400_s0/window_1000kb_step_50/seed_48801_tags_1-2_per_chr/chr22_naive_pert_0_gps_value.tsv"
