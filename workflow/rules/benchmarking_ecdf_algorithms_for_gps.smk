rule benchmark_naive_ecdf_algorithm_for_gps_using_sim_pair:
    input:
        "results/simgwas/simulated_sum_stats/whole_genome_sum_stats/{no_reps}_reps/randomised/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/window_{window}_step_{step}/seed_{seed}_pruned_sum_stats_tags_{tag_A}-{tag_B}.tsv"
    output:
        result_file = "benchmarks/results/gps/simgwas/{no_reps}_reps/randomised/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/window_{window}_step_{step}/pert_0/naive/seed_{seed}_pruned_sum_stats_tags_{tag_A}-{tag_B}_gps_value_rep_{rep}.tsv",
        timing_file = "benchmarks/results/gps/simgwas/{no_reps}_reps/randomised/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/window_{window}_step_{step}/pert_0/naive/seed_{seed}_pruned_sum_stats_tags_{tag_A}-{tag_B}_gps_value_rep_{rep}.benchmark.txt"
    group: "benchmark_naive_ecdf_algorithm"
    params:
        a_colname = lambda wildcards: f"p.{wildcards.tag_A}",
        b_colname = lambda wildcards: f"p.{wildcards.tag_B}"
    threads: 1
    resources:
        runtime = 90
    priority: 1
    shell:
        "workflow/scripts/gps_cpp/build/apps/computeGpsCLI -i {input} -a {params.a_colname} -b {params.b_colname} -c {wildcards.effect_blocks_A} -d {wildcards.effect_blocks_B} -p 0 -f naive -n {threads} -o {output.result_file} -j {output.timing_file}"

rule benchmark_fast_ecdf_algorithm_for_gps_using_sim_pair:
    input:
        "results/simgwas/simulated_sum_stats/whole_genome_sum_stats/{no_reps}_reps/randomised/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/window_{window}_step_{step}/seed_{seed}_pruned_sum_stats_tags_{tag_A}-{tag_B}.tsv"
    output:
        result_file = "benchmarks/results/gps/simgwas/{no_reps}_reps/randomised/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/window_{window}_step_{step}/pert_{no_of_pert_iterations}/{algorithm,pp|lw}/seed_{seed}_pruned_sum_stats_tags_{tag_A}-{tag_B}_gps_value_rep_{rep}.tsv",
        timing_file = "benchmarks/results/gps/simgwas/{no_reps}_reps/randomised/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/window_{window}_step_{step}/pert_{no_of_pert_iterations}/{algorithm,pp|lw}/seed_{seed}_pruned_sum_stats_tags_{tag_A}-{tag_B}_gps_value_rep_{rep}.benchmark.txt"
    group: "benchmark_fast_ecdf_algorithm"
    params:
        a_colname = lambda wildcards: f"p.{wildcards.tag_A}",
        b_colname = lambda wildcards: f"p.{wildcards.tag_B}"
    threads: 1
    resources:
        runtime = 5
    priority: 1
    shell:
        "workflow/scripts/gps_cpp/build/apps/computeGpsCLI -i {input} -a {params.a_colname} -b {params.b_colname} -c {wildcards.effect_blocks_A} -d {wildcards.effect_blocks_B} -p {wildcards.no_of_pert_iterations} -f {wildcards.algorithm} -n {threads} -o {output.result_file} -j {output.timing_file}"

rule benchmark_ecdf_algorithms:
    input:
        [f"benchmarks/results/gps/simgwas/400_reps/randomised/10000_10000_10000_10000/s400_s400_s0/window_1000kb_step_50/pert_0/naive/seed_1_pruned_sum_stats_tags_1-2_gps_value_rep_{x}.tsv" for x in range(1,101)],
        [f"benchmarks/results/gps/simgwas/400_reps/randomised/10000_10000_10000_10000/s400_s400_s0/window_1000kb_step_50/pert_1/pp/seed_1_pruned_sum_stats_tags_1-2_gps_value_rep_{x}.tsv" for x in range(1,101)],
        [f"benchmarks/results/gps/simgwas/400_reps/randomised/10000_10000_10000_10000/s400_s400_s0/window_1000kb_step_50/pert_100/lw/seed_1_pruned_sum_stats_tags_1-2_gps_value_rep_{x}.tsv" for x in range(1,101)]
