tags = [str(x) for x in range(1, 401)]

include: "simgwas_randomised_functions.py"

rule tabulate_randomised_block_sum_stats_file_for_pair:
    input:
        # Much duplication of input files here but it's to avoid calling this time-consuming function more than once
        block_files = ancient(lambda wildcards: get_randomised_block_files_for_pair(wildcards))
    output:
        a_block_file = "results/simgwas/simulated_sum_stats/whole_genome_sum_stats/{no_reps}_reps/randomised/{ncases_A,\d+}_{ncontrols_A,\d+}_{ncases_B,\d+}_{ncontrols_B,\d+}/{effect_blocks_A,[smlvh\d-]+}_{effect_blocks_B,[smlvh\d-]+}_{shared_effect_blocks,[smlvh\d-]+}/seed_{seed,\d+}_sum_stats_A_tags_{tag_A}-{tag_B,\d+}_files.txt",
        b_block_file = "results/simgwas/simulated_sum_stats/whole_genome_sum_stats/{no_reps}_reps/randomised/{ncases_A,\d+}_{ncontrols_A,\d+}_{ncases_B,\d+}_{ncontrols_B,\d+}/{effect_blocks_A,[smlvh\d-]+}_{effect_blocks_B,[smlvh\d-]+}_{shared_effect_blocks,[smlvh\d-]+}/seed_{seed,\d+}_sum_stats_B_tags_{tag_A}-{tag_B,\d+}_files.txt"
    params:
        no_of_blocks_in_genome = block_daf.shape[0]
    threads: 1
    group: "ldsc_hoeffding_sumher_gps_sans_permutation"
    run:
        # Probably too clever by half
        a_block_files = input.block_files[-(2*params.no_of_blocks_in_genome):-params.no_of_blocks_in_genome]
        b_block_files = input.block_files[-params.no_of_blocks_in_genome:]

        with open(output.a_block_file, 'w') as a_out:
            a_out.writelines([f"{x}\n" for x in a_block_files])

        with open(output.b_block_file, 'w') as b_out:
            b_out.writelines([f"{x}\n" for x in b_block_files])

rule combine_randomised_block_sum_stats_for_pair:
    input:
        a_block_file = ancient("results/simgwas/simulated_sum_stats/whole_genome_sum_stats/{no_reps}_reps/randomised/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/seed_{seed}_sum_stats_A_tags_{tag_A}-{tag_B}_files.txt"),
        b_block_file = ancient("results/simgwas/simulated_sum_stats/whole_genome_sum_stats/{no_reps}_reps/randomised/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/seed_{seed}_sum_stats_B_tags_{tag_A}-{tag_B}_files.txt")
    output:
        combined_sum_stats_A = temp("results/simgwas/simulated_sum_stats/whole_genome_sum_stats/{no_reps}_reps/randomised/{ncases_A,\d+}_{ncontrols_A,\d+}_{ncases_B,\d+}_{ncontrols_B,\d+}/{effect_blocks_A,[smlvh\d-]+}_{effect_blocks_B,[smlvh\d-]+}_{shared_effect_blocks,[smlvh\d-]+}/seed_{seed,\d+}_sum_stats_A_tag_{tag_A,\d+}_of_{tag_A}-{tag_B,\d+}.tsv.gz"),
        combined_sum_stats_B = temp("results/simgwas/simulated_sum_stats/whole_genome_sum_stats/{no_reps}_reps/randomised/{ncases_A,\d+}_{ncontrols_A,\d+}_{ncases_B,\d+}_{ncontrols_B,\d+}/{effect_blocks_A,[smlvh\d-]+}_{effect_blocks_B,[smlvh\d-]+}_{shared_effect_blocks,[smlvh\d-]+}/seed_{seed,\d+}_sum_stats_B_tag_{tag_B,\d+}_of_{tag_A,\d+}-{tag_B}.tsv.gz")
    params:
        uncomp_sum_stats_A = "results/simgwas/simulated_sum_stats/whole_genome_sum_stats/{no_reps}_reps/randomised/{ncases_A,\d+}_{ncontrols_A,\d+}_{ncases_B,\d+}_{ncontrols_B,\d+}/{effect_blocks_A,[smlvh\d-]+}_{effect_blocks_B,[smlvh\d-]+}_{shared_effect_blocks,[smlvh\d-]+}/seed_{seed,\d+}_sum_stats_A_tag_{tag_A,\d+}_of_{tag_A}-{tag_B,\d+}.tsv",
        uncomp_sum_stats_B = "results/simgwas/simulated_sum_stats/whole_genome_sum_stats/{no_reps}_reps/randomised/{ncases_A,\d+}_{ncontrols_A,\d+}_{ncases_B,\d+}_{ncontrols_B,\d+}/{effect_blocks_A,[smlvh\d-]+}_{effect_blocks_B,[smlvh\d-]+}_{shared_effect_blocks,[smlvh\d-]+}/seed_{seed,\d+}_sum_stats_B_tag_{tag_B,\d+}_of_{tag_A,\d+}-{tag_B}.tsv"
    log:
        log = "logs/combine_randomised_block_sum_stats_for_pair/{no_reps}_reps/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}_{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}_seed_{seed}_tags_{tag_A}-{tag_B}.log"
    threads: 1
    resources:
        runtime = 270,
        mem_mb = get_mem_mb,
        tmpdir = 'tmp'
    group: "ldsc_hoeffding_sumher_gps_sans_permutation"
    benchmark:
        "results/benchmarks/combine_randomised_block_sum_stats_for_pair/{no_reps}_reps/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}_{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}_seed_{seed}_tags_{tag_A}-{tag_B}.txt"
    script: "../../scripts/simgwas/combine_randomised_block_sum_stats.py"

rule merge_randomised_simulated_sum_stats:
    input:
        sum_stats_file_A = "results/simgwas/simulated_sum_stats/whole_genome_sum_stats/{no_reps}_reps/randomised/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/seed_{seed}_sum_stats_A_tag_{tag_A}_of_{tag_A}-{tag_B}.tsv.gz",
        sum_stats_file_B = "results/simgwas/simulated_sum_stats/whole_genome_sum_stats/{no_reps}_reps/randomised/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/seed_{seed}_sum_stats_B_tag_{tag_B}_of_{tag_A}-{tag_B}.tsv.gz"
    output:
        temp("results/simgwas/simulated_sum_stats/whole_genome_sum_stats/{no_reps}_reps/randomised/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A,[smlvh\d-]+}_{effect_blocks_B,[smlvh\d-]+}_{shared_effect_blocks,[smlvh\d-]+}/seed_{seed,\d+}_merged_sum_stats_tags_{tag_A,\d+}-{tag_B,\d+}.tsv.gz")
    threads: 4
    params:
        file_A_stat_cols = lambda wildcards: f"p.{wildcards.tag_A}",
        file_B_stat_cols = lambda wildcards: f"p.{wildcards.tag_B}"
    resources:
        runtime = 10
    group: "ldsc_hoeffding_sumher_gps_sans_permutation"
    shell:
        "Rscript workflow/scripts/simgwas/merge_sim_sum_stats.R --sum_stats_file_A {input.sum_stats_file_A} --sum_stats_file_B {input.sum_stats_file_B} --file_A_stat_cols {params.file_A_stat_cols} --file_B_stat_cols {params.file_B_stat_cols} -o {output} -nt {threads}"

rule prune_merged_randomised_simulated_sum_stats:
    input:
        bim_file = "resources/1000g/euro/qc/chr1-22_qc.bim",
        pruned_range_file = "resources/plink_ranges/1000g/pruned_ranges/window_{window}_step_{step}/all.prune.in",
        sum_stats_file = "results/simgwas/simulated_sum_stats/whole_genome_sum_stats/{no_reps}_reps/randomised/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/seed_{seed}_merged_sum_stats_tags_{tag_A}-{tag_B}.tsv.gz"
    output:
        temp("results/simgwas/simulated_sum_stats/whole_genome_sum_stats/{no_reps}_reps/randomised/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A,[smlvh\d-]+}_{effect_blocks_B,[smlvh\d-]+}_{shared_effect_blocks,[smlvh\d-]+}/window_{window}_step_{step}/seed_{seed,\d+}_pruned_sum_stats_tags_{tag_A,\d+}-{tag_B,\d+}.tsv.gz")
    threads: 4
    resources:
        runtime = 10
    group: "ldsc_hoeffding_sumher_gps_sans_permutation"
    shell:
        "Rscript workflow/scripts/simgwas/prune_sim_sum_stats.R --sum_stats_file {input.sum_stats_file} --bim_file {input.bim_file} --prune_file {input.pruned_range_file} -o {output} -nt {threads}"

rule unzip_pruned_merged_randomised_simulated_sum_stats:
    input:
        pruned_sum_stats_file = "results/simgwas/simulated_sum_stats/whole_genome_sum_stats/{no_reps}_reps/randomised/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/window_{window}_step_{step}/seed_{seed}_pruned_sum_stats_tags_{tag_A}-{tag_B}.tsv.gz"
    output:
        temp("results/simgwas/simulated_sum_stats/whole_genome_sum_stats/{no_reps}_reps/randomised/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A,[smlvh\d-]+}_{effect_blocks_B,[smlvh\d-]+}_{shared_effect_blocks,[smlvh\d-]+}/window_{window}_step_{step}/seed_{seed,\d+}_pruned_sum_stats_tags_{tag_A,\d+}-{tag_B,\d+}.tsv")
    threads: 1
    resources:
        runtime = 10
    group: "ldsc_hoeffding_sumher_gps_sans_permutation"
    shell:
        "gunzip -c {input} >{output}"
