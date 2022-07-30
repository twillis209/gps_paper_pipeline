from string import ascii_lowercase

tags = [str(x) for x in range(1, 401)]

include: "simgwas_randomised_functions.py"

rule combine_randomised_block_sum_stats_for_pair:
    input:
        block_files = ancient(lambda wildcards: get_randomised_block_files_for_pair(wildcards)[0]),
        a_block_files = ancient(lambda wildcards: get_randomised_block_files_for_pair(wildcards)[1]),
        b_block_files = ancient(lambda wildcards: get_randomised_block_files_for_pair(wildcards)[2])
    output:
        combined_sum_stats_A = temp("results/simgwas/simulated_sum_stats/whole_genome_sum_stats/{no_reps}_reps/randomised/{ncases_A,\d+}_{ncontrols_A,\d+}_{ncases_B,\d+}_{ncontrols_B,\d+}/{effect_blocks_A,[smlvh\d-]+}_{effect_blocks_B,[smlvh\d-]+}_{shared_effect_blocks,[smlvh\d-]+}/{effect_blocks_A}_seed_{seed,\d+}_sum_stats_A_tag_{tag_A,\d+}_of_{tag_A}-{tag_B,\d+}.tsv.gz"),
        combined_sum_stats_B = temp("results/simgwas/simulated_sum_stats/whole_genome_sum_stats/{no_reps}_reps/randomised/{ncases_A,\d+}_{ncontrols_A,\d+}_{ncases_B,\d+}_{ncontrols_B,\d+}/{effect_blocks_A,[smlvh\d-]+}_{effect_blocks_B,[smlvh\d-]+}_{shared_effect_blocks,[smlvh\d-]+}/{effect_blocks_B}_seed_{seed,\d+}_sum_stats_B_tag_{tag_B,\d+}_of_{tag_A,\d+}-{tag_B}.tsv.gz")
    log:
        log = "logs/combine_randomised_block_sum_stats_for_pair/{no_reps}_reps/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}_{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}_seed_{seed}_tags_{tag_A}-{tag_B}.log"
    params:
        no_reps = 400,
        uncomp_sum_stats_A = "results/simgwas/simulated_sum_stats/whole_genome_sum_stats/{no_reps}_reps/randomised/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/{effect_blocks_A}_seed_{seed}_sum_stats_A_tag_{tag_A}_of_{tag_A}-{tag_B}.tsv",
        uncomp_sum_stats_B = "results/simgwas/simulated_sum_stats/whole_genome_sum_stats/{no_reps}_reps/randomised/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/{effect_blocks_B}_seed_{seed}_sum_stats_B_tag_{tag_B}_of_{tag_A}-{tag_B}.tsv",
        tag_A = lambda wildcards: wildcards.tag_A,
        tag_B = lambda wildcards: wildcards.tag_B
    threads: 1
    resources:
        runtime = 180
    group: "combine_randomised_block_sum_stats_for_pair"
    benchmark:
        "results/benchmarks/combine_randomised_block_sum_stats_for_pair/{no_reps}_reps/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}_{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}_seed_{seed}_tags_{tag_A}-{tag_B}.txt"
    run:
        with open(log.log, 'w') as logfile:
            z_column_name_A = f"zsim.{wildcards.tag_A}"
            beta_column_name_A = f"betasim.{wildcards.tag_A}"
            p_column_name_A = f"p.{wildcards.tag_A}"

            header_string_A = "\t".join(["position", "a0", "a1", "id", "block", "TYPE", "EUR", z_column_name_A, beta_column_name_A, p_column_name_A, "ncases", "ncontrols", "chr", "block_effect_size"])

            z_column_index_A = 8 + int(wildcards.tag_A)
            beta_column_index_A = 8 + (params.no_reps * 2) + int(wildcards.tag_A)
            p_column_index_A = 8 + (params.no_reps * 3) + int(wildcards.tag_A)

            ncases_column_index = 8  +  (params.no_reps  *  4) + 2
            ncontrols_column_index = 8 + (params.no_reps * 4) + 3
            chr_column_index = 8 + (params.no_reps * 4) + 4
            block_effect_column_index = 8 + (params.no_reps * 4) + 5

            cut_string_A = f"1-7,{z_column_index_A},{beta_column_index_A},{p_column_index_A},{ncases_column_index},{ncontrols_column_index},{chr_column_index},{block_effect_column_index}"

            shell("echo -e \"{header_string_A}\" > {params.uncomp_sum_stats_A}")

            for x in input.a_block_files:
                shell("zcat {x} | cut -f{cut_string_A} >> {params.uncomp_sum_stats_A}")
                logfile.write(f"{x}\n")

            shell(f"gzip {params.uncomp_sum_stats_A}")

            logfile.write("\n")

            z_column_name_B = f"zsim.{wildcards.tag_A}"
            beta_column_name_B = f"betasim.{wildcards.tag_A}"
            p_column_name_B = f"p.{wildcards.tag_A}"

            header_string_B = "\t".join(["position", "a0", "a1", "id", "block", "TYPE", "EUR", z_column_name_B, beta_column_name_B, p_column_name_B, "ncases", "ncontrols", "chr", "block_effect_size"])

            z_column_index_B = 8 + int(wildcards.tag_B)
            beta_column_index_B = 8 + (params.no_reps * 2) + int(wildcards.tag_B)
            p_column_index_B = 8 + (params.no_reps * 3) + int(wildcards.tag_B)

            cut_string_B = f"1-7,{z_column_index_B},{beta_column_index_B},{p_column_index_B},{ncases_column_index},{ncontrols_column_index},{chr_column_index},{block_effect_column_index}"

            shell("echo -e \"{header_string_B}\" > {params.uncomp_sum_stats_B}")

            for x in input.b_block_files:
                shell("zcat {x} | cut -f{cut_string_B} >> {params.uncomp_sum_stats_B}")
                logfile.write(f"{x}\n")

            shell(f"gzip {params.uncomp_sum_stats_B}")

rule merge_randomised_simulated_sum_stats:
    input:
        sum_stats_file_A = "results/simgwas/simulated_sum_stats/whole_genome_sum_stats/{no_reps}_reps/randomised/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/{effect_blocks_A}_seed_{seed}_sum_stats_A_tag_{tag_A}_of_{tag_A}-{tag_B}.tsv.gz",
        sum_stats_file_B = "results/simgwas/simulated_sum_stats/whole_genome_sum_stats/{no_reps}_reps/randomised/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/{effect_blocks_B}_seed_{seed}_sum_stats_B_tag_{tag_B}_of_{tag_A}-{tag_B}.tsv.gz"
    output:
        temp("results/simgwas/simulated_sum_stats/whole_genome_sum_stats/{no_reps}_reps/randomised/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A,[smlvh\d-]+}_{effect_blocks_B,[smlvh\d-]+}_{shared_effect_blocks,[smlvh\d-]+}/seed_{seed,\d+}_merged_sum_stats_tags_{tag_A,\d+}-{tag_B,\d+}.tsv.gz")
    threads: 4
    params:
        file_A_stat_cols = lambda wildcards: f"p.{wildcards.tag_A}",
        file_B_stat_cols = lambda wildcards: f"p.{wildcards.tag_B}"
    resources:
        runtime = 5
    group: "ldsc_hoeffding_and_gps_sans_permutation"
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
        runtime = 5
    group: "ldsc_hoeffding_and_gps_sans_permutation"
    shell:
        "Rscript workflow/scripts/simgwas/prune_sim_sum_stats.R --sum_stats_file {input.sum_stats_file} --bim_file {input.bim_file} --prune_file {input.pruned_range_file} -o {output} -nt {threads}"
