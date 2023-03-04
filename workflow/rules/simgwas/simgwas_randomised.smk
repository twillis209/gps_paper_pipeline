tags = [str(x) for x in range(1, 401)]

include: "simgwas_randomised_functions.py"

localrules: count_lines_in_sum_stats, count_lines_in_pruned_sum_stats

rule tabulate_randomised_block_sum_stats_file_for_pair:
    input:
        # Much duplication of input files here but it's to avoid calling this time-consuming function more than once
        block_files = ancient(lambda wildcards: get_randomised_block_files_for_pair(wildcards))
    output:
        a_block_file = temp("results/simgwas/simulated_sum_stats/whole_genome_sum_stats/{no_reps}_reps/randomised/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/seed_{seed}_sum_stats_A_tags_{tag_A}-{tag_B}_block_files.txt"),
        b_block_file = temp("results/simgwas/simulated_sum_stats/whole_genome_sum_stats/{no_reps}_reps/randomised/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/seed_{seed}_sum_stats_B_tags_{tag_A}-{tag_B}_block_files.txt")
    params:
        no_of_blocks_in_genome = block_daf.shape[0]
    threads: 1
    resources:
        concurrent_sans_permute_jobs = 1,
        runtime = 2
    group: "ldsc_hoeffding_sumher_gps_sans_permutation"
    run:
        # Probably too clever by half
        a_block_files = input.block_files[-(2*params.no_of_blocks_in_genome):-params.no_of_blocks_in_genome]
        b_block_files = input.block_files[-params.no_of_blocks_in_genome:]

        with open(output.a_block_file, 'w') as a_out:
            a_out.writelines([f"{x}\n" for x in a_block_files])

        with open(output.b_block_file, 'w') as b_out:
            b_out.writelines([f"{x}\n" for x in b_block_files])

rule split_block_files_for_pair:
    input:
        a_block_file = ancient("results/simgwas/simulated_sum_stats/whole_genome_sum_stats/{no_reps}_reps/randomised/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/seed_{seed}_sum_stats_A_tags_{tag_A}-{tag_B}_block_files.txt"),
        b_block_file = ancient("results/simgwas/simulated_sum_stats/whole_genome_sum_stats/{no_reps}_reps/randomised/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/seed_{seed}_sum_stats_B_tags_{tag_A}-{tag_B}_block_files.txt")
    params:
        a_file_prefix = "results/simgwas/simulated_sum_stats/whole_genome_sum_stats/{no_reps}_reps/randomised/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/seed_{seed}_sum_stats_A_tags_{tag_A}-{tag_B}_split_files/a_file_",
        b_file_prefix = "results/simgwas/simulated_sum_stats/whole_genome_sum_stats/{no_reps}_reps/randomised/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/seed_{seed}_sum_stats_B_tags_{tag_A}-{tag_B}_split_files/b_file_",
        suffix = "-of-12.txt",
        # TODO handle this using the config set-scatter value, see issue #1779 and associated fix
        no_of_splits = 12
    output:
        a_files = temp(scatter.split_block_files_for_pair("results/simgwas/simulated_sum_stats/whole_genome_sum_stats/{{no_reps}}_reps/randomised/{{ncases_A}}_{{ncontrols_A}}_{{ncases_B}}_{{ncontrols_B}}/{{effect_blocks_A}}_{{effect_blocks_B}}_{{shared_effect_blocks}}/seed_{{seed}}_sum_stats_A_tags_{{tag_A}}-{{tag_B}}_split_files/a_file_{scatteritem}.txt")),
        b_files = temp(scatter.split_block_files_for_pair("results/simgwas/simulated_sum_stats/whole_genome_sum_stats/{{no_reps}}_reps/randomised/{{ncases_A}}_{{ncontrols_A}}_{{ncases_B}}_{{ncontrols_B}}/{{effect_blocks_A}}_{{effect_blocks_B}}_{{shared_effect_blocks}}/seed_{{seed}}_sum_stats_B_tags_{{tag_A}}-{{tag_B}}_split_files/b_file_{scatteritem}.txt"))
    resources:
        concurrent_sans_permute_jobs = 1,
        runtime = 2
    group: "ldsc_hoeffding_sumher_gps_sans_permutation"
    shell:
        """
        split --numeric-suffixes=1 -nl/{params.no_of_splits} {input.a_block_file} {params.a_file_prefix} --additional-suffix={params.suffix}
        split --numeric-suffixes=1 -nl/{params.no_of_splits} {input.b_block_file} {params.b_file_prefix} --additional-suffix={params.suffix}

        # TODO hard-coding bad
        for i in {{1..9}}; do
            mv {params.a_file_prefix}"0"$i"-of-12.txt"  {params.a_file_prefix}$i"-of-12.txt"
        done

        for i in {{1..9}}; do
            mv {params.b_file_prefix}"0"$i"-of-12.txt"  {params.b_file_prefix}$i"-of-12.txt"
        done
        """

rule cat_split_block_files:
    input:
        a_block_file = ancient("results/simgwas/simulated_sum_stats/whole_genome_sum_stats/{no_reps}_reps/randomised/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/seed_{seed}_sum_stats_A_tags_{tag_A}-{tag_B}_split_files/a_file_{scatteritem}.txt"),
        b_block_file = ancient("results/simgwas/simulated_sum_stats/whole_genome_sum_stats/{no_reps}_reps/randomised/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/seed_{seed}_sum_stats_B_tags_{tag_A}-{tag_B}_split_files/b_file_{scatteritem}.txt")
    output:
        combined_sum_stats_A = temp("results/simgwas/simulated_sum_stats/whole_genome_sum_stats/{no_reps}_reps/randomised/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/seed_{seed}_sum_stats_A_tags_{tag_A}-{tag_B}_split_files/a_stats_{scatteritem}.tsv.gz"),
        combined_sum_stats_B = temp("results/simgwas/simulated_sum_stats/whole_genome_sum_stats/{no_reps}_reps/randomised/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/seed_{seed}_sum_stats_B_tags_{tag_A}-{tag_B}_split_files/b_stats_{scatteritem}.tsv.gz")
    params:
        line_count_bound = 850000
    threads: 1
    resources:
        runtime = lambda wildcards, attempt: 80*attempt,
        mem_mb = get_mem_mb,
    retries: 1
    group: "ldsc_hoeffding_sumher_gps_sans_permutation"
    script: "../../scripts/simgwas/combine_randomised_block_sum_stats.R"

rule gather_split_block_files:
    input:
        a_files = ancient(gather.split_block_files_for_pair("results/simgwas/simulated_sum_stats/whole_genome_sum_stats/{{no_reps}}_reps/randomised/{{ncases_A}}_{{ncontrols_A}}_{{ncases_B}}_{{ncontrols_B}}/{{effect_blocks_A}}_{{effect_blocks_B}}_{{shared_effect_blocks}}/seed_{{seed}}_sum_stats_A_tags_{{tag_A}}-{{tag_B}}_split_files/a_stats_{scatteritem}.tsv.gz")),
        b_files = ancient(gather.split_block_files_for_pair("results/simgwas/simulated_sum_stats/whole_genome_sum_stats/{{no_reps}}_reps/randomised/{{ncases_A}}_{{ncontrols_A}}_{{ncases_B}}_{{ncontrols_B}}/{{effect_blocks_A}}_{{effect_blocks_B}}_{{shared_effect_blocks}}/seed_{{seed}}_sum_stats_B_tags_{{tag_A}}-{{tag_B}}_split_files/b_stats_{scatteritem}.tsv.gz"))
    output:
        combined_sum_stats_A = temp("results/simgwas/simulated_sum_stats/whole_genome_sum_stats/{no_reps}_reps/randomised/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/seed_{seed}_sum_stats_A_tag_{tag_A}_of_{tag_A}-{tag_B}.tsv.gz"),
        combined_sum_stats_B = temp("results/simgwas/simulated_sum_stats/whole_genome_sum_stats/{no_reps}_reps/randomised/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/seed_{seed}_sum_stats_B_tag_{tag_B}_of_{tag_A}-{tag_B}.tsv.gz")
    params:
        expected_line_count = 8998662
    threads: 12
    resources:
        mem_mb = get_mem_mb,
        concurrent_sans_permute_jobs = 1,
        runtime = 10
    group: "ldsc_hoeffding_sumher_gps_sans_permutation"
    script: "../../scripts/simgwas/gather_split_block_files.R"

rule merge_randomised_simulated_sum_stats:
    input:
        sum_stats_file_A = "results/simgwas/simulated_sum_stats/whole_genome_sum_stats/{no_reps}_reps/randomised/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/seed_{seed}_sum_stats_A_tag_{tag_A}_of_{tag_A}-{tag_B}.tsv.gz",
        sum_stats_file_B = "results/simgwas/simulated_sum_stats/whole_genome_sum_stats/{no_reps}_reps/randomised/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/seed_{seed}_sum_stats_B_tag_{tag_B}_of_{tag_A}-{tag_B}.tsv.gz"
    output:
        temp("results/simgwas/simulated_sum_stats/whole_genome_sum_stats/{no_reps}_reps/randomised/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/seed_{seed}_merged_sum_stats_tags_{tag_A}-{tag_B}.tsv.gz")
    threads: 4
    params:
        file_A_stat_cols = lambda wildcards: f"p.{wildcards.tag_A}",
        file_B_stat_cols = lambda wildcards: f"p.{wildcards.tag_B}"
    resources:
        runtime = 5
    priority: 1
    group: "ldsc_hoeffding_sumher_gps_sans_permutation"
    script: "../../scripts/simgwas/merge_sim_sum_stats.R"

rule prune_merged_randomised_simulated_sum_stats:
    input:
        bim_file = "resources/1000g/euro/qc/all/all/all.bim",
        pruned_range_file = "resources/1000g/euro/qc/all/all/ranges/prune/window_{window}_step_{step}_r2_{r2}/all.prune.in",
        sum_stats_file = "results/simgwas/simulated_sum_stats/whole_genome_sum_stats/{no_reps}_reps/randomised/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/seed_{seed}_merged_sum_stats_tags_{tag_A}-{tag_B}.tsv.gz"
    output:
        temp("results/simgwas/simulated_sum_stats/whole_genome_sum_stats/{no_reps}_reps/randomised/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/window_{window}_step_{step}_r2_{r2}/seed_{seed}_pruned_sum_stats_tags_{tag_A}-{tag_B}.tsv.gz")
    threads: 4
    resources:
        runtime = 5
    group: "ldsc_hoeffding_sumher_gps_sans_permutation"
    shell:
        "Rscript workflow/scripts/simgwas/prune_sim_sum_stats.R --sum_stats_file {input.sum_stats_file} --bim_file {input.bim_file} --prune_file {input.pruned_range_file} -o {output} -nt {threads}"

rule unzip_pruned_merged_randomised_simulated_sum_stats:
    input:
        pruned_sum_stats_file = "results/simgwas/simulated_sum_stats/whole_genome_sum_stats/{no_reps}_reps/randomised/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/window_{window}_step_{step}_r2_{r2}/seed_{seed}_pruned_sum_stats_tags_{tag_A}-{tag_B}.tsv.gz"
    output:
        temp("results/simgwas/simulated_sum_stats/whole_genome_sum_stats/{no_reps}_reps/randomised/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/window_{window}_step_{step}_r2_{r2}/seed_{seed}_pruned_sum_stats_tags_{tag_A}-{tag_B}.tsv")
    threads: 1
    resources:
        runtime = 5
    priority: 1
    group: "ldsc_hoeffding_sumher_gps_sans_permutation"
    shell:
        "gunzip -c {input} >{output}"

rule write_out_per_chromosome_merged_randomised_simulated_sum_stats:
    input:
        "results/simgwas/simulated_sum_stats/whole_genome_sum_stats/{no_reps}_reps/randomised/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/window_{window}_step_{step}_r2_{r2}/seed_{seed}_pruned_sum_stats_tags_{tag_A}-{tag_B}.tsv"
    output:
        [f"results/simgwas/simulated_sum_stats/whole_genome_sum_stats/{{no_reps}}_reps/randomised/{{ncases_A}}_{{ncontrols_A}}_{{ncases_B}}_{{ncontrols_B}}/{{effect_blocks_A}}_{{effect_blocks_B}}_{{shared_effect_blocks}}/window_{{window}}_step_{{step}}_r2_{{r2}}/seed_{{seed}}_pruned_sum_stats_tags_{{tag_A}}-{{tag_B}}_per_chr/chr{x}.tsv" for x in range(1, 23)]
    params:
        chr_col = 'CHR38'
    group: "gwas"
    script: "../scripts/write_out_summary_statistics_per_chromosome.R"

rule count_lines_in_sum_stats:
    input:
        "results/simgwas/simulated_sum_stats/whole_genome_sum_stats/{no_reps}_reps/randomised/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/seed_{seed}_sum_stats_{label}_tag_{tag}_of_{tag_A}-{tag_B}.tsv.gz"
    output:
        "results/simgwas/simulated_sum_stats/whole_genome_sum_stats/{no_reps}_reps/randomised/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/seed_{seed}_sum_stats_{label,A|B}_tag_{tag}_of_{tag_A}-{tag_B}_linecount.txt"
    shell:
        "zcat {input} | tail -n +2 | wc -l >{output}"

rule count_lines_in_pruned_sum_stats:
    input:
        "results/simgwas/simulated_sum_stats/whole_genome_sum_stats/{no_reps}_reps/randomised/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/window_{window}_step_{step}_r2_{r2}/seed_{seed}_pruned_sum_stats_tags_{tag_A}-{tag_B}.tsv.gz"
    output:
        "results/simgwas/simulated_sum_stats/whole_genome_sum_stats/{no_reps}_reps/randomised/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/window_{window}_step_{step}_r2_{r2}/seed_{seed}_pruned_sum_stats_tags_{tag_A}-{tag_B}_linecount.txt"
    shell:
        "zcat {input} | tail -n +2 | wc -l >{output}"

rule run_count_jobs:
    input:
        "results/simgwas/simulated_sum_stats/whole_genome_sum_stats/400_reps/randomised/10000_10000_10000_10000/m50_m50_m0/seed_1_sum_stats_A_tag_1_of_1-2_linecount.txt",
        "results/simgwas/simulated_sum_stats/whole_genome_sum_stats/400_reps/randomised/10000_10000_10000_10000/m50_m50_m0/window_1000kb_step_50_r2_0_2/seed_1_pruned_sum_stats_tags_1-2_linecount.txt",
        "results/simgwas/simulated_sum_stats/whole_genome_sum_stats/400_reps/randomised/10000_10000_10000_10000/t1000_t1000_t0/seed_1_sum_stats_A_tag_1_of_1-2_linecount.txt",
        "results/simgwas/simulated_sum_stats/whole_genome_sum_stats/400_reps/randomised/10000_10000_10000_10000/t1000_t1000_t0/window_1000kb_step_50_r2_0_2/seed_1_pruned_sum_stats_tags_1-2_linecount.txt"

