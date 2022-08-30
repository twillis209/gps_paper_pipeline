tags = [str(x) for x in range(1, 401)]

include: "simgwas_randomised_functions.py"

rule tabulate_randomised_block_sum_stats_file_for_pair:
    input:
        # Much duplication of input files here but it's to avoid calling this time-consuming function more than once
        block_files = ancient(lambda wildcards: get_randomised_block_files_for_pair(wildcards))
    output:
        a_block_file = temp("results/simgwas/simulated_sum_stats/whole_genome_sum_stats/{no_reps}_reps/randomised/{ncases_A,\d+}_{ncontrols_A,\d+}_{ncases_B,\d+}_{ncontrols_B,\d+}/{effect_blocks_A,[smlvh\d-]+}_{effect_blocks_B,[smlvh\d-]+}_{shared_effect_blocks,[smlvh\d-]+}/seed_{seed,\d+}_sum_stats_A_tags_{tag_A}-{tag_B,\d+}_block_files.txt"),
        b_block_file = temp("results/simgwas/simulated_sum_stats/whole_genome_sum_stats/{no_reps}_reps/randomised/{ncases_A,\d+}_{ncontrols_A,\d+}_{ncases_B,\d+}_{ncontrols_B,\d+}/{effect_blocks_A,[smlvh\d-]+}_{effect_blocks_B,[smlvh\d-]+}_{shared_effect_blocks,[smlvh\d-]+}/seed_{seed,\d+}_sum_stats_B_tags_{tag_A}-{tag_B,\d+}_block_files.txt")
    params:
        no_of_blocks_in_genome = block_daf.shape[0]
    threads: 1
    group: "tabulate_and_combine_block_files"
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
        a_block_file = "results/simgwas/simulated_sum_stats/whole_genome_sum_stats/{no_reps}_reps/randomised/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/seed_{seed}_sum_stats_A_tags_{tag_A}-{tag_B}_block_files.txt",
        b_block_file = "results/simgwas/simulated_sum_stats/whole_genome_sum_stats/{no_reps}_reps/randomised/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/seed_{seed}_sum_stats_B_tags_{tag_A}-{tag_B}_block_files.txt"
    params:
        a_file_prefix = "results/simgwas/simulated_sum_stats/whole_genome_sum_stats/{no_reps}_reps/randomised/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/seed_{seed}_sum_stats_A_tags_{tag_A}-{tag_B}_split_files/a_file_",
        b_file_prefix = "results/simgwas/simulated_sum_stats/whole_genome_sum_stats/{no_reps}_reps/randomised/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/seed_{seed}_sum_stats_B_tags_{tag_A}-{tag_B}_split_files/b_file_",
        suffix = "-of-12.txt",
        # TODO handle this using the config set-scatter value, see issue #1779 and associated fix
        no_of_splits = 12
    output:
        a_files = temp(scatter.split_block_files_for_pair("results/simgwas/simulated_sum_stats/whole_genome_sum_stats/{{no_reps}}_reps/randomised/{{ncases_A,\d+}}_{{ncontrols_A,\d+}}_{{ncases_B,\d+}}_{{ncontrols_B,\d+}}/{{effect_blocks_A,[smlvh\d-]+}}_{{effect_blocks_B,[smlvh\d-]+}}_{{shared_effect_blocks,[smlvh\d-]+}}/seed_{{seed,\d+}}_sum_stats_A_tags_{{tag_A}}-{{tag_B,\d+}}_split_files/a_file_{scatteritem}.txt")),
        b_files = temp(scatter.split_block_files_for_pair("results/simgwas/simulated_sum_stats/whole_genome_sum_stats/{{no_reps}}_reps/randomised/{{ncases_A,\d+}}_{{ncontrols_A,\d+}}_{{ncases_B,\d+}}_{{ncontrols_B,\d+}}/{{effect_blocks_A,[smlvh\d-]+}}_{{effect_blocks_B,[smlvh\d-]+}}_{{shared_effect_blocks,[smlvh\d-]+}}/seed_{{seed,\d+}}_sum_stats_B_tags_{{tag_A}}-{{tag_B,\d+}}_split_files/b_file_{scatteritem}.txt"))
    group: "tabulate_and_combine_block_files"
    shell:
        """
        # TODO remove leading 0 for single digits
        split --numeric-suffixes=1 -n {params.no_of_splits} {input.a_block_file} {params.a_file_prefix} --additional-suffix={params.suffix}
        split --numeric-suffixes=1 -n {params.no_of_splits} {input.b_block_file} {params.b_file_prefix} --additional-suffix={params.suffix}

        # TODO remove leading 0 for single digits
        # specify "false name" with param?

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
        a_block_file = "results/simgwas/simulated_sum_stats/whole_genome_sum_stats/{no_reps}_reps/randomised/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/seed_{seed}_sum_stats_A_tags_{tag_A}-{tag_B}_split_files/a_file_{scatteritem}.txt",
        b_block_file = "results/simgwas/simulated_sum_stats/whole_genome_sum_stats/{no_reps}_reps/randomised/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/seed_{seed}_sum_stats_B_tags_{tag_A}-{tag_B}_split_files/b_file_{scatteritem}.txt"
    output:
        combined_sum_stats_A = temp("results/simgwas/simulated_sum_stats/whole_genome_sum_stats/{no_reps}_reps/randomised/{ncases_A,\d+}_{ncontrols_A,\d+}_{ncases_B,\d+}_{ncontrols_B,\d+}/{effect_blocks_A,[smlvh\d-]+}_{effect_blocks_B,[smlvh\d-]+}_{shared_effect_blocks,[smlvh\d-]+}/seed_{seed,\d+}_sum_stats_A_tags_{tag_A}-{tag_B,\d+}_split_files/a_stats_{scatteritem}.tsv.gz"),
        combined_sum_stats_B = temp("results/simgwas/simulated_sum_stats/whole_genome_sum_stats/{no_reps}_reps/randomised/{ncases_A,\d+}_{ncontrols_A,\d+}_{ncases_B,\d+}_{ncontrols_B,\d+}/{effect_blocks_A,[smlvh\d-]+}_{effect_blocks_B,[smlvh\d-]+}_{shared_effect_blocks,[smlvh\d-]+}/seed_{seed,\d+}_sum_stats_B_tags_{tag_A}-{tag_B,\d+}_split_files/b_stats_{scatteritem}.tsv.gz")
    params:
        uncomp_sum_stats_A = "results/simgwas/simulated_sum_stats/whole_genome_sum_stats/{no_reps}_reps/randomised/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/seed_{seed}_sum_stats_A_tags_{tag_A}-{tag_B}_split_files/a_stats_{scatteritem}.tsv",
        uncomp_sum_stats_B = "results/simgwas/simulated_sum_stats/whole_genome_sum_stats/{no_reps}_reps/randomised/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/seed_{seed}_sum_stats_B_tags_{tag_A}-{tag_B}_split_files/b_stats_{scatteritem}.tsv"
    group: "tabulate_and_combine_block_files"
    resources:
        runtime = 270,
        mem_mb = get_mem_mb,
        tmpdir = 'tmp'
    script: "../../scripts/simgwas/combine_randomised_block_sum_stats.py"

rule gather_split_block_files:
    input:
        a_files = gather.split_block_files_for_pair("results/simgwas/simulated_sum_stats/whole_genome_sum_stats/{{no_reps}}_reps/randomised/{{ncases_A}}_{{ncontrols_A}}_{{ncases_B}}_{{ncontrols_B}}/{{effect_blocks_A}}_{{effect_blocks_B}}_{{shared_effect_blocks}}/seed_{{seed}}_sum_stats_A_tags_{{tag_A}}-{{tag_B}}_split_files/a_stats_{scatteritem}.tsv.gz"),
        b_files = gather.split_block_files_for_pair("results/simgwas/simulated_sum_stats/whole_genome_sum_stats/{{no_reps}}_reps/randomised/{{ncases_A}}_{{ncontrols_A}}_{{ncases_B}}_{{ncontrols_B}}/{{effect_blocks_A}}_{{effect_blocks_B}}_{{shared_effect_blocks}}/seed_{{seed}}_sum_stats_B_tags_{{tag_A}}-{{tag_B}}_split_files/b_stats_{scatteritem}.tsv.gz")
    output:
        combined_sum_stats_A = temp("results/simgwas/simulated_sum_stats/whole_genome_sum_stats/{no_reps}_reps/randomised/{ncases_A,\d+}_{ncontrols_A,\d+}_{ncases_B,\d+}_{ncontrols_B,\d+}/{effect_blocks_A,[smlvh\d-]+}_{effect_blocks_B,[smlvh\d-]+}_{shared_effect_blocks,[smlvh\d-]+}/seed_{seed,\d+}_sum_stats_A_tag_{tag_A,\d+}_of_{tag_A}-{tag_B,\d+}.tsv.gz"),
        combined_sum_stats_B = temp("results/simgwas/simulated_sum_stats/whole_genome_sum_stats/{no_reps}_reps/randomised/{ncases_A,\d+}_{ncontrols_A,\d+}_{ncases_B,\d+}_{ncontrols_B,\d+}/{effect_blocks_A,[smlvh\d-]+}_{effect_blocks_B,[smlvh\d-]+}_{shared_effect_blocks,[smlvh\d-]+}/seed_{seed,\d+}_sum_stats_B_tag_{tag_B,\d+}_of_{tag_A,\d+}-{tag_B}.tsv.gz")
    params:
        uncomp_sum_stats_A = temp("results/simgwas/simulated_sum_stats/whole_genome_sum_stats/{no_reps}_reps/randomised/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/seed_{seed}_sum_stats_A_tag_{tag_A}_of_{tag_A}-{tag_B}.tsv"),
        uncomp_sum_stats_B = temp("results/simgwas/simulated_sum_stats/whole_genome_sum_stats/{no_reps}_reps/randomised/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/seed_{seed}_sum_stats_B_tag_{tag_B}_of_{tag_A}-{tag_B}.tsv")
    threads: 1
    group: "tabulate_and_combine_block_files"
    run:
        z_column_name_A = f"zsim.{wildcards.tag_A}"
        beta_column_name_A = f"betasim.{wildcards.tag_A}"
        p_column_name_A = f"p.{wildcards.tag_A}"

        for x in a_files:
            shell("zcat {x} >> {params.uncomp_sum_stats_A}")

        header_string_A = "\t".join(["position", "a0", "a1", "id", "block", "TYPE", "EUR", z_column_name_A, beta_column_name_A, p_column_name_A, "ncases", "ncontrols", "chr", "block_effect_size"])

        shell("echo -e \"{header_string_A}\" > {params.uncomp_sum_stats_A}")

        shell("gzip {params.uncomp_sum_stats_A}")

        z_column_name_B = f"zsim.{wildcards.tag_B}"
        beta_column_name_B = f"betasim.{wildcards.tag_B}"
        p_column_name_B = f"p.{wildcards.tag_B}"

        for x in b_files:
            shell("zcat {x} >> {params.uncomp_sum_stats_B}")

        header_string_B = "\t".join(["position", "a0", "a1", "id", "block", "TYPE", "EUR", z_column_name_B, beta_column_name_B, p_column_name_B, "ncases", "ncontrols", "chr", "block_effect_size"])

        shell("echo -e \"{header_string_B}\" > {params.uncomp_sum_stats_B}")

        shell("gzip {params.uncomp_sum_stats_B}")

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
