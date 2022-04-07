from string import ascii_lowercase

tags = list(ascii_lowercase[:20])
tag_pvalue_dict = dict(zip(tags, [f"p.{x}" for x in range(1,21)]))

def get_randomised_block_files(wildcards):
        block_files = []

        if wildcards.effect_blocks != 'null':
            block_match = re.match('([smlvh])(\d+)', wildcards.effect_blocks)

            if not block_match:
                raise ValueError("Invalid block format: %s" % wildcards.effect_blocks)

            effect = effect_size_dict[block_match.group(1)]
            no_of_blocks = int(block_match.group(2))

            effect_block_files = []

            for i in range(no_of_blocks):
                chrom = random.choice(list(block_dict.keys()))
                block_no = random.choice(block_dict[chrom])
                effect_block_files.append(f"results/simgwas/simulated_sum_stats/chr{chrom}/block_sum_stats/{effect}/{wildcards.ncases}_{wildcards.ncontrols}/block_{block_no}_sum_stats.tsv.gz")

            for chrom in block_dict.keys():
                for block_no in block_dict[chrom]:
                    if f"results/simgwas/simulated_sum_stats/chr{chrom}/block_sum_stats/{effect}/{wildcards.ncases}_{wildcards.ncontrols}/block_{block_no}_sum_stats.tsv.gz" not in effect_block_files:
                        block_files.append(f"results/simgwas/simulated_sum_stats/chr{chrom}/block_sum_stats/null/{wildcards.ncases}_{wildcards.ncontrols}/block_{block_no}_sum_stats.tsv.gz")

            block_files += effect_block_files
        else:
            for chrom in block_dict.items():
                for block_no in block_dict[chrom]:
                    block_files.append(f"results/simgwas/simulated_sum_stats/chr{chrom}/block_sum_stats/null/{wildcards.ncases}_{wildcards.ncontrols}/block_{block_no}_sum_stats.tsv.gz")

        with open(f"results/simgwas/simulated_sum_stats/whole_genome_sum_stats/randomised_blocks_logs/{wildcards.ncases}_{wildcards.ncontrols}/{wildcards.effect_blocks}_sum_stats.log", 'w') as outfile:
            for x in block_files:
                outfile.write(f"{x}\n")

        return block_files

def get_randomised_block_files_for_pair(wildcards):
    random.seed(wildcards.seed)

    shared_chrom_block_nos = []

    if wildcards.shared_effect_blocks != 'null':
        block_match = re.match('([smlvh])(\d+)', wildcards.shared_effect_blocks)

        if not block_match:
            raise ValueError("Invalid block format: %s" % wildcards.shared_effect_blocks)

        effect = effect_size_dict[block_match.group(1)]
        no_of_shared_blocks = int(block_match.group(2))

        for i in range(no_of_shared_blocks):
            chrom = random.choice(list(block_dict.keys()))
            block_no = random.choice(block_dict[chrom])

            if (chrom, block_no) not in shared_chrom_block_nos:
                shared_chrom_block_nos.append((chrom, block_no))

    a_chrom_block_nos = []

    if wildcards.effect_blocks_A != 'null':
        block_match_a = re.match('([smlvh])(\d+)', wildcards.effect_blocks_A)

        if not block_match_a:
            raise ValueError("Invalid block format: %s" % wildcards.effect_blocks_A)

        effect_a = effect_size_dict[block_match_a.group(1)]
        no_of_blocks_a = int(block_match_a.group(2))

        i = 0

        while i < max(no_of_blocks_a-no_of_shared_blocks, 0):
            chrom = random.choice(list(block_dict.keys()))
            block_no = random.choice(block_dict[chrom])

            if (chrom, block_no) not in shared_chrom_block_nos and (chrom, block_no) not in a_chrom_block_nos:
                a_chrom_block_nos.append((chrom, block_no))
                i += 1

    b_chrom_block_nos = []

    if wildcards.effect_blocks_B != 'null':
        block_match_b = re.match('([smlvh])(\d+)', wildcards.effect_blocks_B)

        if not block_match_b:
            raise ValueError("Invalid block format: %s" % wildcards.effect_blocks_B)

        effect_b = effect_size_dict[block_match_b.group(1)]
        no_of_blocks_b = int(block_match_b.group(2))

        i = 0

        while i < max(no_of_blocks_b-no_of_shared_blocks, 0):
            chrom = random.choice(list(block_dict.keys()))
            block_no = random.choice(block_dict[chrom])

            if (chrom, block_no) not in shared_chrom_block_nos and (chrom, block_no) not in a_chrom_block_nos and (chrom, block_no) not in b_chrom_block_nos:
                b_chrom_block_nos.append((chrom, block_no))
                i += 1

    block_files = []

    a_block_files = []
    b_block_files = []

    for x in shared_chrom_block_nos:
        block_files.append(f"results/simgwas/simulated_sum_stats/chr{x[0]}/block_sum_stats/{effect}/{wildcards.ncases_A}_{wildcards.ncontrols_A}/block_{x[1]}_sum_stats.tsv.gz")
        a_block_files.append(f"results/simgwas/simulated_sum_stats/chr{x[0]}/block_sum_stats/{effect}/{wildcards.ncases_A}_{wildcards.ncontrols_A}/block_{x[1]}_sum_stats.tsv.gz")
        b_block_files.append(f"results/simgwas/simulated_sum_stats/chr{x[0]}/block_sum_stats/{effect}/{wildcards.ncases_B}_{wildcards.ncontrols_B}/block_{x[1]}_sum_stats.tsv.gz")
        if wildcards.ncases_A != wildcards.ncases_B or wildcards.ncontrols_A != wildcards.ncontrols_B:
            block_files.append(f"results/simgwas/simulated_sum_stats/chr{x[0]}/block_sum_stats/{effect}/{wildcards.ncases_B}_{wildcards.ncontrols_B}/block_{x[1]}_sum_stats.tsv.gz")

    for x in a_chrom_block_nos:
        block_files.append(f"results/simgwas/simulated_sum_stats/chr{x[0]}/block_sum_stats/{effect_a}/{wildcards.ncases_A}_{wildcards.ncontrols_A}/block_{x[1]}_sum_stats.tsv.gz")
        a_block_files.append(f"results/simgwas/simulated_sum_stats/chr{x[0]}/block_sum_stats/{effect_a}/{wildcards.ncases_A}_{wildcards.ncontrols_A}/block_{x[1]}_sum_stats.tsv.gz")

    for x in b_chrom_block_nos:
        b_block_files.append(f"results/simgwas/simulated_sum_stats/chr{x[0]}/block_sum_stats/{effect_b}/{wildcards.ncases_B}_{wildcards.ncontrols_B}/block_{x[1]}_sum_stats.tsv.gz")
        if f"results/simgwas/simulated_sum_stats/chr{x[0]}/block_sum_stats/{effect_b}/{wildcards.ncases_B}_{wildcards.ncontrols_B}/block_{x[1]}_sum_stats.tsv.gz" not in block_files:
            block_files.append(f"results/simgwas/simulated_sum_stats/chr{x[0]}/block_sum_stats/{effect_b}/{wildcards.ncases_B}_{wildcards.ncontrols_B}/block_{x[1]}_sum_stats.tsv.gz")

    for chrom in block_dict.keys():
        for block_no in block_dict[chrom]:
            if (chrom, block_no) not in a_chrom_block_nos and (chrom, block_no) not in shared_chrom_block_nos:
                block_files.append(f"results/simgwas/simulated_sum_stats/chr{chrom}/block_sum_stats/null/{wildcards.ncases_A}_{wildcards.ncontrols_A}/block_{block_no}_sum_stats.tsv.gz")
                a_block_files.append(f"results/simgwas/simulated_sum_stats/chr{chrom}/block_sum_stats/null/{wildcards.ncases_A}_{wildcards.ncontrols_A}/block_{block_no}_sum_stats.tsv.gz")
            if (chrom, block_no) not in b_chrom_block_nos and (chrom, block_no) not in shared_chrom_block_nos:
                b_block_files.append(f"results/simgwas/simulated_sum_stats/chr{chrom}/block_sum_stats/null/{wildcards.ncases_B}_{wildcards.ncontrols_B}/block_{block_no}_sum_stats.tsv.gz")
                if f"results/simgwas/simulated_sum_stats/chr{chrom}/block_sum_stats/null/{wildcards.ncases_B}_{wildcards.ncontrols_B}/block_{block_no}_sum_stats.tsv.gz" not in block_files:
                    block_files.append(f"results/simgwas/simulated_sum_stats/chr{chrom}/block_sum_stats/null/{wildcards.ncases_B}_{wildcards.ncontrols_B}/block_{block_no}_sum_stats.tsv.gz")

    return (block_files, a_block_files, b_block_files)

rule combine_randomised_block_sum_stats_for_pair:
    input:
        block_files = lambda wildcards: get_randomised_block_files_for_pair(wildcards)[0],
        a_block_files = lambda wildcards: get_randomised_block_files_for_pair(wildcards)[1],
        b_block_files = lambda wildcards: get_randomised_block_files_for_pair(wildcards)[2]
    output:
        combined_sum_stats_A = temp("results/simgwas/simulated_sum_stats/whole_genome_sum_stats/randomised/{ncases_A,\d+}_{ncontrols_A,\d+}_{ncases_B,\d+}_{ncontrols_B,\d+}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/{effect_blocks_A}_seed_{seed}_sum_stats_A_tag_{tag_A}_of_{tag_A}{tag_B}.tsv.gz"),
        combined_sum_stats_B = temp("results/simgwas/simulated_sum_stats/whole_genome_sum_stats/randomised/{ncases_A,\d+}_{ncontrols_A,\d+}_{ncases_B,\d+}_{ncontrols_B,\d+}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/{effect_blocks_B}_seed_{seed}_sum_stats_B_tag_{tag_B}_of_{tag_A}{tag_B}.tsv.gz")
    log:
        log = "logs/combine_randomised_block_sum_stats_for_pair/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}_{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}_seed_{seed}_tags_{tag_A}{tag_B}.log"
    params:
        no_reps = 20,
        uncomp_sum_stats_A = "results/simgwas/simulated_sum_stats/whole_genome_sum_stats/randomised/{ncases_A,\d+}_{ncontrols_A,\d+}_{ncases_B,\d+}_{ncontrols_B,\d+}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/{effect_blocks_A}_seed_{seed}_sum_stats_A_tag_{tag_A}_of_{tag_A}{tag_B}.tsv",
        uncomp_sum_stats_B = "results/simgwas/simulated_sum_stats/whole_genome_sum_stats/randomised/{ncases_A,\d+}_{ncontrols_A,\d+}_{ncases_B,\d+}_{ncontrols_B,\d+}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/{effect_blocks_B}_seed_{seed}_sum_stats_B_tag_{tag_B}_of_{tag_A}{tag_B}.tsv",
        # NB: 1-indexed
        tag_index_A = lambda wildcards: tags.index(wildcards.tag_A)+1,
        tag_index_B = lambda wildcards: tags.index(wildcards.tag_B)+1
    run:
        with open(log.log, 'w') as logfile:
            p_column_name_A = f"p.{params.tag_index_A}"
            beta_column_name_A = f"betasim.{params.tag_index_A}"
            header_string_A = "\t".join(["position", "a0", "a1", "id", "block", "TYPE", "EUR", beta_column_name_A, p_column_name_A, "ncases", "ncontrols", "chr"])

            beta_column_index_A = 8+(params.no_reps*2)+params.tag_index_A
            p_column_index_A = 8+(params.no_reps*3)+params.tag_index_A

            ncases_column_index = 8+(params.no_reps*4)+2
            ncontrols_column_index = 8+(params.no_reps*4)+3
            chr_column_index = 8+(params.no_reps*4)+4

            cut_string_A = f"1-7,{beta_column_index_A},{p_column_index_A},{ncases_column_index},{ncontrols_column_index},{chr_column_index}"

            shell("echo -e \"{header_string_A}\" > {params.uncomp_sum_stats_A}")

            for x in input.a_block_files:
                shell("zcat {x} | cut -f{cut_string_A} >> {params.uncomp_sum_stats_A}")
                logfile.write(f"{x}\n")

            shell(f"gzip {params.uncomp_sum_stats_A}")

            logfile.write("\n")

            p_column_name_B = f"p.{params.tag_index_B}"
            beta_column_name_B = f"betasim.{params.tag_index_B}"
            header_string_B = "\t".join(["position", "a0", "a1", "id", "block", "TYPE", "EUR", beta_column_name_B, p_column_name_B, "ncases", "ncontrols", "chr"])

            beta_column_index_B = 8+(params.no_reps*2)+params.tag_index_B
            p_column_index_B = 8+(params.no_reps*3)+params.tag_index_B

            cut_string_B = f"1-7,{beta_column_index_B},{p_column_index_B},{ncases_column_index},{ncontrols_column_index},{chr_column_index}"

            shell("echo -e \"{header_string_B}\" > {params.uncomp_sum_stats_B}")

            for x in input.b_block_files:
                shell("zcat {x} | cut -f{cut_string_B} >> {params.uncomp_sum_stats_B}")
                logfile.write(f"{x}\n")

            shell(f"gzip {params.uncomp_sum_stats_B}")

# TODO Uh-oh: 8,998,661 in m1_seed_111_sum_stats_A.tsv.gz, 9,000,062 in seed_111_merged_sum_stats.tsv.gz
rule merge_randomised_simulated_sum_stats:
    input:
        sum_stats_file_A = "results/simgwas/simulated_sum_stats/whole_genome_sum_stats/randomised/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/{effect_blocks_A}_seed_{seed}_sum_stats_A_tag_{tag_A}_of_{tag_A}{tag_B}.tsv.gz",
        sum_stats_file_B = "results/simgwas/simulated_sum_stats/whole_genome_sum_stats/randomised/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/{effect_blocks_B}_seed_{seed}_sum_stats_B_tag_{tag_B}_of_{tag_A}{tag_B}.tsv.gz"
    output:
        temp("results/simgwas/simulated_sum_stats/whole_genome_sum_stats/randomised/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/seed_{seed}_merged_sum_stats_tags_{tag_A}{tag_B}.tsv.gz")
    threads: 4
    params:
        file_A_stat_cols = lambda wildcards: tag_pvalue_dict[wildcards.tag_A],
        file_B_stat_cols = lambda wildcards: tag_pvalue_dict[wildcards.tag_B]
    shell:
        "Rscript workflow/scripts/simgwas/merge_sim_sum_stats.R --sum_stats_file_A {input.sum_stats_file_A} --sum_stats_file_B {input.sum_stats_file_B} --file_A_stat_cols {params.file_A_stat_cols} --file_B_stat_cols {params.file_B_stat_cols} -o {output} -nt {threads}"

rule prune_merged_randomised_simulated_sum_stats:
    input:
        bim_file = "resources/1000g/euro/qc/chr1-22_qc.bim",
        pruned_range_file = "resources/plink_ranges/simgwas/pruned_ranges/window_{window}_step_{step}/all.prune.in",
        sum_stats_file = "results/simgwas/simulated_sum_stats/whole_genome_sum_stats/randomised/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/seed_{seed}_merged_sum_stats_tags_{tag_A}{tag_B}.tsv.gz"
    output:
        temp("results/simgwas/simulated_sum_stats/whole_genome_sum_stats/randomised/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/window_{window}_step_{step}/seed_{seed}_pruned_sum_stats_tags_{tag_A}{tag_B}.tsv.gz")
    threads: 4
    shell:
        "Rscript workflow/scripts/simgwas/prune_sim_sum_stats.R --sum_stats_file {input.sum_stats_file} --bim_file {input.bim_file} --prune_file {input.pruned_range_file} -o {output} -nt {threads}"
