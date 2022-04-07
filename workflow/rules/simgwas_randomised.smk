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

            if (chrom, block_no) not in shared_chrom_block_nos:
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

            if (chrom, block_no) not in shared_chrom_block_nos and (chrom, block_no) not in a_chrom_block_nos:
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
        # TODO get metadata indices
        # metadata (but we don't need all of these): [1,7], [7+(no_reps*4)+1, 7+(no_reps*4)+5]
    # beta: 7+(no_reps*2)+index_of_replicate
    # p-value: 7+(no_reps*3)+index_of_replicate
    # ncases: 7+(no_reps*4)+2
    # ncontrols: 7+(no_reps*4)+3
    # chr: 7+(no_reps*4)+4
    # block_effect_size: 7+(no_reps*4)+5
    # OLD
#        zsim_string = "\t".join([f"zsim.{x}" for x in range(1,21)])
#        p_string = "\t".join([f"p.{x}" for x in range(1,21)])
#        beta_string = "\t".join([f"betasim.{x}" for x in range(1,21)])
#        vbeta_string = "\t".join([f"vbetasim.{x}" for x in range(1,21)])
#        header_string = f"position\ta0\ta1\tid\tblock\tTYPE\tEUR\tzexp\t{zsim_string}\t{vbeta_string}\t{beta_string}\t{p_string}\tchosen_or\tncases\tncontrols\tchr\tblock_effect_size"

        with open(log.log, 'w') as logfile:
            p_column_name_A = f"p.{params.tag_index_A}"
            beta_column_name_A = f"betasim.{params.tag_index_A}"
            header_string_A = f"position\ta0\ta1\tid\tblock\tTYPE\tEUR\tzexp\t{beta_column_name_A}\t{p_column_name_A}\tncases\tncontrols\tchr"

            beta_column_index_A = 7+(params.no_reps*2)+params.tag_index_A
            p_column_index_A = 7+(params.no_reps*3)+params.tag_index_A

            ncases_column_index = 7+(params.no_reps*4)+2
            ncontrols_column_index = 7+(params.no_reps*4)+3
            chr_column_index = 7+(params.no_reps*4)+4

            awk_string_A = f"$1\t$2\t$3\t$4\t$5\t$6\t$7\t${beta_column_index_A}\t${p_column_index_A}\t${ncases_column_index}\t${ncontrols_column_index}\t${chr_column_index}"

            with open(params.uncomp_sum_stats_A, 'w') as file_A:
                file_A.write(f"{header_string_A}\n")

            for x in input.a_block_files:
                shell("zcat {x} | awk '{{ print {awk_string_A} }}' >> {params.uncomp_sum_stats_A}")
                logfile.write(f"{x}\n")

            shell(f"gzip {params.uncomp_sum_stats_A}")

            logfile.write("\n")

            p_column_name_B = f"p.{params.tag_index_B}"
            beta_column_name_B = f"betasim.{params.tag_index_B}"
            header_string_B = f"position\ta0\ta1\tid\tblock\tTYPE\tEUR\tzexp\t{p_column_name_B}\t{beta_column_name_B}\tncases\tncontrols\tchr"

            beta_column_index_B = 7+(params.no_reps*2)+params.tag_index_B
            p_column_index_B = 7+(params.no_reps*3)+params.tag_index_B

            awk_string_B = f"$1\t$2\t$3\t$4\t$5\t$6\t$7\t${beta_column_index_B}\t${p_column_index_B}\t${ncases_column_index}\t${ncontrols_column_index}\t${chr_column_index}"

            with open(params.uncomp_sum_stats_B, 'w') as file_B:
                file_B.write(f"{header_string_B}\n")

            for x in input.b_block_files:
                shell("zcat {x} | awk '{{ print {awk_string_B} }}' >> {params.uncomp_sum_stats_B}")
                logfile.write(f"{x}\n")

            shell(f"gzip {params.uncomp_sum_stats_B}")

# TODO Uh-oh: 8,998,661 in m1_seed_111_sum_stats_A.tsv.gz, 9,000,062 in seed_111_merged_sum_stats.tsv.gz
rule merge_randomised_simulated_sum_stats:
    input:
        sum_stats_file_A = "results/simgwas/simulated_sum_stats/whole_genome_sum_stats/randomised/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/{effect_blocks_A}_seed_{seed}_sum_stats_A.tsv.gz",
        sum_stats_file_B = "results/simgwas/simulated_sum_stats/whole_genome_sum_stats/randomised/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/{effect_blocks_B}_seed_{seed}_sum_stats_B.tsv.gz"
    output:
        temp("results/simgwas/simulated_sum_stats/whole_genome_sum_stats/randomised/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/seed_{seed}_merged_sum_stats.tsv.gz")
    params:
        no_reps = 20
    threads: 12
    shell:
        "Rscript workflow/scripts/simgwas/merge_sim_sum_stats.R --sum_stats_file_A {input.sum_stats_file_A} --sum_stats_file_B {input.sum_stats_file_B} --no_reps {params.no_reps} -o {output} -nt {threads}"

# TODO fix column order based on code in merge_sim_sum_stats.R and test column handling
rule prune_merged_randomised_simulated_sum_stats:
    input:
        bim_file = "resources/1000g/euro/qc/chr1-22_qc.bim",
        pruned_range_file = "resources/plink_ranges/simgwas/pruned_ranges/window_{window}_step_{step}/all.prune.in",
        sum_stats_file = "results/simgwas/simulated_sum_stats/whole_genome_sum_stats/randomised/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/seed_{seed}_merged_sum_stats.tsv.gz"
    output:
        temp("results/simgwas/simulated_sum_stats/whole_genome_sum_stats/randomised/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/window_{window}_step_{step}/seed_{seed}_pruned_sum_stats.tsv.gz")
    params:
        no_reps = 20
    threads: 4
    shell:
        "Rscript workflow/scripts/simgwas/prune_sim_sum_stats.R --sum_stats_file {input.sum_stats_file} --bim_file {input.bim_file} --prune_file {input.pruned_range_file} --no_reps {params.no_reps} -o {output} -nt {threads}"
