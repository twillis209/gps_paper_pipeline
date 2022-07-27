import re as re
import pandas as pd
import random

block_daf = pd.read_csv("resources/ldetect/available_blocks.tsv", sep = "\t")

block_dict = {}

for i in range(1,23):
    block_dict[i] = list(block_daf[block_daf["chr"] == i][block_daf["available"] == True]["block"])

effect_size_dict = {'s': 'small', 'm': 'medium', 'l': 'large', 'v': 'vlarge', 'h': 'huge', 'r': 'random', 'i': 'intermediate'}

def get_null_block_files(wildcards):
    null_block_file_format = "results/simgwas/simulated_sum_stats/block_sum_stats/{no_reps}_reps/null/{ncases}_{ncontrols}/chr%d/block_%s_seed_{seed}_sum_stats.tsv.gz"

    effect_block_files = get_effect_block_files(wildcards)

    null_block_files_to_omit = [re.sub('small|medium|large|vlarge|huge|intermediate|random_\d+-\d+', 'null', x) for x in effect_block_files]

    null_block_files = []

    for chrom, blocks in block_dict.items():
        null_block_files += [null_block_file_format % (chrom, y) for y in blocks]

    for x in null_block_files_to_omit:
        if x in null_block_files:
            null_block_files.remove(x)

    return null_block_files

def get_effect_block_files(wildcards):
    block_file_format = "results/simgwas/simulated_sum_stats/block_sum_stats/{no_reps}_reps/%s/{ncases}_{ncontrols}/chr%d/block_%s_seed_{seed}_sum_stats.tsv.gz"

    effect_block_files = []

    if wildcards.effect_blocks == 'null':
        return []

    for x in wildcards.effect_blocks.split('/')[-1].split('+'):
        block_match = re.match('^(\d+)-(.+)', x)

        chrom = int(block_match.group(1))

        if ':' in block_match.group(2):
            range_match = re.match('([smlvhi]|r\d+-\d+-)(\d+):(\d+)', block_match.group(2))

            if 'r' in range_match.group(1):
                effect = 'random_' + range_match.group(1)[1:-1]

                effect_block_files += [block_file_format % (chrom, effect, y) for y in range(int(range_match.group(2)), int(range_match.group(3))+1) if y in block_dict[chrom]]
            else:
                effect = effect_size_dict[range_match.group(1)]

                effect_block_files += [block_file_format % (chrom, effect, y) for y in range(int(range_match.group(2)), int(range_match.group(3))+1) if y in block_dict[chrom]]
        else:
            singleton_match = re.match('([smlvhi]|r\d+-\d+-)(\d+)', x)

            if 'r' in x:
                effect = 'random_' + singleton_match.group(1)[1:-1]

                if int(singleton_match.group(2)) in block_dict[chrom]:
                    effect_block_files += [block_file_format % (chrom, effect, int(singleton_match.group(2)))]
            else:
                effect = effect_size_dict[singleton_match.group(1)]

                if int(singleton_match.group(2)) in block_dict[chrom]:
                    effect_block_files.append(block_file_format % (chrom, effect, int(singleton_match.group(2))))

    return effect_block_files

# NB: Implements a messy anti-pattern of behaviour, shouldn't be called outside the context of the rule combine_randomised_effect_blocks
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

# NB: Implements a messy anti-pattern of behaviour, shouldn't be called outside the context of the rule combine_randomised_block_sum_stats_for_pair
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

            if (chrom, block_no) not in shared_chrom_block_nos:
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
            # if not in a block files or not in b block files, handle case where ncases/ncontrols are not the same
            if (chrom, block_no) not in a_chrom_block_nos and (chrom, block_no) not in shared_chrom_block_nos:
                block_files.append(f"results/simgwas/simulated_sum_stats/chr{chrom}/block_sum_stats/null/{wildcards.ncases_A}_{wildcards.ncontrols_A}/block_{block_no}_sum_stats.tsv.gz")
                a_block_files.append(f"results/simgwas/simulated_sum_stats/chr{chrom}/block_sum_stats/null/{wildcards.ncases_A}_{wildcards.ncontrols_A}/block_{block_no}_sum_stats.tsv.gz")
            if (chrom, block_no) not in b_chrom_block_nos and (chrom, block_no) not in shared_chrom_block_nos:
                b_block_files.append(f"results/simgwas/simulated_sum_stats/chr{chrom}/block_sum_stats/null/{wildcards.ncases_B}_{wildcards.ncontrols_B}/block_{block_no}_sum_stats.tsv.gz")
                if f"results/simgwas/simulated_sum_stats/chr{chrom}/block_sum_stats/null/{wildcards.ncases_B}_{wildcards.ncontrols_B}/block_{block_no}_sum_stats.tsv.gz" not in block_files:
                    block_files.append(f"results/simgwas/simulated_sum_stats/chr{chrom}/block_sum_stats/null/{wildcards.ncases_B}_{wildcards.ncontrols_B}/block_{block_no}_sum_stats.tsv.gz")

    return (block_files, a_block_files, b_block_files)

rule get_legend_files_with_euro_common_maf:
    input:
        "resources/simgwas/1000g/1000GP_Phase3_chr{ch}.legend.gz"
    output:
        "resources/simgwas/1000g/1000GP_Phase3_chr{ch}_eur_common_maf.legend.gz"
    params:
        "resources/simgwas/1000g/1000GP_Phase3_chr{ch}_eur_common_maf.legend"
    shell:
        """
        zcat {input} | awk -F' ' '($9 >= 0.01 && $9 <= 0.99) || NR == 1' >{params};
        gzip {params}
        """

rule get_euro_hap_files_with_metadata:
    input:
        legend_file = "resources/simgwas/1000g/1000GP_Phase3_chr{ch}.legend.gz",
        hap_file = "resources/simgwas/1000g/1000GP_Phase3_chr{ch}.hap.gz",
        samples_file = "resources/simgwas/1000g/chr{ch}.samples"
    output:
        output_hap_file = temp("resources/simgwas/1000g/1000GP_Phase3_chr{ch}_with_meta_eur_common_maf.hap.gz")
    params:
        uncomp_hap_file =  "resources/simgwas/1000g/1000GP_Phase3_chr{ch}_with_meta_eur_common_maf.hap",
        temp_hap_with_maf_file =  "resources/simgwas/1000g/1000GP_Phase3_chr{ch}_with_eur_maf.hap",
        temp_hap_filtered_maf_file =  "resources/simgwas/1000g/1000GP_Phase3_chr{ch}_filtered_eur_maf.hap",
        temp_hap_file =  "resources/simgwas/1000g/1000GP_Phase3_chr{ch}_with_meta.hap",
        temp_eur_cols_file =  "resources/simgwas/1000g/chr{ch}_EUR_cols.csv"
    shell:
        """
        paste -d' ' <(zcat {input.hap_file}) <(zcat {input.legend_file} | cut -d' ' -f9 | tail -n +2)  >{params.temp_hap_with_maf_file}

        awk -F' ' '$5009 >= 0.01 && $5009 <= 0.99' {params.temp_hap_with_maf_file} | cut -d' ' -f1-5008  >>{params.temp_hap_filtered_maf_file}

        for j in {{1..4}}; do
            for i in $(tail -n +2 {input.samples_file} | cut -f"$j" -d' ' | tr '\n' ' '); do echo -n "$i $i "; done >>{params.temp_hap_file}
        echo "" >>{params.temp_hap_file};
        done

        cat {params.temp_hap_filtered_maf_file} >>{params.temp_hap_file}

        head -n 3 {params.temp_hap_file} | tail -n 1 | awk -F' ' '{{for(i=1; i<=NF; i++){{ if($i == "EUR") {{printf "%s ", i}} }}}}'| sed 's/ /,/g' | sed 's/,$//g' > {params.temp_eur_cols_file}

        cut -d' ' -f$(cat {params.temp_eur_cols_file}) {params.temp_hap_file} >{params.uncomp_hap_file}; gzip {params.uncomp_hap_file};

        rm {params.temp_eur_cols_file} {params.temp_hap_file} {params.temp_hap_filtered_maf_file} {params.temp_hap_with_maf_file}
        """

rule write_out_available_blocks:
    output:
        "resources/simgwas/available_blocks.tsv"
    run:
        with open(output[0], 'w') as outfile:
            outfile.write("chr\tblock\n")
            for x in sorted(block_dict.items()):
                for y in x[1]:
                    outfile.write(f"{x[0]}\t{y}\n")

# Generates rules to generate LD block files for whole genome
# After https://stackoverflow.com/a/49004234
for chrom, blocks in block_dict.items():
    rule:
        input:
            haplotype_file = "resources/simgwas/1000g/1000GP_Phase3_chr{ch}_with_meta_eur_common_maf.hap.gz".format(ch = chrom),
            legend_file = "resources/simgwas/1000g/1000GP_Phase3_chr{ch}_eur_common_maf.legend.gz".format(ch = chrom),
            block_file = "resources/ldetect/blocks.txt"
        output:
            ["resources/simgwas/1000g/blockwise/chr{ch}/block_{block_no}.hap.gz".format(ch = chrom, block_no = x) for x in blocks],
            ["resources/simgwas/1000g/blockwise/chr{ch}/block_{block_no}.legend.gz".format(ch = chrom, block_no = x) for x in blocks]
        params:
            output_root = "resources/simgwas/1000g/blockwise/chr{ch}".format(ch = chrom),
            chr_no = chrom
        threads: 4
        shell:
            "Rscript workflow/scripts/simgwas/write_out_ld_block_files.R --hap_file {input.haplotype_file} --leg_file {input.legend_file} -b {input.block_file} --chr_no {params.chr_no} -o {params.output_root} -nt {threads}"

rule compute_block_ld_matrix:
    input:
        block_haplotype_file = ancient("resources/simgwas/1000g/blockwise/chr{ch}/block_{block}.hap.gz"),
        block_legend_file = ancient("resources/simgwas/1000g/blockwise/chr{ch}/block_{block}.legend.gz"),
        block_file = "resources/ldetect/blocks.txt"
    output:
        "results/simgwas/chr{ch}_ld_matrices/block_{block}_ld_matrix.RData"
    threads: 4
    resources:
        runtime = 60,
        mem_mb=get_mem_mb
    shell:
        "Rscript workflow/scripts/simgwas/compute_block_ld_matrix.R --hap_file {input.block_haplotype_file} --leg_file {input.block_legend_file} --output_file {output} -nt {threads}"

rule simulate_sum_stats_by_ld_block:
    input:
        bim_file = ancient("resources/1000g/chr{ch}.bim"),
        block_haplotype_file = ancient("resources/simgwas/1000g/blockwise/chr{ch}/block_{block}.hap.gz"),
        block_legend_file = ancient("resources/simgwas/1000g/blockwise/chr{ch}/block_{block}.legend.gz"),
        ld_mat_file = ancient("results/simgwas/chr{ch}_ld_matrices/block_{block}_ld_matrix.RData")
    output:
        "results/simgwas/simulated_sum_stats/block_sum_stats/{no_reps}_reps/{effect_size}/{ncases,\d+}_{ncontrols,\d+}/chr{ch}/block_{block,\d+}_seed_{seed,\d+}_sum_stats.tsv.gz"
    threads: 1
    resources:
        mem_mb = get_mem_mb,
        runtime = 30
    group: 'simulate'
    benchmark: 'results/benchmarks/simulate_sum_stats_by_ld_block/{no_reps}_reps/{effect_size}/{ncases}_{ncontrols}/chr{ch}/block_{block}_seed_{seed}_sum_stats.txt'
    shell:
        "Rscript workflow/scripts/simgwas/simulate_sum_stats_by_ld_block.R --hap_file {input.block_haplotype_file} --leg_file {input.block_legend_file} --bim_file {input.bim_file} --ld_mat_file {input.ld_mat_file} --chr_no {wildcards.ch} --causal_variant_ind 2000 --effect_size {wildcards.effect_size} --no_controls {wildcards.ncontrols} --no_cases {wildcards.ncases} --no_reps {wildcards.no_reps} --seed {wildcards.seed} -o {output} -nt {threads}"

rule get_causal_variant_by_ld_block:
    input:
        bim_file = ancient("resources/1000g/chr{ch}.bim"),
        block_haplotype_file = ancient("resources/simgwas/1000g/blockwise/chr{ch}/block_{block}.hap.gz"),
        block_legend_file = ancient("resources/simgwas/1000g/blockwise/chr{ch}/block_{block}.legend.gz"),
        ld_mat_file = ancient("results/simgwas/chr{ch}_ld_matrices/block_{block}_ld_matrix.RData")
    output:
        temp("results/simgwas/simulated_sum_stats/chr{ch}/block_sum_stats/null/10000_10000/block_{block}_causal_variant.tsv")
    threads: 4
    resources:
        mem_mb=get_mem_mb,
        runtime = 1
    shell:
        "Rscript workflow/scripts/simgwas/get_causal_variant_by_ld_block.R --hap_file {input.block_haplotype_file} --leg_file {input.block_legend_file} --bim_file {input.bim_file} --ld_mat_file {input.ld_mat_file} --chr_no {wildcards.ch} --causal_variant_ind 2000 -o {output} -nt {threads}"

rule combine_block_sum_stats:
    input:
        null_block_files = ancient(get_null_block_files),
        effect_block_files = ancient(get_effect_block_files)
    output:
        temp("results/simgwas/simulated_sum_stats/whole_genome_sum_stats/{ncases,\d+}_{ncontrols,\d+}/{effect_blocks}_sum_stats.tsv.gz")
    run:
        for x in input:
            shell(f"cat {x} >> {output}")

rule make_whole_genome_metadata_file:
    input:
        [["results/simgwas/simulated_sum_stats/chr%d/block_sum_stats/null/10000_10000/block_%d_sum_stats.tsv.gz" % (chrom, block) for block in block_dict[chrom]] for chrom in range(1,23)]
    output:
        temp("results/simgwas/simulated_sum_stats/whole_genome_sum_stats/metadata_only.tsv.gz")
    params:
        uncomp_output = "results/simgwas/simulated_sum_stats/whole_genome_sum_stats/metadata_only.tsv"
    run:
        header_string = "position\ta0\ta1\tid\tTYPE\tEUR\tchr"

        shell(f"echo -e \"{header_string}\" > {params.uncomp_output}")

        for x in input:
            print(x)
            shell(f"zcat {x} | cut -f1-4,6-7,92 >> {params.uncomp_output}")

        shell("gzip {params.uncomp_output}")

rule make_simgwas_plink_ranges:
    input:
        sum_stats_metadata_file = "results/simgwas/simulated_sum_stats/whole_genome_sum_stats/metadata_only.tsv.gz",
        bim_file = "resources/1000g/euro/qc/chr1-22_qc.bim",
    output:
        [("resources/plink_ranges/1000g/chr%d.txt" % x for x in range(1,23))]
    threads: 4
    shell:
        "Rscript workflow/scripts/simgwas/make_simgwas_plink_ranges.R --sum_stats_file {input.sum_stats_metadata_file} --input_bim_file {input.bim_file} --output_range_files {output} -nt {threads} --bp_pos 1 --chr_pos 92 --a0_pos 3 --a1_pos 4"

rule merge_simulated_sum_stats:
    input:
        sum_stats_file_A = "results/simgwas/simulated_sum_stats/whole_genome_sum_stats/{ncases_A}_{ncontrols_A}/{effect_blocks_A}_sum_stats.tsv.gz",
        sum_stats_file_B = "results/simgwas/simulated_sum_stats/whole_genome_sum_stats/{ncases_B}_{ncontrols_B}/{effect_blocks_B}_sum_stats.tsv.gz"
    output:
        temp("results/simgwas/simulated_sum_stats/merged/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_sum_stats.tsv.gz")
    params:
        no_reps = 20
    threads: 4
    shell:
        "Rscript workflow/scripts/simgwas/merge_sim_sum_stats.R --sum_stats_file_A {input.sum_stats_file_A} --sum_stats_file_B {input.sum_stats_file_B} --no_reps {params.no_reps} -o {output} -nt {threads}"

rule prune_merged_sim_sum_stats:
    input:
        bim_file = "resources/1000g/euro/qc/chr1-22_qc.bim",
        pruned_range_file = "resources/plink_ranges/1000g/pruned_ranges/window_{window}_step_{step}/all.prune.in",
        sum_stats_file = "results/simgwas/simulated_sum_stats/merged/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_sum_stats.tsv.gz"
    output:
        temp("results/simgwas/simulated_sum_stats/pruned/window_{window}_step_{step}/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_sum_stats.tsv")
    params:
        no_reps = 20
    threads: 4
    shell:
        "Rscript workflow/scripts/simgwas/prune_sim_sum_stats.R --sum_stats_file {input.sum_stats_file} --bim_file {input.bim_file} --prune_file {input.pruned_range_file} -o {output} -nt {threads} --no_reps {params.no_reps}"

# TODO chr6:block86 CV duplication problem
rule get_causal_variants:
    input:
        [["results/simgwas/simulated_sum_stats/chr%d/block_sum_stats/null/10000_10000/block_%d_causal_variant.tsv" % (x,y) for y in block_dict[x]] for x in range(1,23)]
    output:
        "results/simgwas/combined_causal_variants.tsv"
    run:
        for i,x in enumerate(input):
            if i == 0:
                shell("cat %s > %s" % (x, output[0]))
            else:
                shell("cat %s | tail -n +2  >> %s" % (x, output[0]))
