import re as re
import pandas as pd

# NB: Doesn't work when I specify the chromosome no. wildcard as '{chr}'
chr1_ld_matrices = ["results/simgwas/chr1_ld_matrices/chr1_block_%d_ld_matrix.RData" % i for i in range(122)]
chr21_ld_matrices = ["results/simgwas/chr21_ld_matrices/chr21_block_%d_ld_matrix.RData" % i for i in range(23)]

block_daf = pd.read_csv("resources/ldetect/blocks.txt", sep = " ", names = ["chr_block", "chr", "start", "stop"])

block_dict = {}

for i in range(1,23):
    block_dict[i] = list(block_daf[block_daf["chr"] == i]["chr_block"].apply(lambda x: int(x.split('_block')[1])))

# No SNPs in chr5:block108
block_dict[5] = range(108)
# No SNPs?
block_dict[17].remove(42)
block_dict[18].remove(45)
block_dict[18].remove(46)
block_dict[20].remove(35)
# Can't compute this one inside four hours; maybe it's the MHC?
block_dict[6].remove(25)

effect_size_dict = {'s': 'small', 'm': 'medium', 'l': 'large', 'v': 'vlarge', 'h': 'huge', 'r': 'random'}

def get_whole_genome_null_block_files(wildcards):
    null_block_file_format = "results/simgwas/simulated_sum_stats/chr%d/block_sum_stats/null/{ncases}_{ncontrols}/block_%s_sum_stats.tsv.gz"

    effect_block_files = get_whole_genome_effect_block_files(wildcards)

    null_block_files_to_omit = [re.sub('small|medium|large|vlarge|huge|random_\d+-\d+', 'null', x) for x in effect_block_files]

    null_block_files = []

    for chrom, blocks in block_dict.items():
        null_block_files += [null_block_file_format % (chrom, y) for y in blocks]

    for x in null_block_files_to_omit:
        if x in null_block_files:
            null_block_files.remove(x)

    return null_block_files

def get_whole_genome_effect_block_files(wildcards):
    block_file_format = "results/simgwas/simulated_sum_stats/chr%d/block_sum_stats/%s/{ncases}_{ncontrols}/block_%s_sum_stats.tsv.gz"

    effect_block_files = []

    for x in wildcards.effect_blocks.split('/')[-1].split('+'):
        block_match = re.match('^(\d+)-(.+)', x)

        chrom = int(block_match.group(1))

        if ':' in block_match.group(2):
            range_match = re.match('([smlvh]|r\d+-\d+-)(\d+):(\d+)', block_match.group(2))

            if 'r' in range_match.group(1):
                effect = 'random_' + range_match.group(1)[1:-1]

                effect_block_files += [block_file_format % (chrom, effect, y) for y in range(int(range_match.group(2)), int(range_match.group(3))+1)]
            else:
                effect = effect_size_dict[range_match.group(1)]

                effect_block_files += [block_file_format % (chrom, effect, y) for y in range(int(range_match.group(2)), int(range_match.group(3))+1)]
        else:
            singleton_match = re.match('([smlvh]|r\d+-\d+-)(\d+)', x)

            if 'r' in x:
                effect = 'random_' + singleton_match.group(1)[1:-1]

                effect_block_files += [block_file_format % (chrom, effect, int(singleton_match.group(2)))]

            else:
                effect = effect_size_dict[singleton_match.group(1)]

                effect_block_files.append(block_file_format % (chrom, effect, int(singleton_match.group(2))))

    return effect_block_files

# TODO add check for overlap in effect blocks, e.g. r3-123-1\+r1-0-2
def get_effect_block_files(wildcards):
    block_file_format = "results/simgwas/simulated_sum_stats/chr{ch}/block_sum_stats/%s/{ncases}_{ncontrols}/block_%s_sum_stats.tsv.gz"

    effect_block_files = []

    for x in wildcards.effect_blocks.split('/')[-1].split('+'):
        if ':' in x:
            range_match = re.match('([smlvh]|r\d+-\d+-)(\d+):(\d+)', x)

            if 'r' in x:
                effect = 'random_' + range_match.group(1)[1:-1]

                effect_block_files += [block_file_format % (effect, y) for y in range(int(range_match.group(2)), int(range_match.group(3))+1)]

            else:
                effect = effect_size_dict[range_match.group(1)]
                effect_block_files += [block_file_format % (effect, y) for y in range(int(range_match.group(2)), int(range_match.group(3))+1)]
        else:
            singleton_match = re.match('([smlvh]|r\d+-\d+-)(\d+)', x)

            if 'r' in x:
                effect = 'random_' + singleton_match.group(1)[1:-1]

                effect_block_files += [block_file_format % (effect, int(singleton_match.group(2)))]

            else:
                effect = effect_size_dict[singleton_match.group(1)]

                effect_block_files.append(block_file_format % (effect, int(singleton_match.group(2))))

    return effect_block_files

def get_null_block_files(wildcards):
    null_block_file_format = "results/simgwas/simulated_sum_stats/chr{ch}/block_sum_stats/null/{ncases}_{ncontrols}/block_%d_sum_stats.tsv.gz"

    effect_block_numbers = [int(re.search('block_(\d+)_sum_stats.tsv.gz', x).group(1)) for x in get_effect_block_files(wildcards)]

    return [null_block_file_format % i for i in block_dict[int(wildcards.ch)] if i not in effect_block_numbers]

# TODO make the legend and hap whole files temp
rule get_legend_files_with_euro_common_maf:
    input:
        "resources/simgwas/1000g/1000GP_Phase3_chr{ch}.legend.gz"
    output:
        temp("resources/simgwas/1000g/1000GP_Phase3_chr{ch}_eur_common_maf.legend.gz")
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
        time = "8:00:00",
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
        "results/simgwas/simulated_sum_stats/chr{ch}/block_sum_stats/{effect_size}/{ncases}_{ncontrols}/block_{block}_sum_stats.tsv.gz"
    threads: 5
    resources:
        mem_mb=get_mem_mb,
        time = "00:15:00"
    shell:
        "Rscript workflow/scripts/simgwas/simulate_sum_stats_by_ld_block.R --hap_file {input.block_haplotype_file} --leg_file {input.block_legend_file} --bim_file {input.bim_file} --ld_mat_file {input.ld_mat_file} --chr_no {wildcards.ch} --causal_variant_ind 2000 --effect_size {wildcards.effect_size} --no_controls {wildcards.ncontrols} --no_cases {wildcards.ncases} --no_reps 5 -o {output} -nt {threads}"

rule combine_block_sum_stats_for_single_chr:
    input:
        null_block_files = get_null_block_files,
        effect_block_files = get_effect_block_files
    output:
        temp("results/simgwas/simulated_sum_stats/chr{ch}/whole_chr_sum_stats/{ncases}_{ncontrols}/{effect_blocks}_sum_stats.tsv.gz")
    run:
        for i,x in enumerate(input):
            if i == 0:
                shell("zcat %s >> %s" % (x, output[0].replace('.gz', '')))
            else:
                shell("zcat %s | tail -n +2  >> %s" % (x, output[0].replace('.gz', '')))
        shell("gzip %s" % output[0].replace('.gz', ''))

rule combine_block_sum_stats_whole_genome:
    input:
        null_block_files = ancient(get_whole_genome_null_block_files),
        effect_block_files = ancient(get_whole_genome_effect_block_files)
    output:
        temp("results/simgwas/simulated_sum_stats/whole_genome_sum_stats/{ncases}_{ncontrols}/{effect_blocks}_sum_stats.tsv.gz")
    run:
        for i,x in enumerate(input):
            if i == 0:
                shell("zcat %s >> %s" % (x, output[0].replace('.gz', '')))
            else:
                shell("zcat %s | tail -n +2  >> %s" % (x, output[0].replace('.gz', '')))
        shell("gzip %s" % output[0].replace('.gz', ''))

# TODO finish me then change the prune rule
rule merge_simulated_sum_stats:
    input:
        sum_stats_file_A = "results/simgwas/simulated_sum_stats/chr{ch_A}/whole_chr_sum_stats/{ncases_A}_{ncontrols_A}/{effect_blocks_A}_sum_stats.tsv.gz",
        sum_stats_file_B = "results/simgwas/simulated_sum_stats/chr{ch_B}/whole_chr_sum_stats/{ncases_B}_{ncontrols_B}/{effect_blocks_B}_sum_stats.tsv.gz"
    output:
        temp("results/simgwas/simulated_sum_stats/merged/chr{ch_A}/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_sum_stats.tsv")
    threads: 4
    shell:
        "Rscript workflow/scripts/simgwas/merge_sim_sum_stats.R --sum_stats_file_A {input.sum_stats_file_A} --sum_stats_file_B {input.sum_stats_file_B} --no_reps 5 -o {output} -nt {threads}"

# TODO generate new prune files specific to simulated sum stats
rule prune_merged_sim_sum_stats:
    input:
      bim_file = ancient("resources/1000g/euro/qc/chr{ch_A}_qc.bim"),
      pruned_range_file = ancient("resources/plink_ranges/ukbb/pruned_ranges/window_{window}_step_{step}/all.prune.in"),
      sum_stats_file = "results/simgwas/simulated_sum_stats/merged/chr{ch_A}/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_sum_stats.tsv"
    output:
        temp("results/simgwas/simulated_sum_stats/pruned/window_{window}_step_{step}/chr{ch_A}/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_sum_stats.tsv")
    threads: 4
    shell:
      "Rscript workflow/scripts/simgwas/prune_sim_sum_stats.R --sum_stats_file {input.sum_stats_file} --bim_file {input.bim_file} --prune_file {input.pruned_range_file} -o {output} -nt {threads}"
