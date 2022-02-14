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

#def get_block_hap_files(wildcards):
#    return ["resources/simgwas/1000g/blockwise/chr{ch}/block_%s.hap.gz" % i for i in block_dict[wildcards.ch]]
#
#def get_block_leg_files(wildcards):
#    return ["resources/simgwas/1000g/blockwise/chr{ch}/block_%s.legend.gz" % i for i in block_dict[wildcards.ch]]

def get_null_block_files(wildcards):
    return ["results/simgwas/simulated_sum_stats/chr{ch}/block_sum_stats/null/{ncases}_{ncontrols}/block_%d_sum_stats.tsv.gz" % i for i in block_dict[int(wildcards.ch)] if i not in [int(x[1:]) for x in wildcards.effect_block.split(':')]]

def get_effect_block_files(wildcards):
    effect_size_dict = {'s': 'small', 'm': 'medium', 'l': 'large', 'v': 'vlarge', 'h': 'huge'}

    token_re = re.compile('([smlv])(\d+)')

    tokens = wildcards.effect_block.split(':')

    matches = [token_re.match(x) for x in tokens]

    return ["results/simgwas/simulated_sum_stats/chr{ch}/block_sum_stats/%s/{ncases}_{ncontrols}/block_%s_sum_stats.tsv.gz" % (effect_size_dict[x.group(1)], x.group(2)) for x in matches]

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
            "Rscript workflow/scripts/write_out_ld_block_files.R --hap_file {input.haplotype_file} --leg_file {input.legend_file} -b {input.block_file} --chr_no {params.chr_no} -o {params.output_root} -nt {threads}"

# TODO remove me
rule produce_all_blocks:
    input:
        ["resources/simgwas/1000g/blockwise/chr%d/block_0.hap.gz" % i for i in range(1,23)]

rule compute_block_ld_matrix:
    input:
        block_haplotype_file = ancient("resources/simgwas/1000g/blockwise/chr{ch}/block_{block}.hap.gz"),
        block_legend_file = ancient("resources/simgwas/1000g/blockwise/chr{ch}/block_{block}.legend.gz"),
        block_file = "resources/ldetect/blocks.txt"
    output:
        "results/simgwas/chr{ch}_ld_matrices/block_{block}_ld_matrix.RData"
    threads: 4
    resources:
        time = "4:00:00",
        mem_mb=get_mem_mb
    shell:
        "Rscript workflow/scripts/compute_block_ld_matrix.R --hap_file {input.block_haplotype_file} --leg_file {input.block_legend_file} --output_file {output} -nt {threads}"

rule simulate_sum_stats_by_ld_block:
    input:
        bim_file = ancient("resources/1000g/chr{ch}.bim"),
        block_haplotype_file = ancient("resources/simgwas/1000g/blockwise/chr{ch}/block_{block}.hap.gz"),
        block_legend_file = ancient("resources/simgwas/1000g/blockwise/chr{ch}/block_{block}.legend.gz"),
        ld_mat_file = ancient("results/simgwas/chr{ch}_ld_matrices/block_{block}_ld_matrix.RData")
    output:
        "results/simgwas/simulated_sum_stats/chr{ch}/block_sum_stats/{effect_size}/{ncases}_{ncontrols}/block_{block}_sum_stats.tsv.gz"
    threads: 2
    resources:
        mem_mb=get_mem_mb,
        time = "00:15:00"
    shell:
        "Rscript workflow/scripts/simulate_sum_stats_by_ld_block.R --hap_file {input.block_haplotype_file} --leg_file {input.block_legend_file} --bim_file {input.bim_file} --ld_mat_file {input.ld_mat_file} --chr_no {wildcards.ch} --causal_variant_ind 2000 --effect_size {wildcards.effect_size} --no_controls {wildcards.ncontrols} --no_cases {wildcards.ncases} --no_reps 2 -o {output} -nt {threads}"

"""
rule simulate_all_null_blocks:
    input:
        [["results/simgwas/simulated_sum_stats/chr{ch}/block_sum_stats/null/10000_10000/block_{block}_sum_stats.tsv.gz".format(ch = x, block = y) for y in block_dict[x]] for x in range(1,23)]
"""

rule combine_block_sum_stats:
    input:
        null_block_files = ancient(get_null_block_files),
        effect_block_files = ancient(get_effect_block_files)
    output:
        temp("results/simgwas/simulated_sum_stats/chr{ch}/whole_chr_sum_stats/{ncases}_{ncontrols}/blocks_{first}_{last}_{effect_block,[smlv]\d+(:[smlv]\d+)*}_sum_stats.tsv.gz")
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
        sum_stats_file_A = "results/simgwas/simulated_sum_stats/chr{ch_A}/whole_chr_sum_stats/{ncases_A}_{ncontrols_A}/blocks_{first_A}_{last_A}_{effect_block_A}_sum_stats.tsv.gz",
        sum_stats_file_B = "results/simgwas/simulated_sum_stats/chr{ch_B}/whole_chr_sum_stats/{ncases_B}_{ncontrols_B}/blocks_{first_B}_{last_B}_{effect_block_B}_sum_stats.tsv.gz"
    output:
        temp("results/simgwas/simulated_sum_stats/merged/chr{ch_A}_{ch_B}/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/blocks_{first_A}_{last_A}_{effect_block_A}_{first_B}_{last_B}_{effect_block_B}_sum_stats.tsv")
    threads: 4
    shell:
        "Rscript workflow/scripts/merge_sim_sum_stats.R --sum_stats_file_A {input.sum_stats_file_A} --sum_stats_file_B {input.sum_stats_file_B} --no_reps 2 -o {output} -nt {threads}"

# TODO generate new prune files specific to simulated sum stats
rule prune_merged_sim_sum_stats:
    input:
      bim_file = ancient("resources/1000g/euro/qc/chr{ch_A}_qc.bim"),
      pruned_range_file = ancient("resources/plink_ranges/ukbb/pruned_ranges/window_{window}_step_{step}/all.prune.in"),
      sum_stats_file = "results/simgwas/simulated_sum_stats/merged/chr{ch_A}_{ch_B}/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/blocks_{first_A}_{last_A}_{effect_block_A}_{first_B}_{last_B}_{effect_block_B}_sum_stats.tsv"
    output:
        temp("results/simgwas/simulated_sum_stats/pruned/window_{window}_step_{step}/chr{ch_A}_{ch_B}/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/blocks_{first_A}_{last_A}_{effect_block_A}_{first_B}_{last_B}_{effect_block_B}_sum_stats.tsv")
    threads: 4
    shell:
      "Rscript workflow/scripts/prune_sim_sum_stats.R --sum_stats_file {input.sum_stats_file} --bim_file {input.bim_file} --prune_file {input.pruned_range_file} -o {output} -nt {threads}"
