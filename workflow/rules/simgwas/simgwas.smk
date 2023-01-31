import re as re
import pandas as pd
import random

localrules: get_causal_variants

effect_size_dict = {"s": "small", "m": "medium", "l": "large", "v": "vlarge", "h": "huge", "r": "random", "i": "infinitesimal", "t": "tiny"}
cv_per_block_dict = {"small": 1, "medium": 1, "large": 1, "vlarge": 1, "huge": 1, "tiny": 2, "null": 0}

cv_index_dict = {1: [2000], 2: [1000, 2000]}

odds_ratio_dict = {"small": 1.05, "medium": 1.2, "null": 1, "tiny": 1.01}

include: "simgwas_functions.py"

rule get_legend_files_with_euro_common_maf:
    input:
        "resources/simgwas/1000g/1000GP_Phase3_{chr}.legend.gz"
    output:
        "resources/simgwas/1000g/1000GP_Phase3_{chr}_eur_common_maf.legend.gz"
    params:
        "resources/simgwas/1000g/1000GP_Phase3_{chr}_eur_common_maf.legend"
    shell:
        """
        zcat {input} | awk -F' ' '($9 >= 0.01 && $9 <= 0.99) || NR == 1' >{params};
        gzip {params}
        """

rule get_euro_hap_files_with_metadata:
    input:
        legend_file = "resources/simgwas/1000g/1000GP_Phase3_{chr}.legend.gz",
        hap_file = "resources/simgwas/1000g/1000GP_Phase3_{chr}.hap.gz",
        samples_file = "resources/simgwas/1000g/{chr}.samples"
    output:
        output_hap_file = temp("resources/simgwas/1000g/1000GP_Phase3_{chr}_with_meta_eur_common_maf.hap.gz")
    params:
        uncomp_hap_file =  "resources/simgwas/1000g/1000GP_Phase3_{chr}_with_meta_eur_common_maf.hap",
        temp_hap_with_maf_file =  "resources/simgwas/1000g/1000GP_Phase3_{chr}_with_eur_maf.hap",
        temp_hap_filtered_maf_file =  "resources/simgwas/1000g/1000GP_Phase3_{chr}_filtered_eur_maf.hap",
        temp_hap_file =  "resources/simgwas/1000g/1000GP_Phase3_{chr}_with_meta.hap",
        temp_eur_cols_file =  "resources/simgwas/1000g/{chr}_EUR_cols.csv"
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
for chrom in range(1,23):
    blocks = block_daf.query('chr == @chrom')['block']

    rule:
        input:
            haplotype_file = f"resources/simgwas/1000g/1000GP_Phase3_chr{chrom}_with_meta_eur_common_maf.hap.gz",
            legend_file = f"resources/simgwas/1000g/1000GP_Phase3_chr{chrom}_eur_common_maf.legend.gz",
            block_file = "resources/ldetect/blocks.txt"
        output:
            [f"resources/simgwas/1000g/blockwise/chr{chrom}/block_{block}.hap.gz" for block in blocks],
            [f"resources/simgwas/1000g/blockwise/chr{chrom}/block_{block}.legend.gz" for block in blocks]
        params:
            output_root = f"resources/simgwas/1000g/blockwise/chr{chrom}",
            chr_no = chrom
        threads: 4
        shell:
            "Rscript workflow/scripts/simgwas/write_out_ld_block_files.R --hap_file {input.haplotype_file} --leg_file {input.legend_file} -b {input.block_file} --chr_no {params.chr_no} -o {params.output_root} -nt {threads}"

rule compute_block_ld_matrix:
    input:
        block_haplotype_file = ancient("resources/simgwas/1000g/blockwise/{chr}/block_{block}.hap.gz"),
        block_legend_file = ancient("resources/simgwas/1000g/blockwise/{chr}/block_{block}.legend.gz"),
        block_file = "resources/ldetect/blocks.txt"
    output:
        "results/simgwas/{chr}_ld_matrices/block_{block}_ld_matrix.RData"
    threads: 4
    resources:
        runtime = 60,
        mem_mb=get_mem_mb
    shell:
        "Rscript workflow/scripts/simgwas/compute_block_ld_matrix.R --hap_file {input.block_haplotype_file} --leg_file {input.block_legend_file} --output_file {output} -nt {threads}"

rule simulate_sum_stats_by_ld_block:
    input:
        bim_file = ancient("resources/1000g/{chr}.bim"),
        block_haplotype_file = ancient("resources/simgwas/1000g/blockwise/{chr}/block_{block}.hap.gz"),
        block_legend_file = ancient("resources/simgwas/1000g/blockwise/{chr}/block_{block}.legend.gz"),
        ld_mat_file = ancient("results/simgwas/{chr}_ld_matrices/block_{block}_ld_matrix.RData")
    output:
        "results/simgwas/simulated_sum_stats/block_sum_stats/{no_reps}_reps/{effect_size}/{no_cvariants}_cv/{ncases,\d+}_{ncontrols,\d+}/{chr}/block_{block,\d+}_seed_{seed,\d+}_sum_stats.tsv.gz"
    params:
        chr_no = lambda wildcards: wildcards.chr.replace('chr', ''),
        causal_variant_indices = lambda wildcards: cv_index_dict[int(wildcards.no_cvariants)],
        odds_ratio = lambda wildcards: odds_ratio_dict[wildcards.effect_size],
        no_cases = lambda wildcards: int(wildcards.ncases),
        no_controls = lambda wildcards: int(wildcards.ncontrols),
        no_reps = lambda wildcards: int(wildcards.no_reps)
    threads: 3
    resources:
        mem_mb = get_mem_mb,
        runtime = get_simulation_runtime
    group: 'simulate'
    script:
        "../../scripts/simgwas/simulate_sum_stats_by_ld_block.R"

rule get_causal_variants_by_ld_block:
    input:
        bim_file = ancient("resources/1000g/{chr}.bim"),
        block_haplotype_file = ancient("resources/simgwas/1000g/blockwise/{chr}/block_{block}.hap.gz"),
        block_legend_file = ancient("resources/simgwas/1000g/blockwise/{chr}/block_{block}.legend.gz"),
        ld_mat_file = ancient("results/simgwas/{chr}_ld_matrices/block_{block}_ld_matrix.RData")
    params:
        chr_no = lambda wildcards: wildcards.chr.replace('chr', ''),
        causal_variant_indices = lambda wildcards: cv_index_dict[2]
    output:
        temp("results/simgwas/simulated_sum_stats/block_sum_stats/{chr}/null/10000_10000/block_{block}_causal_variants.tsv")
    threads: 4
    resources:
        mem_mb = get_mem_mb,
        runtime = 5
    group: 'causal_variants'
    script:
        "../../scripts/simgwas/get_causal_variants_by_ld_block.R"

rule get_causal_variants:
    input:
        [[f"results/simgwas/simulated_sum_stats/block_sum_stats/chr{chrom}/null/10000_10000/block_{block}_causal_variants.tsv" for block in block_daf.query('chr == @chrom')['block']] for chrom in range(1,23)]
    output:
        "results/simgwas/combined_causal_variants.tsv"
    run:
        for i,x in enumerate(input):
            if i == 0:
                shell("cat %s > %s" % (x, output[0]))
            else:
                shell("cat %s | tail -n +2  >> %s" % (x, output[0]))

rule make_whole_genome_metadata_file:
    input:
        [["results/simgwas/simulated_sum_stats/chr%d/block_sum_stats/null/10000_10000/block_%d_sum_stats.tsv.gz" % (chrom, block) for block in block_daf.query('chr == @chrom')['block']] for chrom in range(1,23)]
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
