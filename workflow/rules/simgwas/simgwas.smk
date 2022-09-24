import re as re
import pandas as pd
import random

effect_size_dict = {'s': 'small', 'm': 'medium', 'l': 'large', 'v': 'vlarge', 'h': 'huge', 'r': 'random', 'i': 'intermediate'}

include: "simgwas_functions.py"

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
        temp("results/simgwas/simulated_sum_stats/block_sum_stats/{no_reps}_reps/{effect_size}/{ncases,\d+}_{ncontrols,\d+}/chr{ch}/block_{block,\d+}_seed_{seed,\d+}_sum_stats.tsv.gz")
    threads: 3
    resources:
        mem_mb = get_mem_mb,
        runtime = get_simulation_runtime
    group: 'simulate'
    benchmark: 'results/benchmarks/simulate_sum_stats_by_ld_block/{no_reps}_reps/{effect_size}/{ncases}_{ncontrols}/chr{ch}/block_{block}_seed_{seed}_sum_stats.txt'
    shell:
        "Rscript workflow/scripts/simgwas/simulate_sum_stats_by_ld_block.R --hap_file {input.block_haplotype_file} --leg_file {input.block_legend_file} --bim_file {input.bim_file} --ld_mat_file {input.ld_mat_file} --chr_no {wildcards.ch} --causal_variant_ind 2000 --effect_size {wildcards.effect_size} --no_controls {wildcards.ncontrols} --no_cases {wildcards.ncases} --no_reps {wildcards.no_reps} --seed {wildcards.seed} -o {output} -nt {threads}"

rule unzip_block_file:
    input:
        ancient("results/simgwas/simulated_sum_stats/block_sum_stats/{no_reps}_reps/{effect_size}/{ncases}_{ncontrols}/chr{ch}/block_{block}_seed_{seed}_sum_stats.tsv.gz")
    output:
        "results/simgwas/simulated_sum_stats/block_sum_stats/{no_reps}_reps/{effect_size}/{ncases,\d+}_{ncontrols,\d+}/chr{ch}/block_{block,\d+}_seed_{seed,\d+}_sum_stats.tsv"
    threads: 1
    resources:
        runtime = 1
    group: "simulate"
    shell:
        """
        zcat {input} >{output}
        """

rule get_causal_variant_by_ld_block:
    input:
        bim_file = ancient("resources/1000g/chr{ch}.bim"),
        block_haplotype_file = ancient("resources/simgwas/1000g/blockwise/chr{ch}/block_{block}.hap.gz"),
        block_legend_file = ancient("resources/simgwas/1000g/blockwise/chr{ch}/block_{block}.legend.gz"),
        ld_mat_file = ancient("results/simgwas/chr{ch}_ld_matrices/block_{block}_ld_matrix.RData")
    output:
        temp("results/simgwas/simulated_sum_stats/block_sum_stats/chr{ch}/null/10000_10000/block_{block}_causal_variant.tsv")
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
        temp("results/simgwas/simulated_sum_stats/whole_genome_sum_stats/{no_reps}_reps/{ncases,\d+}_{ncontrols,\d+}/{effect_blocks}_sum_stats.tsv.gz")
    run:
        for x in input:
            shell(f"cat {x} >> {output}")

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
        [[f"results/simgwas/simulated_sum_stats/chr{chrom}/block_sum_stats/null/10000_10000/block_{block}_causal_variant.tsv" for block in block_daf.query('chr == @chrom')['block']] for chrom in range(1,23)]
    output:
        "results/simgwas/combined_causal_variants.tsv"
    run:
        for i,x in enumerate(input):
            if i == 0:
                shell("cat %s > %s" % (x, output[0]))
            else:
                shell("cat %s | tail -n +2  >> %s" % (x, output[0]))
