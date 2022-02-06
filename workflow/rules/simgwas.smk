# NB: Doesn't work when I specify the chromosome no. wildcard as '{chr}'
rule compute_ld_matrices:
    input:
        haplotype_file = "resources/simgwas/1000g/1000GP_Phase3_chr{ch}_with_meta_eur_common_maf.hap.gz",
        legend_file = "resources/simgwas/1000g/1000GP_Phase3_chr{ch}_eur_common_maf.legend.gz",
        block_file = "resources/ldetect/blocks.txt",
    output:
        "results/simgwas/chr{ch}_block_ld_matrices.RData"
    threads: 20
    resources:
        time = "10:00:00"
    shell:
        "Rscript workflow/scripts/compute_ld_matrices.R --hap_file {input.haplotype_file} --leg_file {input.legend_file} --output_file {output} -nt {threads} --block_file {input.block_file} --chr_no {wildcards.ch}"

rule compute_block_ld_matrix:
    input:
        haplotype_file = "resources/simgwas/1000g/1000GP_Phase3_chr{ch}_with_meta_eur_common_maf.hap.gz",
        legend_file = "resources/simgwas/1000g/1000GP_Phase3_chr{ch}_eur_common_maf.legend.gz",
        block_file = "resources/ldetect/blocks.txt",
    output:
        "results/simgwas/chr{ch}_ld_matrices/chr{ch}_block_{block}_ld_matrix.RData"
    threads: 4
    resources:
        time = "2:00:00"
    shell:
        "Rscript workflow/scripts/compute_block_ld_matrix.R --hap_file {input.haplotype_file} --leg_file {input.legend_file} --output_file {output} -nt {threads} --block_file {input.block_file} --chr_no {wildcards.ch} --block_no {wildcards.block}"

rule run_simgwas:
    input:
     bim_file = ancient("resources/1000g/chr{ch}.bim"),
     haplotype_file = "resources/simgwas/1000g/1000GP_Phase3_chr{ch}_with_meta_eur.hap.gz",
     legend_file = "resources/simgwas/1000g/1000GP_Phase3_chr{ch}.legend.gz",
     block_file = "resources/ldetect/blocks.txt",
     ld_mats_file = "results/simgwas/chr{ch}_block_ld_matrices.RData"
    output:
        "results/simgwas/simulated_sum_stats/chr{ch}_{no_blocks}_blocks_sum_stats.tsv.gz"
    benchmark:
        "results/benchmark/simgwas/chr{ch}_{no_blocks}_blocks_benchmark.tsv"
    resources:
        threads = lambda wildcards, attempt: wildcards.no_blocks
    shell:
        "Rscript workflow/scripts/simulate_sum_stats.R --hap_file {input.haplotype_file} --leg_file {input.legend_file} --bim_file {input.bim_file} --ld_mats_file {input.ld_mats_file} -b {input.block_file} --no_blocks {wildcards.no_blocks} --chr_no {wildcards.ch} --no_controls 10000 --no_cases 10000 --no_reps 1 -o {output} -nt {threads}"

        # TODO fix this
rule prune_simulated_sum_stats:
    input:
      bim_file = ancient("resources/1000g/euro/qc/chr21_qc.bim"),
      pruned_range_file = ancient("resources/plink_ranges/ukbb/pruned_ranges/window_{window}_step_{step}/all.prune.in"),
      sum_stats_file = "results/simgwas/simulated_sum_stats/{no_snps}_snps_{no_cv}_cv_sum_stats.tsv.gz"
    output:
      "results/simgwas/pruned_simulated_sum_stats/window_{window}_step_{step}/{no_snps}_snps_{no_cv}_cv_sum_stats.tsv.gz"
    threads: 8
    shell:
      "Rscript workflow/scripts/prune_sim_sum_stats.R --sum_stats_file {input.sum_stats_file} --bim_file {input.bim_file} --prune_file {input.pruned_range_file} -o {output} -nt {threads}"
