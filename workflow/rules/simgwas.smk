rule compute_ld_matrices:
    input:
        bim_file = ancient("resources/1000g/chr21.bim"),
        haplotype_file = "resources/simgwas/1000g/1000GP_Phase3_chr21_with_meta_eur.hap.gz",
        legend_file = "resources/simgwas/1000g/1000GP_Phase3_chr21.legend.gz",
        block_file = "resources/ldetect/blocks.txt"
    output:
        "results/simgwas/chr21_block_ld_matrices_alt.RData"
    benchmark:
        "results/benchmark/simgwas/chr21_block_ld_matrices_benchmark.tsv"
    threads: 8
    resources:
        time = "10:00:00"
    shell:
        "Rscript workflow/scripts/compute_ld_matrices.R --hap_file {input.haplotype_file} --leg_file {input.legend_file} --bim_file {input.bim_file} -o {output} -nt {threads} --block_file {input.block_file}"

rule run_simgwas:
    input:
     bim_file = ancient("resources/1000g/chr21.bim"),
     haplotype_file = "resources/simgwas/1000g/1000GP_Phase3_chr21_with_meta_eur.hap.gz",
     legend_file = "resources/simgwas/1000g/1000GP_Phase3_chr21.legend.gz",
     block_file = "resources/ldetect/blocks.txt",
     ld_mats_file = "results/simgwas/chr21_block_ld_matrices.RData"
    output:
     "results/simgwas/simulated_sum_stats/{no_snps}_snps_{no_cv}_cv_sum_stats.tsv.gz"
    benchmark:
     "results/benchmark/simgwas/{no_snps}_snps_{no_cv}_cv_benchmark.tsv"
    threads: 8
    resources:
        time = "10:00:00"
    shell:
        "Rscript workflow/scripts/simulate_sum_stats.R --hap_file {input.haplotype_file} --leg_file {input.legend_file} --bim_file {input.bim_file} --block-file {input.block_file} --ld_mats_file {input.ld_mats_file} --no_snps {wildcards.no_snps} --no_causal_variants {wildcards.no_cv} --odds_ratios 1.3 --no_controls 10000 --no_cases 10000 --no_reps 1 -o {output} -nt {threads}"

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
