# NB: Doesn't work when I specify the chromosome no. wildcard as '{chr}'
chr1_ld_matrices = ["results/simgwas/chr1_ld_matrices/chr1_block_%d_ld_matrix.RData" % i for i in range(122)]
chr21_ld_matrices = ["results/simgwas/chr21_ld_matrices/chr21_block_%d_ld_matrix.RData" % i for i in range(23)]

"""
rule compute_block_ld_matrix:
    input:
        haplotype_file = "resources/simgwas/1000g/1000GP_Phase3_chr{ch}_with_meta_eur_common_maf.hap.gz",
        legend_file = "resources/simgwas/1000g/1000GP_Phase3_chr{ch}_eur_common_maf.legend.gz",
        block_file = "resources/ldetect/blocks.txt",
    output:
        "results/simgwas/chr{ch}_ld_matrices/chr{ch}_block_{block}_ld_matrix.RData"
    threads: 4
    resources:
        time = "4:00:00",
        mem_mb=get_mem_mb
    shell:
        "Rscript workflow/scripts/compute_block_ld_matrix.R --hap_file {input.haplotype_file} --leg_file {input.legend_file} --output_file {output} -nt {threads} --block_file {input.block_file} --chr_no {wildcards.ch} --block_no {wildcards.block}"
"""

rule simulate_sum_stats_by_ld_block:
    input:
        bim_file = ancient("resources/1000g/chr{ch}.bim"),
        haplotype_file = "resources/simgwas/1000g/1000GP_Phase3_chr{ch}_with_meta_eur_common_maf.hap.gz",
        legend_file = "resources/simgwas/1000g/1000GP_Phase3_chr{ch}_eur_common_maf.legend.gz",
        block_file = "resources/ldetect/blocks.txt",
        ld_mat_file = ancient("results/simgwas/chr{ch}_ld_matrices/chr{ch}_block_{block}_ld_matrix.RData")
    output:
        "results/simgwas/simulated_sum_stats/chr{ch}/{effect_size}/block_{block}_sum_stats.tsv.gz"
    threads: 2
    resources:
        mem_mb=get_mem_mb
    shell:
        "Rscript workflow/scripts/simulate_sum_stats_by_ld_block.R --hap_file {input.haplotype_file} --leg_file {input.legend_file} --bim_file {input.bim_file} --ld_mat_file {input.ld_mat_file} -b {input.block_file} --block_no {wildcards.block} --chr_no {wildcards.ch} --causal_variant_ind 2000 --effect_size {wildcards.effect_size} --no_controls 10000 --no_cases 10000 --no_reps 1 -o {output} -nt {threads}"

rule simulate_null_chr1_by_ld_block:
    input:
        ["results/simgwas/simulated_sum_stats/chr1/null/block_%d_sum_stats.tsv.gz" % i for i in range(122)]


"""
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
"""
