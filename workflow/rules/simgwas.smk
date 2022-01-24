# TODO REPLACE AFR with EUR files, amend simulate_sum_stats to cope with metadata rows
rule run_simgwas:
    input:
     haplotype_file = "resources/simgwas/1000g/chr21_afr.hap.gz",
     legend_file = "resources/simgwas/1000g/1000GP_Phase3_chr21.legend.gz"
    output:
     "results/simgwas/simulated_sum_stats/{no_snps}_snps_{no_cv}_cv_sum_stats.tsv.gz"
    threads: 4
    shell:
        "Rscript workflow/scripts/simulate_sum_stats.R --hap_file {input.haplotype_file} --leg_file {input.legend_file} --no_snps {wildcards.no_snps} --no_causal_variants {wildcards.no_cv} --odds_ratios 3 --no_controls 10000 --no_cases 10000 --no_reps 1 -o {output} -nt {threads}"

        # TODO write script to downsample output, maybe just add a column
rule prune_simulated_sum_stats:
    input:
        pruned_range_file = ancient("resources/plink_ranges/{join}/pruned_ranges/window_{window}_step_{step}/all.prune.in")
        # TODO put this in a directory labelled by pruning parameter values
        sum_stats_file = "results/simgwas/pruned_simulated_sum_stats/window_{window}_step_{step}/{no_snps}_snps_{no_cv}_cv_sum_stats.tsv.gz"
    output:
        pruned_haplotype_file = "resources/simgwas/1000g/chr21_afr.hap.gz",

