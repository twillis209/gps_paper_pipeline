from scipy.stats import chi2

rule process_combined_simgwas_sum_stats_for_chrom:
    input:
        "results/simgwas/simulated_sum_stats/per_chrom_sum_stats/{no_reps}_reps/randomised/{chr}/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/seed_{seed}_sum_stats_{pair_label}_tag_{tag}_of_{tag_A}-{tag_B}.tsv.gz"
    output:
        "results/simgwas/simulated_sum_stats/per_chrom_sum_stats/{no_reps}_reps/randomised/{chr}/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/seed_{seed}_sum_stats_{pair_label}_tag_{tag}_of_{tag_A}-{tag_B}.assoc"
    params:
        z_colname = lambda wildcards: f'zsim.{wildcards.tag}',
        chr_colname = 'chr',
        bp_colname = 'position',
        a1_colname = 'a1',
        a2_colname = 'a0',
        sample_size = lambda wildcards: int(wildcards.ncases_A)+int(wildcards.ncontrols_A)
    threads: 2
    resources:
        runtime = 5
    priority: 1
    group: "ldsc_hoeffding_sumher_gps_sans_permutation"
    script:
        "../../scripts/process_combined_simgwas_sum_stats_for_sumher.R"

rule estimate_rg_with_ldak_thin_for_simgwas_for_chrom:
    input:
        tagging_file = "results/ldak/ldak-thin/taggings/all/{chr}.tagging",
        sum_stats_file_A = "results/simgwas/simulated_sum_stats/per_chrom_sum_stats/{no_reps}_reps/randomised/{chr}/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/seed_{seed}_sum_stats_A_tag_{tag_A}_of_{tag_A}-{tag_B}.assoc",
        sum_stats_file_B = "results/simgwas/simulated_sum_stats/per_chrom_sum_stats/{no_reps}_reps/randomised/{chr}/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/seed_{seed}_sum_stats_B_tag_{tag_B}_of_{tag_A}-{tag_B}.assoc"
    output:
#        progress_file = "results/ldak/ldak-thin/simgwas/{no_reps}_reps/randomised/rg/{chr}/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/seed_{seed}_tags_{tag_A}-{tag_B}.progress",
        cors_file = "results/ldak/ldak-thin/simgwas/{no_reps}_reps/randomised/rg/{chr}/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/seed_{seed}_tags_{tag_A}-{tag_B}.cors"
#        cors_full_file = "results/ldak/ldak-thin/simgwas/{no_reps}_reps/randomised/rg/{chr}/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/seed_{seed}_tags_{tag_A}-{tag_B}.cors.full",
#        labels_file = "results/ldak/ldak-thin/simgwas/{no_reps}_reps/randomised/rg/{chr}/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/seed_{seed}_tags_{tag_A}-{tag_B}.labels",
#        overlap_file = "results/ldak/ldak-thin/simgwas/{no_reps}_reps/randomised/rg/{chr}/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/seed_{seed}_tags_{tag_A}-{tag_B}.overlap"
    log:
        log_file = "results/ldak/ldak-thin/simgwas/{no_reps}_reps/randomised/rg/{chr}/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/seed_{seed}_tags_{tag_A}-{tag_B}.log"
    params:
        output_stem = "results/ldak/ldak-thin/simgwas/{no_reps}_reps/randomised/rg/{chr}/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/seed_{seed}_tags_{tag_A}-{tag_B}"
    resources:
        runtime = 5
    priority: 1
    group: "ldsc_hoeffding_sumher_gps_sans_permutation"
    shell:
        """
        $ldakRoot/ldak --sum-cors {params.output_stem} --tagfile {input.tagging_file} --summary {input.sum_stats_file_A} --summary2 {input.sum_stats_file_B} --allow-ambiguous YES --check-sums NO --cutoff 0.01 > {log.log_file}
        """
