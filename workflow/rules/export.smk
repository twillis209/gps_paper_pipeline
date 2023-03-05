rule ukbb_sans_mhc:
    input:
        "results/gps/ukbb_sans_mhc/snps_only/window_1000kb_step_50_r2_0_2/compiled_top_maximands.tsv",
        "results/gps/combined/ukbb_sans_mhc/snps_only/window_1000kb_step_50_r2_0_2/gps_pvalues_3000_permutations.tsv",
        "results/ldsc/rg/ukbb_sans_mhc/fixed_h2_free_rg_intercept/compiled_results.tsv",
        "results/combined/sans_mhc/window_1000kb_step_50_r2_0_2/hoeffdings_results.tsv",
        "results/ldak/ldak-thin/ukbb/sans_mhc/rg/compiled_ukbb_sumher_results.tsv"

rule ukbb_with_mhc:
    input:
        "results/gps/ukbb/all/window_1000kb_step_50_r2_0_2/compiled_top_maximands.tsv",
        "results/gps/combined/all_pruned_snps/window_1000kb_step_50_r2_0_2/gps_pvalues_3000_permutations.tsv",
        "results/ldsc/rg/ukbb/all/fixed_h2_free_rg_intercept/compiled_results.tsv",
        "results/combined/all_pruned_snps/window_1000kb_step_50_r2_0_2/hoeffdings_results.tsv",
        "results/ldak/ldak-thin/ukbb/all/rg/compiled_ukbb_sumher_results.tsv"

rule run_simulated_analyses:
    input:
        "results/s400_400_reps.done",
        "results/m25_400_reps.done",
        "results/m50_400_reps.done",
        "results/s200-m25_400_reps.done"

rule compile_all_simulation_results:
    input:
        [[f"results/ldsc/simgwas/400_reps/randomised/compiled_{x}_results.tsv",
        f"results/ldak/ldak-thin/simgwas/400_reps/randomised/rg/compiled_{x}_results.tsv",
        f"results/hoeffdings/simgwas/400_reps/randomised/compiled_{x}_results.tsv",
        f"results/gps/simgwas/400_reps/randomised/compiled_{x}_results.tsv",
         f"results/ldsc/simgwas/400_reps/randomised/compiled_{x}_theo_rg.tsv"] for x in ['s400', 'm25', 'm50', 's200-m25']]
    output:
        "results/all_sim.tsv"
    shell:
        """
        touch {output}
        """
rule compile_real_results:
    input:
        "results/gps/combined/sans_mhc/window_1000kb_step_50_r2_0_2/gps_pvalues_3000_permutations.tsv",
        "results/ldsc/rg/ukbb/sans_mhc/fixed_h2_free_rg_intercept/compiled_results.tsv",
        "results/combined/sans_mhc/window_1000kb_step_50_r2_0_2/hoeffdings_results.tsv",
        "results/ldak/ldak-thin/ukbb/sans_mhc/rg/compiled_ukbb_sumher_results.tsv",
        "results/gps/combined/all/window_1000kb_step_50_r2_0_2/gps_pvalues_3000_permutations.tsv",
        "results/ldsc/rg/ukbb/all/fixed_h2_free_rg_intercept/compiled_results.tsv",
        "results/combined/all_pruned_snps/window_1000kb_step_50_r2_0_2/hoeffdings_results.tsv",
        "results/ldak/ldak-thin/ukbb/all/rg/compiled_ukbb_sumher_results.tsv"
    script: "../scripts/compile_ukbb_results.R"

rule compile_intermediate_values_across_s400:
    input:
        "results/gps/simgwas/400_reps/randomised/500_10000_500_10000/s400_s400_s0/window_1000kb_step_50_r2_0_2/seed_1_tags_1-2_gps_intermediates.tsv",
        "results/gps/simgwas/400_reps/randomised/1000_10000_1000_10000/s400_s400_s0/window_1000kb_step_50_r2_0_2/seed_12201_tags_1-2_gps_intermediates.tsv",
        "results/gps/simgwas/400_reps/randomised/5000_10000_5000_10000/s400_s400_s0/window_1000kb_step_50_r2_0_2/seed_36601_tags_1-2_gps_intermediates.tsv",
        "results/gps/simgwas/400_reps/randomised/10000_10000_10000_10000/s400_s400_s0/window_1000kb_step_50_r2_0_2/seed_48801_tags_1-2_gps_intermediates.tsv",
        "results/gps/simgwas/400_reps/randomised/100000_100000_100000_100000/s400_s400_s0/window_1000kb_step_50_r2_0_2/seed_73201_tags_1-2_gps_intermediates.tsv"

rule test_one_chrom_r2_sweep:
    input:
        [f"results/simgwas/done/400_reps/randomised/chr1/10000_10000_10000_10000/m20_m20_m0/window_1000kb_step_50_r2_{x}/seed_1_tags_1-2.done" for x in ["0_2", "0_5", "0_8"]]


