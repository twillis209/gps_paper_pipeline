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

rule compile_intermediate_values_across_s400:
    input:
        "results/gps/simgwas/400_reps/randomised/500_10000_500_10000/s400_s400_s0/window_1000kb_step_50_r2_0_2/seed_1_tags_1-2_gps_intermediates.tsv",
        "results/gps/simgwas/400_reps/randomised/1000_10000_1000_10000/s400_s400_s0/window_1000kb_step_50_r2_0_2/seed_12201_tags_1-2_gps_intermediates.tsv",
        "results/gps/simgwas/400_reps/randomised/5000_10000_5000_10000/s400_s400_s0/window_1000kb_step_50_r2_0_2/seed_36601_tags_1-2_gps_intermediates.tsv",
        "results/gps/simgwas/400_reps/randomised/10000_10000_10000_10000/s400_s400_s0/window_1000kb_step_50_r2_0_2/seed_48801_tags_1-2_gps_intermediates.tsv",
        "results/gps/simgwas/400_reps/randomised/100000_100000_100000_100000/s400_s400_s0/window_1000kb_step_50_r2_0_2/seed_73201_tags_1-2_gps_intermediates.tsv"

rule run_m20_one_chrom_r2_sweep:
    input:
        get_one_chrom_simulation_done_files("results/simgwas/one_chrom_simulation_parameters.tsv", reps = 400)

rule compile_li_gps_ukbb_results:
    input:
        "results/gps/combined/ukbb_with_mhc/all/window_1000kb_step_50_r2_0_2/li_gps_pvalues.tsv",
        "results/gps/combined/ukbb_sans_mhc/all/window_1000kb_step_50_r2_0_2/li_gps_pvalues.tsv"

