rule run_new_gps_algorithm:
    input:
        "results/gps/simgwas/400_reps/randomised/10000_10000_10000_10000/s0_s0_s0/window_1000kb_step_50/seed_600000_tags_1-2_li_gps_pvalue.tsv"

rule run_new_pp_algorithm:
    input:
        "results/gps/simgwas/400_reps/randomised/10000_10000_10000_10000/s0_s0_s0/window_1000kb_step_50/3000_permutations/seed_600000_tags_1-2.tsv"
