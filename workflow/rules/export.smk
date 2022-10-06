rule ukbb_sans_mhc:
    input:
        "results/gps/ukbb/sans_mhc/window_1000kb_step_50/compiled_top_maximands.tsv",
        "results/gps/combined/sans_mhc/window_1000kb_step_50/gps_pvalues_3000_permutations.tsv",
        "results/ldsc/rg/ukbb/sans_mhc/fixed_h2_free_rg_intercept/compiled_results.tsv",
        "results/combined/sans_mhc/window_1000kb_step_50/hoeffdings_results.tsv",
        "results/ldak/ldak-thin/ukbb/sans_mhc/rg/compiled_ukbb_sumher_results.tsv"

rule ukbb_with_mhc:
    input:
        "results/gps/ukbb/all/window_1000kb_step_50/compiled_top_maximands.tsv",
        "results/gps/combined/all/window_1000kb_step_50/gps_pvalues_3000_permutations.tsv",
        "results/ldsc/rg/ukbb/all/fixed_h2_free_rg_intercept/compiled_results.tsv",
        "results/combined/all/window_1000kb_step_50/hoeffdings_results.tsv",
        "results/ldak/ldak-thin/ukbb/all/rg/compiled_ukbb_sumher_results.tsv"
