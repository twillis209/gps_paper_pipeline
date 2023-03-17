include: "ukbb_export_functions.py"

localrules: compile_ukbb_gps_results, compile_ukbb_li_gps_results, compile_ukbb_ldsc_results, compile_ukbb_sumher_results, compile_ukbb_hoeffdings_results

# TODO still haven't changed the paths for the existing results to match the new snp_set/variant_set convention
rule compile_ukbb_gps_results:
    input:
        [f"results/gps/{{snp_set}}/{{variant_set}}/window_{{window}}_step_{{step}}_r2_{{r2}}/{x}_{{draws}}_permutations_gps_pvalue.tsv" for x in ukbb_trait_pairs]
    output:
        "results/gps/combined/{snp_set}/{variant_set}/window_{window}_step_{step}_r2_{r2}/{draws}_permutations/gps_pvalues.tsv"
    run:
        daf = compile_ukbb_gps_results_into_daf(input)
        daf.to_csv(output[0], sep = '\t', index = False)

rule compile_ukbb_li_gps_results:
    input:
        [f"results/gps/{{snp_set}}/{{variant_set}}/window_{{window}}_step_{{step}}_r2_{{r2}}/{x}_li_gps_pvalue.tsv" for x in ukbb_trait_pairs]
    output:
        "results/gps/combined/{snp_set}/{variant_set}/window_{window}_step_{step}_r2_{r2}/li_gps_pvalues.tsv"
    run:
        daf = compile_ukbb_li_gps_results_into_daf(input)
        daf.to_csv(output[0], sep = '\t', index = False)

rule compile_ukbb_ldsc_results:
    input:
        [f"results/ldsc/rg/ukbb/{{snp_set}}/fixed_h2_free_rg_intercept/{trait_pair}.log" for trait_pair in ukbb_trait_pairs]
    output:
        "results/ldsc/rg/ukbb/{snp_set}/fixed_h2_free_rg_intercept/compiled_results.tsv"
    run:
        daf = compile_ukbb_ldsc_results_into_daf(input)
        daf.to_csv(output[0], sep = '\t', index = False)

rule compile_ukbb_sumher_results:
    input:
        [f"results/ldak/ldak-thin/{{snp_set}}/rg/{trait_pair}.cors.full" for trait_pair in ukbb_trait_pairs]
    output:
        "results/ldak/ldak-thin/{snp_set}/rg/compiled_ukbb_sumher_results.tsv"
    run:
        daf = compile_ukbb_sumher_results_into_daf(input)
        daf.to_csv(output[0], sep = '\t', index = False)

rule compile_ukbb_hoeffdings_results:
    input:
        ["results/{{snp_set}}/{{variant_set}}/window_{{window}}_step_{{step}}_r2_{{r2}}/{x}_hoeffdings.tsv" for x in ukbb_trait_pairs]
    output:
        "results/combined/{snp_set}/{variant_set}/window_{window}_step_{step}_r2_{r2}/hoeffdings_results.tsv"
    run:
        daf = compile_ukbb_hoeffdings_results_into_daf(input)
        daf.to_csv(output[0], sep = '\t', index = False)

# TODO need to fix top maximands file paths
rule ukbb_sans_mhc:
    input:
        "results/gps/ukbb_sans_mhc/snps_only/window_1000kb_step_50_r2_0_2/compiled_top_maximands.tsv",
        "results/gps/combined/ukbb_sans_mhc/snps_only/window_1000kb_step_50_r2_0_2/3000_permutations.tsv",
        "results/gps/combined/ukbb_sans_mhc/snps_only/window_1000kb_step_50_r2_0_2/li_gps_pvalues.tsv",
        "results/ldsc/rg/ukbb/ukbb_sans_mhc/fixed_h2_free_rg_intercept/compiled_results.tsv",
        "results/combined/ukbb_sans_mhc/snps_only/window_1000kb_step_50_r2_0_2/hoeffdings_results.tsv",
        "results/ldak/ldak-thin/ukbb_sans_mhc/rg/compiled_ukbb_sumher_results.tsv"

rule ukbb_with_mhc:
    input:
        "results/gps/ukbb_sans_mhc/snps_only/window_1000kb_step_50_r2_0_2/compiled_top_maximands.tsv",
        "results/gps/combined/ukbb_sans_mhc/snps_only/window_1000kb_step_50_r2_0_2/3000_permutations.tsv",
        "results/gps/combined/ukbb_sans_mhc/snps_only/window_1000kb_step_50_r2_0_2/li_gps_pvalues.tsv",
        "results/ldsc/rg/ukbb/ukbb_sans_mhc/fixed_h2_free_rg_intercept/compiled_results.tsv",
        "results/combined/ukbb_sans_mhc/snps_only/window_1000kb_step_50_r2_0_2/hoeffdings_results.tsv",
        "results/ldak/ldak-thin/ukbb_sans_mhc/rg/compiled_ukbb_sumher_results.tsv"

# TODO get 
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
