rule compute_gps_for_trait_pair_and_write_out_intermediate_values:
    input:
        ancient("resources/{trait_A}"),
        ancient("resources/{trait_B}"),
        sum_stats_file = ancient("resources/ukbb_sum_stats/{snp_set}/{variant_set}/window_{window}_step_{step}_r2_{r2}/pruned_merged_sum_stats.tsv"),
    output:
        temp("results/gps/{snp_set}/{variant_set}/window_{window}_step_{step}_r2_{r2}/{trait_A}-{trait_B}_intermediates.tsv")
    threads: 4
    resources:
        runtime = 2,
        mem_mb = get_mem_mb
    group: "gps"
    shell:
        "workflow/scripts/gps_cpp/build/apps/fitAndEvaluateEcdfsCLI -i {input.sum_stats_file} -a {wildcards.trait_A} -b {wildcards.trait_B} -n {threads} -o {output} -f pp"

rule annotate_intermediate_gps_output:
    input:
        intermediates_file = "results/gps/{snp_set}/{variant_set}/window_{window}_step_{step}_r2_{r2}/{trait_A}-{trait_B}_intermediates.tsv",
        sum_stats_file = ancient("resources/ukbb_sum_stats/{snp_set}/{variant_set}/window_{window}_step_{step}_r2_{r2}/pruned_merged_sum_stats.tsv")
    output:
        "results/gps/{snp_set}/{variant_set}/window_{window}_step_{step}_r2_{r2}/{trait_A}-{trait_B}_intermediates_annot.tsv"
    threads: 4
    resources:
        runtime = 2
    group: "gps"
    script: "../../scripts/gps/annotate_intermediate_gps_output.R"

rule plot_denominator_heatmap:
    input:
        "results/gps/{snp_set}/{variant_set}/window_{window}_step_{step}_r2_{r2}/{trait_A}-{trait_B}_ecdf.tsv"
    output:
        "results/gps/{snp_set}/{variant_set}/gps_heatmaps/{trait_A}-{trait_B}.png"
    group: "gps"
    shell:
        "workflow/scripts/gps/plot_gps_denom_heatmap.R -i {input} -o {output}"

rule compile_top_maximands_for_ukbb_traits:
    input:
        annot_files = [f"results/gps/{{snp_set}}/{{variant_set}}/window_{{window}}_step_{{step}}_r2_{{r2}}/{trait_pair}_intermediates_annot.tsv" for trait_pair in ukbb_trait_pairs],
        pvalue_files = [f"results/gps/{{snp_set}}/{{variant_set}}/window_{{window}}_step_{{step}}_r2_{{r2}}/{trait_pair}_3000_permutations_gps_pvalue.tsv" for trait_pair in ukbb_trait_pairs]
    output:
        "results/gps/combined/{snp_set}/{variant_set}/window_{window}_step_{step}_r2_{r2}/compiled_top_maximands.tsv"
    threads: 4
    group: "gps"
    script: "../../scripts/gps/compile_top_maximands_for_ukbb_traits.R"
