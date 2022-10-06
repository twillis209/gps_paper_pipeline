
rule compute_gps_for_trait_pair_and_write_out_intermediate_values:
    input:
        ancient("resources/{trait_A}"),
        ancient("resources/{trait_B}"),
        sum_stats_file = ancient("resources/pruned_sum_stats/{join}/{snp_set}/window_{window}_step_{step}/pruned_merged_sum_stats.tsv"),
    output:
        temp("results/gps/{join}/{snp_set}/window_{window}_step_{step}/{trait_A}-{trait_B}_intermediates.tsv")
    threads: 12
    resources:
        runtime = 30
    group: "gps"
    shell:
        "workflow/scripts/gps_cpp/build/apps/fitAndEvaluateEcdfsCLI -i {input.sum_stats_file} -a {wildcards.trait_A} -b {wildcards.trait_B} -n {threads} -o {output}"

rule annotate_intermediate_gps_output:
    input:
        intermediates_file = "results/gps/{join}/{snp_set}/window_{window}_step_{step}/{trait_A}-{trait_B}_intermediates.tsv",
        sum_stats_file = ancient("resources/pruned_sum_stats/{join}/{snp_set}/window_{window}_step_{step}/pruned_merged_sum_stats.tsv")
    output:
        "results/gps/{join}/{snp_set}/window_{window}_step_{step}/{trait_A}-{trait_B}_intermediates_annot.tsv"
    group: "gps"
    threads: 4
    script: "../../scripts/annotate_intermediate_gps_output.R"

rule plot_denominator_heatmap:
    input:
        "results/gps/{join}/{snp_set}/window_{window}_step_{step}/{trait_A}-{trait_B}_ecdf.tsv"
    output:
        "results/plots/{join}/gps_heatmaps/{trait_A}-{trait_B}.png"
    shell:
        "workflow/scripts/plot_gps_denom_heatmap.R -i {input} -o {output}"

rule compile_top_maximands_for_ukbb_traits:
    input:
        annot_files = [f"results/gps/ukbb/{{snp_set}}/window_1000kb_step_50/{trait_pair}_intermediates_annot.tsv" for trait_pair in ukbb_trait_pairs],
        pvalue_files = [f"results/gps/ukbb/{{snp_set}}/window_1000kb_step_50/{trait_pair}_3000_permutations_gps_pvalue.tsv" for trait_pair in ukbb_trait_pairs]
    output:
        "results/gps/ukbb/{snp_set}/window_1000kb_step_50/compiled_top_maximands.tsv"
    threads: 4
    script: "../../scripts/compile_top_maximands_for_ukbb_traits.R"
