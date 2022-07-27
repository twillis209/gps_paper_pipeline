rule compute_gps_for_trait_pair:
    input:
      ancient("resources/ukbb_sum_stats/{trait_A}.done"),
      ancient("resources/ukbb_sum_stats/{trait_B}.done"),
      sum_stats_file = ancient("resources/pruned_sum_stats/{join}/{snp_set}/window_{window}_step_{step}/pruned_merged_sum_stats.tsv"),
    output:
        temp("results/{join}/{snp_set}/window_{window}_step_{step}/{trait_A}-{trait_B}_{draws}_permutations_gps_value.tsv")
    log:
        "results/{join}/{snp_set}/window_{window}_step_{step}/{trait_A}-{trait_B}_{draws}_permutations_gps_value.log"
    params:
        no_of_perturbations = 1
    group: "gps"
    resources:
        runtime = 10
    shell:
      "workflow/scripts/gps_cpp/build/apps/computeGpsCLI -i {input.sum_stats_file} -a {wildcards.trait_A} -b {wildcards.trait_B} -c {wildcards.trait_A} -d {wildcards.trait_B} -p {params.no_of_perturbations} -l -o {output} -g {log}"

rule compute_gps_for_trait_pair_with_naive_ecdf_algo:
    input:
        ancient("resources/ukbb_sum_stats/{trait_A}.done"),
        ancient("resources/ukbb_sum_stats/{trait_B}.done"),
        sum_stats_file = ancient("resources/pruned_sum_stats/{join}/{snp_set}/window_{window}_step_{step}/pruned_merged_sum_stats.tsv"),
    output:
        "results/{join}/{snp_set}/window_{window}_step_{step}/{trait_A}-{trait_B}_naive_gps_value.tsv"
    params:
        no_of_pert_iterations = 0
    threads: 10
    resources:
        runtime = 30
    group: "gps"
    shell:
        "workflow/scripts/gps_cpp/build/apps/computeGpsCLI -i {input.sum_stats_file} -a {wildcards.trait_A} -b {wildcards.trait_B} -c {wildcards.trait_A} -d {wildcards.trait_B} -n {threads} -p {params.no_of_pert_iterations} -o {output}"

rule permute_trait_pair:
    input:
      ancient("resources/ukbb_sum_stats/{trait_A}.done"),
      ancient("resources/ukbb_sum_stats/{trait_B}.done"),
      sum_stats_file = ancient("resources/pruned_sum_stats/{join}/{snp_set}/window_{window}_step_{step}/pruned_merged_sum_stats.tsv"),
    output:
        "results/{join}/{snp_set,all_pruned_snps|sans_mhc}/window_{window}_step_{step}/{draws}_permutations/{trait_A}-{trait_B}.tsv"
    params:
        no_of_perturbations = 1
    threads: 8
    resources:
        mem_mb = get_mem_mb,
        runtime = get_permute_time,
    group: "gps"
    shell:
      "workflow/scripts/gps_cpp/build/apps/permuteTraitsCLI -i {input.sum_stats_file} -o {output} -a {wildcards.trait_A} -b {wildcards.trait_B} -c {threads} -n {wildcards.draws} -p {params.no_of_perturbations}"

rule fit_gev_and_compute_gps_pvalue_for_trait_pair:
    input:
        gps_file = "results/{join}/{snp_set}/window_{window}_step_{step}/{trait_A}-{trait_B}_{draws}_permutations_gps_value.tsv",
        perm_file = ancient("results/{join}/{snp_set}/window_{window}_step_{step}/{draws}_permutations/{trait_A}-{trait_B}.tsv")
    output:
        "results/gps/{join}/{snp_set}/window_{window}_step_{step}/{trait_A}-{trait_B}_{draws}_permutations_gps_pvalue.tsv"
    resources:
        runtime = 5
    group: "gps"
    shell:
      "Rscript workflow/scripts/fit_gev_and_compute_gps_pvalue.R -g {input.gps_file} -p {input.perm_file} -a {wildcards.trait_A} -b {wildcards.trait_B} -o {output}"

rule collate_gps_pvalue_data:
    input:
        pvalue_files = ["results/gps/ukbb/{snp_set}/window_{window}_step_{step}/%s_{draws}_permutations_gps_pvalue.tsv" % x for x in ukbb_trait_pairs],
        lookup_file = "resources/ukbb_sum_stats/trait_metadata.tsv"
    output:
        combined_pvalue_file = "results/gps/combined/{snp_set}/window_{window}_step_{step}/gps_pvalues_{draws}_permutations_with_labels.tsv"
    resources:
        runtime = 5
    run:
        print(output)
        with open(output.combined_pvalue_file, 'w') as outfile:
            outfile.write(("\t".join(["trait_A", "trait_B", "gps", "n", "loc", "loc.sd", "scale", "scale.sd", "shape", "shape.sd", "pval"]))+"\n")
            for i,x in enumerate(input.pvalue_files):
                with open(x, 'r') as infile:
                    line = infile.readline()
                    line = infile.readline()

                    m = re.match("results/gps/ukbb/%s/window_%s_step_%s/(\w+)-(\w+)_%s_permutations_gps_pvalue.tsv" % (wildcards.snp_set, wildcards.window, wildcards.step, wildcards.draws), x)

                outfile.write(("\t".join([m[1], m[2], line])))
        shell("Rscript workflow/scripts/add_trait_labels_to_gps_results.R -p {output} -l {input.lookup_file} -o {output}")

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
    script: "../scripts/annotate_intermediate_gps_output.R"

rule plot_denominator_heatmap:
    input:
        "results/gps/{join}/{snp_set}/window_{window}_step_{step}/{trait_A}-{trait_B}_ecdf.tsv"
    output:
        "results/plots/{join}/gps_heatmaps/{trait_A}-{trait_B}.png"
    shell:
        "workflow/scripts/plot_gps_denom_heatmap.R -i {input} -o {output}"

rule plot_null_dists_to_compare:
    input:
        perm_file = "results/{join}/{snp_set}/window_{window}_step_{step}/{draws}_permutations/{trait_A}-{trait_B}.tsv",
        fit_file = "results/{join}/{snp_set}/window_{window}_step_{step}/{trait_A}-{trait_B}_{draws}_permutations_gps_pvalue.tsv"
    output:
        exp1_null = "results/plots/{join}/{snp_set}/window_{window}_step_{step}/{trait_A}-{trait_B}_{draws}_permutations_exp1_null_dist.png",
        gev_null = "results/plots/{join}/{snp_set}/window_{window}_step_{step}/{trait_A}-{trait_B}_{draws}_permutations_gev_null_dist.png",
        exp1_gev_combined = "results/plots/{join}/{snp_set}/window_{window}_step_{step}/{trait_A}-{trait_B}_{draws}_permutations_null_dists.png"
    shell:
      "Rscript workflow/scripts/plot_null_dists_to_compare.R -f {input.fit_file} -p {input.perm_file} --exp1_null {output.exp1_null} --gev_null {output.gev_null} --exp1_gev_combined {output.exp1_gev_combined}"

rule plot_all_null_dists_to_compare:
    input:
        ["results/plots/ukbb/all_pruned_snps/window_{window}_step_{step}/%s_3000_permutations_null_dists.png" % x for x in ukbb_trait_pairs]+
        ["results/plots/pid_ukbb/all_pruned_snps/window_{window}_step_{step}/%s_3000_permutations_null_dists.png" % x for x in pid_ukbb_trait_pairs]

rule plot_gof_plots:
    input:
        perm_file = "results/{join}/{snp_set}/window_{window}_step_{step}/{draws}_permutations/{trait_A}-{trait_B}.tsv",
        fit_file = "results/{join}/{snp_set}/window_{window}_step_{step}/{trait_A}-{trait_B}_{draws}_permutations_gps_pvalue.tsv"
    output:
        "results/plots/{join}/{snp_set}/window_{window}_step_{step}/{trait_A}-{trait_B}_{draws}_permutations_gof_plots.png"
    shell:
      "Rscript workflow/scripts/plot_gev_gof_plots.R -f {input.fit_file} -p {input.perm_file} -l {wildcards.trait_A}-{wildcards.trait_B} -o {output}"

rule plot_all_gof_plots:
    input:
        ["results/plots/ukbb/all_pruned_snps/window_1000kb_step_50/%s_3000_permutations_gof_plots.png" % x for x in ukbb_trait_pairs]+
        ["results/plots/pid_ukbb/all_pruned_snps/window_1000kb_step_50/%s_3000_permutations_gof_plots.png" % x for x in pid_ukbb_trait_pairs]

# TODO remove me
rule write_freq_map:
    input:
        ancient("resources/{trait}.temp"),
        sum_stats_file = ancient("resources/pruned_sum_stats/{join}/{snp_set}/window_{window}_step_{step}/pruned_merged_sum_stats.tsv"),
    output:
        "results/{join}/{snp_set}/window_{window}_step_{step}/freq_map/add_epsilon/{trait}.tsv"
    shell:
        "workflow/scripts/gps_cpp/build/apps/writeFreqMapCLI -i {input.sum_stats_file} -a {wildcards.trait} -b {output}"

# TODO remove me
rule freq_maps:
    input:
        ["results/ukbb/all_pruned_snps/window_{window}_step_{step}/freq_map/add_epsilon/%s.tsv" % x for x in ukbb_trait_codes]

rule fit_gev_to_increasing_n:
    input:
        "results/{join}/{snp_set}/window_{window}_step_{step}/10000_permutations/{trait_A}-{trait_B}.tsv"
    output:
        "results/{join}/{snp_set}/window_{window}_step_{step}/{trait_A}-{trait_B}_1000-10000_permutations_estimates.tsv"
    shell:
      "Rscript workflow/scripts/fit_gev_to_increasing_n.R -a {wildcards.trait_A} -b {wildcards.trait_B} -p {input} -n 500 1000 2000 3000 4000 5000 6000 7000 8000 9000 10000 -o {output}"

rule fit_gev_to_increasing_n_for_selected_ukbb_trait_pairs:
    input:
        ["results/ukbb/all_pruned_snps/window_1000kb_step_50/%s_1000-10000_permutations_estimates.tsv" % x for x in trait_pairs_for_increasing_perm_fits]
    output:
        "results/ukbb/all_pruned_snps/window_1000kb_step_50/1000-10000_combined_permutations_estimates.tsv"
    shell:
        """
        echo -e 'trait_A\ttrait_B\tn\tloc\tloc.sd\tscale\tscale.sd\tshape\tshape.sd' >> {output}
        for x in {input}; do
        cat <(tail -n +2 $x) >> {output}
        done
        """

rule plot_gev_estimates_for_increasing_n:
    input:
        "results/{join}/{snp_set}/window_{window}_step_{step}/{trait_A}-{trait_B}_1000-10000_permutations_estimates.tsv"
    output:
        "results/plots/{join}/{snp_set}/window_{window}_step_{step}/{trait_A}-{trait_B}_1000-10000_permutations_estimates.png"
    shell:
        "Rscript workflow/scripts/plot_gev_estimates_for_increasing_n.R -f {input} -o {output}"

rule plot_gev_estimates_for_increasing_n_for_selected_ukbb_trait_pairs:
    input:
        ["results/plots/ukbb/all_pruned_snps/window_1000kb_step_50/%s_1000-10000_permutations_estimates.png" % x for x in trait_pairs_for_increasing_perm_fits]

rule plot_gev_estimates_for_increasing_no_snps:
    input:
        ["results/{join}/%s_snps/window_{window}_step_{step}/{trait_A}-{trait_B}_{draws}_permutations_estimates.tsv" % x for x in [10000, 50000, 100000, 200000, 300000, 400000]]
    output:
        "results/plots/{join}/window_{window}_step_{step}/{trait_A}-{trait_B}_{draws}_permutations_variable_no_snps_estimates.png"
    params:
        fit_file_string = " ".join(["results/{join}/%s_snps/window_{window}_step_{step}/{trait_A}-{trait_B}_{draws}_permutations_estimates.tsv" % x for x in [10000, 50000, 100000, 200000, 300000, 400000]])
    shell:
      "Rscript workflow/scripts/plot_gev_estimates_for_increasing_no_snps.R -i {params.fit_file_string} -n 10000 50000 100000 200000 300000 400000 -o {output}"

rule plot_gev_estimates_for_increasing_no_snps_for_selected_ukbb_traits:
    input:
        ["results/plots/ukbb/window_1000kb_step_50/%s_3000_permutations_variable_no_snps_estimates.png" % x for x in trait_pairs_for_increasing_perm_fits]

rule compile_naive_gps_values:
    input:
        [f"results/ukbb/{{snp_set}}/window_1000kb_step_50/{trait_pair}_naive_gps_value.tsv" for trait_pair in ukbb_trait_pairs]
    output:
        "results/ukbb/{snp_set}/window_1000kb_step_50/compiled_naive_gps_values.tsv"
    shell:
        """
        echo -e "Trait_A\tTrait_B\tGPS" >{output}

        for x in {input}; do
            tail -n 1 $x >>{output}
        done
        """

rule compile_top_maximands_for_ukbb_traits:
    input:
        annot_files = [f"results/gps/ukbb/{{snp_set}}/window_1000kb_step_50/{trait_pair}_intermediates_annot.tsv" for trait_pair in ukbb_trait_pairs],
        pvalue_files = [f"results/gps/ukbb/{{snp_set}}/window_1000kb_step_50/{trait_pair}_3000_permutations_gps_pvalue.tsv" for trait_pair in ukbb_trait_pairs]
    output:
        "results/gps/ukbb/{snp_set}/window_1000kb_step_50/compiled_top_maximands.tsv"
    threads: 4
    script: "../scripts/compile_top_maximands_for_ukbb_traits.R"
