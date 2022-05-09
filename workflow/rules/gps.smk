rule compute_gps_for_trait_pair:
    input:
      ancient("resources/{trait_A}.temp"),
      ancient("resources/{trait_B}.temp"),
      sum_stats_file = ancient("resources/pruned_sum_stats/{join}/{snp_set}/window_{window}_step_{step}/pruned_merged_sum_stats.tsv"),
    output:
        temp("results/gps/{join}/{snp_set}/window_{window}_step_{step}/{trait_A}-{trait_B}_{draws}_permutations_gps_value.tsv")
    shell:
      "workflow/scripts/gps_cpp/build/apps/computeGpsCLI -i {input.sum_stats_file} -a {wildcards.trait_A} -b {wildcards.trait_B} -c {wildcards.trait_A} -d {wildcards.trait_B} -o {output}"

rule permute_trait_pair:
    input:
      ancient("resources/{trait_A}.temp"),
      ancient("resources/{trait_B}.temp"),
      sum_stats_file = ancient("resources/pruned_sum_stats/{join}/{snp_set}/window_{window}_step_{step}/pruned_merged_sum_stats.tsv"),
    output:
        "results/gps/{join}/{snp_set,!randomised}/window_{window}_step_{step}/{draws}_permutations/{trait_A}-{trait_B}.tsv"
    threads: 8
    resources:
        mem_mb = get_mem_mb,
        time = get_permute_time,
    shell:
      "workflow/scripts/gps_cpp/build/apps/permuteTraitsCLI -i {input.sum_stats_file} -o {output} -a {wildcards.trait_A} -b {wildcards.trait_B} -c {threads} -n {wildcards.draws}"

rule fit_gev_and_compute_gps_pvalue_for_trait_pair:
    input:
        gps_file = "results/{join}/{snp_set}/window_{window}_step_{step}/{trait_A}-{trait_B}_{draws}_permutations_gps_value.tsv",
        perm_file = ancient("results/{join}/{snp_set}/window_{window}_step_{step}/{draws}_permutations/{trait_A}-{trait_B}.tsv")
    output:
        "results/gps/{join}/{snp_set}/window_{window}_step_{step}/{trait_A}-{trait_B}_{draws}_permutations_gps_pvalue.tsv"
    shell:
      "Rscript workflow/scripts/fit_gev_and_compute_gps_pvalue.R -g {input.gps_file} -p {input.perm_file} -a {wildcards.trait_A} -b {wildcards.trait_B} -o {output}"

rule collate_gps_pvalue_data:
    input:
        pvalue_files = ["results/ukbb/{snp_set}/window_{window}_step_{step}/%s_{draws}_permutations_gps_pvalue.tsv" % x for x in ukbb_trait_pairs]+
        ["results/pid_ukbb/{snp_set}/window_{window}_step_{step}/%s_{draws}_permutations_gps_pvalue.tsv" % x for x in pid_ukbb_trait_pairs],
        lookup_file = "resources/ukbb_sum_stats/trait_metadata.tsv"
    output:
        "results/gps/combined/{snp_set}/window_{window}_step_{step}/gps_pvalues_{draws}_permutations_with_labels.tsv"
    run:
        with open(output[0], 'w') as outfile:
            outfile.write(("\t".join(["trait_A", "trait_B", "gps", "n", "loc", "loc.sd", "scale", "scale.sd", "shape", "shape.sd", "pval"]))+"\n")
            for i,x in enumerate(input.pvalue_files):
                with open(x, 'r') as infile:
                    data_line = infile.readlines()[1]

                    m = re.match("results/(pid_){0,1}ukbb/%s/window_%s_step_%s/(\w+)-(\w+)_%s_permutations_gps_pvalue.tsv" % (wildcards.snp_set, wildcards.window, wildcards.step, wildcards.draws), x)

                outfile.write(("\t".join([m[2], m[3], data_line])))
        shell("Rscript workflow/scripts/add_trait_labels_to_gps_results.R -p {output} -l {input.lookup_file} -o {output}")

rule collate_gps_pvalue_data_for_new_prune:
    input:
        pvalue_files = ["results/pid_ukbb/all_pruned_snps/window_50_step_5/%s_3000_permutations_gps_pvalue.tsv" % x for x in pid_ukbb_trait_pairs]
    output:
        "results/gps/combined/all_pruned_snps/window_50_step_5/gps_pvalues_3000_permutations_with_labels.tsv"
    run:
        with open(output[0], 'w') as outfile:
            outfile.write(("\t".join(["trait_A", "trait_B", "gps", "n", "loc", "loc.sd", "scale", "scale.sd", "shape", "shape.sd", "pval"]))+"\n")
            for i,x in enumerate(input.pvalue_files):
                with open(x, 'r') as infile:
                    data_line = infile.readlines()[1]

                    m = re.match("results/(pid_){0,1}ukbb/%s/window_%s_step_%s/(\w+)-(\w+)_%s_permutations_gps_pvalue.tsv" % (wildcards.snp_set, wildcards.window, wildcards.step, wildcards.draws), x)

                outfile.write(("\t".join([m[2], m[3], data_line])))
        shell("Rscript workflow/scripts/add_trait_labels_to_gps_results.R -p {output} -l {input.lookup_file} -o {output}")

rule generate_ecdf_values_for_trait_pair:
    input:
        ancient("resources/{trait_A}.temp"),
        ancient("resources/{trait_B}.temp"),
        sum_stats_file = ancient("resources/pruned_sum_stats/{join}/{snp_set}/window_{window}_step_{step}/pruned_merged_sum_stats.tsv"),
    output:
        "results/gps/{join}/{snp_set}/window_{window}_step_{step}/{trait_A}-{trait_B}_ecdf.tsv"
    shell:
        "workflow/scripts/gps_cpp/build/apps/fitAndEvaluateEcdfsCLI -i {input.sum_stats_file} -a {wildcards.trait_A} -b {wildcards.trait_B} -o {output}"

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

#rule fit_gev_to_increasing_no_snps_for_selected_ukbb_trait_pairs:
#    input:
#        list(chain(*[["results/ukbb/%s_snps/%s_3000_permutations_estimates.tsv" % (x,y) for y in trait_pairs_for_increasing_perm_fits] for x in [10000, 50000, 100000, 200000, 300000, 400000]]))

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
