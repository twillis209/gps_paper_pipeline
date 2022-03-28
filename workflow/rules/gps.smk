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
        "results/gps/{join}/{snp_set}/window_{window}_step_{step}/{draws}_permutations/{trait_A}-{trait_B}.tsv"
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
