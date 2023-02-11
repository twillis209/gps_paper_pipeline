rule compute_gps_for_trait_pair:
    input:
      "resources/ukbb_sum_stats/{trait_A}.done",
      "resources/ukbb_sum_stats/{trait_B}.done",
      sum_stats_file = "resources/ukbb_sum_stats/{snp_set}/{variant_set}/window_{window}_step_{step}_r2_{r2}/pruned_merged_sum_stats.tsv",
    output:
        temp("results/gps/{snp_set}/{variant_set}/window_{window}_step_{step}_r2_{r2}/{trait_A}-{trait_B}_{draws}_permutations_gps_value.tsv")
    log:
        "results/gps/{snp_set}/{variant_set}//window_{window}_step_{step}_r2_{r2}/{trait_A}-{trait_B}_{draws}_permutations_gps_value.log"
    params:
        no_of_pert_iterations = 0
    threads: 10
    group: "gps"
    resources:
        runtime = 10
    shell:
      "workflow/scripts/gps_cpp/build/apps/computeGpsCLI -i {input.sum_stats_file} -a {wildcards.trait_A} -b {wildcards.trait_B} -c {wildcards.trait_A} -d {wildcards.trait_B} -p {params.no_of_perturbations} -n {threads} -o {output}"

rule permute_trait_pair:
    input:
      "resources/ukbb_sum_stats/{trait_A}.done",
      "resources/ukbb_sum_stats/{trait_B}.done",
      sum_stats_file = "resources/ukbb_sum_stats/{snp_set}/{variant_set}/window_{window}_step_{step}_r2_{r2}/pruned_merged_sum_stats.tsv",
    output:
        "results/gps/{snp_set}/{variant_set}/window_{window}_step_{step}_r2_{r2}/{draws}_permutations/{trait_A}-{trait_B}.tsv"
    params:
        no_of_perturbations = 1
    threads: 12
    resources:
        mem_mb = get_mem_mb,
        runtime = get_permute_time,
    group: "gps"
    shell:
      "workflow/scripts/gps_cpp/build/apps/permuteTraitsCLI -i {input.sum_stats_file} -o {output} -a {wildcards.trait_A} -b {wildcards.trait_B} -c {threads} -n {wildcards.draws} -p {params.no_of_perturbations}"

rule fit_gev_and_compute_gps_pvalue_for_trait_pair:
    input:
        gps_file = "results/gps/{snp_set}/{variant_set}/window_{window}_step_{step}_r2_{r2}/{trait_A}-{trait_B}_{draws}_permutations_gps_value.tsv",
        perm_file = "results/gps/{snp_set}/{variant_set}/window_{window}_step_{step}_r2_{r2}/{draws}_permutations/{trait_A}-{trait_B}.tsv"
    output:
        "results/gps/{snp_set}/{variant_set}/window_{window}_step_{step}_r2_{r2}/{trait_A}-{trait_B}_{draws}_permutations_gps_pvalue.tsv"
    resources:
        runtime = 5
    group: "gps"
    shell:
      "Rscript workflow/scripts/fit_gev_and_compute_gps_pvalue.R -g {input.gps_file} -p {input.perm_file} -a {wildcards.trait_A} -b {wildcards.trait_B} -o {output}"

rule collate_gps_pvalue_data:
    input:
        pvalue_files = ["results/gps/{snp_set}/{variant_set}/window_{window}_step_{step}_r2_{r2}/%s_{draws}_permutations_gps_pvalue.tsv" % x for x in ukbb_trait_pairs]
    output:
        combined_pvalue_file = "results/gps/combined/{snp_set}/{variant_set}/window_{window}_step_{step}_r2_{r2}/gps_pvalues_{draws}_permutations.tsv"
    resources:
        runtime = 5
    run:
        with open(output.combined_pvalue_file, 'w') as outfile:
            outfile.write(("\t".join(["trait_A", "trait_B", "gps", "n", "loc", "loc.sd", "scale", "scale.sd", "shape", "shape.sd", "pval"]))+"\n")
            for i,x in enumerate(input.pvalue_files):
                with open(x, 'r') as infile:
                    line = infile.readline()
                    line = infile.readline()

                    m = re.match("results/gps/ukbb/%s/window_%s_step_%s_r2_%s/(\w+)-(\w+)_%s_permutations_gps_pvalue.tsv" % (wildcards.snp_set, wildcards.window, wildcards.step, wildcards.r2, wildcards.draws), x)

                outfile.write(("\t".join([m[1], m[2], line])))
