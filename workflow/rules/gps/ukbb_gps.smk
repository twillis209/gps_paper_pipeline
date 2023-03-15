localrules: fit_gev_and_compute_gps_pvalue_for_trait_pair

rule compute_gps_for_trait_pair:
    input:
      "resources/ukbb_sum_stats/{trait_A}.done",
      "resources/ukbb_sum_stats/{trait_B}.done",
      sum_stats_file = "resources/ukbb_sum_stats/{snp_set}/{variant_set}/window_{window}_step_{step}_r2_{r2}/pruned_merged_sum_stats.tsv",
    output:
        temp("results/gps/{snp_set}/{variant_set}/window_{window}_step_{step}_r2_{r2}/{trait_A}-{trait_B}_gps_value.tsv")
    params:
        no_of_pert_iterations = 1
    threads: 10
    group: "gps"
    resources:
        runtime = 10
    shell:
      "workflow/scripts/gps_cpp/build/apps/computeGpsCLI -i {input.sum_stats_file} -a {wildcards.trait_A} -b {wildcards.trait_B} -c {wildcards.trait_A} -d {wildcards.trait_B} -p {params.no_of_perturbations} -n {threads} -f pp -o {output}"

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
        gps_file = "results/gps/{snp_set}/{variant_set}/window_{window}_step_{step}_r2_{r2}/{trait_A}-{trait_B}_gps_value.tsv",
        perm_file = "results/gps/{snp_set}/{variant_set}/window_{window}_step_{step}_r2_{r2}/{draws}_permutations/{trait_A}-{trait_B}.tsv"
    output:
        "results/gps/{snp_set}/{variant_set}/window_{window}_step_{step}_r2_{r2}/{trait_A}-{trait_B}_{draws}_permutations_gps_pvalue.tsv"
    resources:
        runtime = 5
    group: "gps"
    shell:
      "Rscript workflow/scripts/fit_gev_and_compute_gps_pvalue.R -g {input.gps_file} -p {input.perm_file} -a {wildcards.trait_A} -b {wildcards.trait_B} -o {output}"

rule compute_li_gps_pvalue_for_trait_pair:
    input:
        "results/gps/{snp_set}/{variant_set}/window_{window}_step_{step}_r2_{r2}/{trait_A}-{trait_B}_gps_value.tsv"
    output:
        "results/gps/{snp_set}/{variant_set}/window_{window}_step_{step}_r2_{r2}/{trait_A}-{trait_B}_li_gps_pvalue.tsv"
    group: "gps"
    script:
        "../../scripts/compute_li_gps_pvalue.R"
