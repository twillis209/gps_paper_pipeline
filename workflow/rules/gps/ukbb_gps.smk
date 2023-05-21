localrules: fit_gev_and_compute_gps_pvalue_for_trait_pair, compute_li_gps_pvalue_for_trait_pair

rule compute_gps_for_trait_pair:
    input:
      "resources/ukbb_sum_stats/{trait_A}.done",
      "resources/ukbb_sum_stats/{trait_B}.done",
      sum_stats_file = "resources/ukbb_sum_stats/{snp_set}/{variant_set}/window_{window}_step_{step}_r2_{r2}/pruned_merged_sum_stats.tsv",
    output:
        temp("results/gps/{snp_set}/{variant_set}/window_{window}_step_{step}_r2_{r2}/{trait_A}-{trait_B}_gps_value.tsv")
    # Need extra threads for the memory
    threads: 8
    resources:
        runtime = 10,
        mem_mb = get_mem_mb
    group: "gps"
    shell:
      "workflow/scripts/gps_cpp/build/apps/computeGpsCLI -i {input.sum_stats_file} -a {wildcards.trait_A} -b {wildcards.trait_B} -c {wildcards.trait_A} -d {wildcards.trait_B} -n {threads} -f pp -o {output}"

rule permute_trait_pair:
    input:
      "resources/ukbb_sum_stats/{trait_A}.done",
      "resources/ukbb_sum_stats/{trait_B}.done",
      sum_stats_file = "resources/ukbb_sum_stats/{snp_set}/{variant_set}/window_{window}_step_{step}_r2_{r2}/pruned_merged_sum_stats.tsv",
    output:
        "results/gps/{snp_set}/{variant_set}/window_{window}_step_{step}_r2_{r2}/{draws}_permutations/{trait_A}-{trait_B}.tsv"
    threads: 12
    resources:
        mem_mb = get_mem_mb,
        runtime = get_permute_time
    shell:
      "workflow/scripts/gps_cpp/build/apps/permuteTraitsCLI -i {input.sum_stats_file} -o {output} -a {wildcards.trait_A} -b {wildcards.trait_B} -c {threads} -n {wildcards.draws}"

rule fit_gev_and_compute_gps_pvalue_for_trait_pair:
    input:
        gps_file = "results/gps/{snp_set}/{variant_set}/window_{window}_step_{step}_r2_{r2}/{trait_A}-{trait_B}_gps_value.tsv",
        perm_file = "results/gps/{snp_set}/{variant_set}/window_{window}_step_{step}_r2_{r2}/{draws}_permutations/{trait_A}-{trait_B}.tsv"
    output:
        "results/gps/{snp_set}/{variant_set}/window_{window}_step_{step}_r2_{r2}/{trait_A}-{trait_B}_{draws}_permutations_gps_pvalue.tsv"
    params:
        trait_A = lambda wildcards: wildcards.trait_A,
        trait_B = lambda wildcards: wildcards.trait_B
    resources:
        runtime = 1
    script:
      "../../scripts/fit_gev_and_compute_gps_pvalue.R"

rule compute_li_gps_pvalue_for_trait_pair:
    input:
        "results/gps/{snp_set}/{variant_set}/window_{window}_step_{step}_r2_{r2}/{trait_A}-{trait_B}_gps_value.tsv"
    output:
        "results/gps/{snp_set}/{variant_set}/window_{window}_step_{step}_r2_{r2}/{trait_A}-{trait_B}_li_gps_pvalue.tsv"
    resources:
        runtime = 1
    script:
        "../../scripts/compute_li_gps_pvalue.R"
