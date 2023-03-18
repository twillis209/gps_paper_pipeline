rule compute_hoeffdings_for_trait_pair:
    input:
        sum_stats_file = ancient("resources/ukbb_sum_stats/{snp_set}/{variant_set}/window_{window}_step_{step}_r2_{r2}/pruned_merged_sum_stats.tsv"),
    output:
        "results/{snp_set}/{variant_set}/window_{window}_step_{step}_r2_{r2}/{trait_A}-{trait_B}_hoeffdings.tsv"
    group: "gps"
    script:
        "../../scripts/ukbb/compute_hoeffdings.R"
