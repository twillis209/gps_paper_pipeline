rule compute_hoeffdings_for_trait_pair:
    input:
        sum_stats_file = ancient("resources/ukbb_sum_stats/{snp_set}/{variant_set}/window_{window}_step_{step}_r2_{r2}/pruned_merged_sum_stats.tsv"),
    output:
        "results/{snp_set}/{variant_set}/window_{window}_step_{step}_r2_{r2}/{trait_A}-{trait_B}_hoeffdings.tsv"
    script:
        "../scripts/ukbb/compute_hoeffdings.R"

rule collate_hoeffdings_results:
    input:
        ["results/{{snp_set}}/{{variant_set}}/window_{{window}}_step_{{step}}_r2_{{r2}}/{x}_hoeffdings.tsv" for x in ukbb_trait_pairs]
    output:
        "results/combined/{snp_set}/{variant_set}/window_{window}_step_{step}_r2_{r2}/hoeffdings_results.tsv"
    shell:
       """
       echo -e 'trait_A\ttrait_B\tn\tDn\tscaled\tp.value' >> {output}
       for x in {input}; do
          cat <(tail -n 1 $x) >> {output}
       done
       """
