rule compute_hoeffdings_for_trait_pair:
    input:
        sum_stats_file = ancient("resources/pruned_sum_stats/{join}/{snp_set}/window_{window}_step_{step}/pruned_merged_sum_stats.tsv"),
    output:
        "results/{join}/{snp_set}/window_{window}_step_{step}/{trait_A}-{trait_B}_hoeffdings.tsv"
    shell:
        "Rscript workflow/scripts/compute_hoeffdings.R -i {input.sum_stats_file} -a {wildcards.trait_A} -b {wildcards.trait_B} -o {output} -nt 1"

rule collate_hoeffdings_results:
    input:
      ["results/ukbb/{snp_set}/window_{window}_step_{step}/%s_hoeffdings.tsv" % x for x in ukbb_trait_pairs]
    output:
      "results/combined/{snp_set}/window_{window}_step_{step}/hoeffdings_results.tsv"
    shell:
       """
       echo -e 'trait_A\ttrait_B\tn\tDn\tscaled\tp.value' >> {output}
       for x in {input}; do
          cat <(tail -n 1 $x) >> {output}
       done
       """
