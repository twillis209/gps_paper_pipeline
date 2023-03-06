from scipy.stats import chi2
import pandas as pd

rule process_ukbb_sum_stats:
    input:
        "resources/ukbb_sum_stats/{snp_set}/merged_ukbb_sum_stats.tsv.gz"
    output:
        "resources/ukbb_sum_stats/{snp_set}/{trait_code}.assoc"
    params:
        z_colname = lambda wildcards: '.'.join(['tstat', wildcards.trait_code]),
        sample_size_colname = lambda wildcards: '.'.join(['n_complete_samples', wildcards.trait_code])
    threads: 12
    resources:
        runtime = 8
    group: "sumher"
    script:
        "../../scripts/process_ukbb_sum_stats.R"

rule estimate_rg_with_ldak_thin_for_ukbb:
    input:
        wg_tagging_file = "results/ldak/ldak-thin/{snp_set}/whole_genome.tagging",
        sum_stats_file_A = "resources/ukbb_sum_stats/{snp_set}/{trait_A}.assoc",
        sum_stats_file_B = "resources/ukbb_sum_stats/{snp_set}/{trait_B}.assoc"
    output:
        cors_full_file = "results/ldak/ldak-thin/{snp_set}/rg/{trait_A}-{trait_B}.cors.full"
    log:
        log_file = "results/ldak/ldak-thin/{snp_set}/rg/{trait_A}-{trait_B}.log"
    params:
        output_stem = "results/ldak/ldak-thin/{snp_set}/rg/{trait_A}-{trait_B}"
    resources:
        runtime = 5
    group: "sumher"
    shell:
        """
        $ldakRoot/ldak --sum-cors {params.output_stem} --tagfile {input.wg_tagging_file} --summary {input.sum_stats_file_A} --summary2 {input.sum_stats_file_B} --allow-ambiguous YES --check-sums NO --cutoff 0.01 > {log.log_file}
        """
