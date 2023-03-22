rule download_ukbb_sum_stats:
    output:
        "resources/ukbb_sum_stats/{id}.gwas.imputed_v3.both_sexes.tsv.bgz"
    resources:
        runtime = 5
    group: 'ukbb'
    shell:
        """
        wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/{wildcards.id}.gwas.imputed_v3.both_sexes.tsv.bgz -O resources/ukbb_sum_stats/{wildcards.id}.gwas.imputed_v3.both_sexes.tsv.bgz
        """

rule decompress_ukbb_sum_stats:
    input:
        "resources/ukbb_sum_stats/{id}.gwas.imputed_v3.both_sexes.tsv.bgz"
    output:
        decomp_file = temp("resources/ukbb_sum_stats/{id}.gwas.imputed_v3.both_sexes.tsv"),
        flag_file = "resources/ukbb_sum_stats/{id}.done"
    resources:
        runtime = 5
    group: 'ukbb'
    shell:
        """
        gunzip -c {input} >{output.decomp_file} && touch {output.flag_file}
        """

rule merge_ukbb_sum_stats:
    input:
        ukbb_files = ["resources/ukbb_sum_stats/%s.gwas.imputed_v3.both_sexes.tsv" % x for x in ukbb_trait_codes]
    output:
        merged_file = "resources/ukbb_sum_stats/{snp_set}/all/merged_ukbb_sum_stats.tsv.gz"
    params:
        ukbb_trait_codes = ukbb_trait_codes,
        sans_mhc = lambda wildcards: True if wildcards.snp_set == 'sans_mhc' else False
    threads: 12
    resources:
        runtime = 20
    group: 'ukbb'
    script: "../../scripts/ukbb/merge_ukbb_sum_stats.R"

rule prune_merged_sum_stats:
    input:
        sum_stats_file = "resources/ukbb_sum_stats/{snp_set}/all/merged_ukbb_sum_stats.tsv.gz",
        bim_file = "resources/1000g/euro/qc/{snp_set}/{variant_set}/all.bim",
        pruned_range_file = "resources/1000g/euro/qc/{snp_set}/{variant_set}/pruned_ranges/window_{window}_step_{step}_r2_{r2}/all.prune.in"
    output:
        "resources/ukbb_sum_stats/{snp_set}/{variant_set}/window_{window}_step_{step}_r2_{r2}/pruned_merged_sum_stats.tsv"
    threads: 12
    resources:
        runtime = 30
    group: 'ukbb'
    script: "../../scripts/ukbb/prune_merged_sum_stats.R"

rule downsample_pruned_merged_sum_stats:
    input:
        "resources/ukbb_sum_stats/{snp_set}/{variant_set}/window_{window}_step_{step}_r2_{r2}/pruned_merged_sum_stats.tsv"
    output:
        temp("resources/ukbb_sum_stats/{snp_set}/{variant_set}/window_{window}_step_{step}_r2_{r2}/{no_snps}_snps/pruned_merged_sum_stats.tsv")
    threads: 12
    group: 'ukbb'
    script:
     "../../scripts/ukbb/downsample_sum_stats.R"

rule count_lines_in_ukbb_sum_stats:
    input:
        "resources/ukbb_sum_stats/{snp_set}/merged_ukbb_sum_stats.tsv.gz"
    output:
        "resources/ukbb_sum_stats/{snp_set}/merged_ukbb_sum_stats_linecount.txt"
    shell:
        "tail -n +2 {input} | wc -l >{output}"

rule count_lines_in_ukbb_pruned_sum_stats:
    input:
        "results/ukbb_sum_stats/{snp_set}/{variant_set}/window_{window}_step_{step}_r2_{r2}/pruned_merged_sum_stats.tsv"
    output:
        "results/ukbb_sum_stats/{snp_set}/{variant_set}/window_{window}_step_{step}_r2_{r2}/pruned_merged_sum_stats_linecount.txt"
    shell:
        "tail -n +2 {input} | wc -l >{output}"
