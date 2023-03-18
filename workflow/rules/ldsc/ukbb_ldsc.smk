rule preprocess_ukbb_sum_stats_trait:
    input:
        flag_file = "resources/ukbb_sum_stats/{trait}.done",
        sum_stats_file = "resources/ukbb_sum_stats/{snp_set}/merged_ukbb_sum_stats.tsv.gz",
        snplist_file = "resources/ldsc/ldsc_snps.tsv.gz"
    params:
        pval_col = lambda wildcards: f"pval.{wildcards.trait}",
        tstat_col = lambda wildcards: f"tstat.{wildcards.trait}",
        n_col = lambda wildcards: f"n_complete_samples.{wildcards.trait}",
    output:
        temp("results/ldsc/munged_sum_stats/ukbb/{snp_set}/{trait}_preprocessed_sum_stats.tsv.gz")
    threads: 12
    resources:
        tmpdir = 'tmp'
    group: "ukbb_ldsc"
    script: "../../scripts/ldsc/preprocess_ukbb_sum_stats_for_ldsc.R"

rule munge_ukbb_sum_stats:
    input:
        "results/ldsc/munged_sum_stats/ukbb/{snp_set}/{trait}_preprocessed_sum_stats.tsv.gz"
    output:
        "results/ldsc/munged_sum_stats/ukbb/{snp_set}/{trait}.tsv.sumstats.gz"
    params:
        output_filename = "results/ldsc/munged_sum_stats/ukbb/{snp_set}/{trait}.tsv",
        # NB: The '0' below gives the null value for beta
        signed_sumstats_col = lambda wildcards: f"tstat.{wildcards.trait},0",
        pvalue_col = lambda wildcards: f"pval.{wildcards.trait}",
        n_col = lambda wildcards: f"n_complete_samples.{wildcards.trait}",
        id_col = "SNP"
    log:
        log = "results/ldsc/munged_sum_stats/ukbb/{snp_set}/{trait}.log"
    threads: 1
    resources:
        runtime = 10
    conda:
        "envs/ldsc.yaml"
    group: "ukbb_ldsc"
    shell:
        """
        python $ldsc/munge_sumstats.py --sumstats {input} --N-col {params.n_col} --snp {params.id_col} --out {params.output_filename} --signed-sumstats {params.signed_sumstats_col} --p {params.pvalue_col} --a1 a1 --a2 a2 --frq EUR;
        """

rule estimate_rg_for_ukbb_sum_stats:
    input:
        ["resources/ldsc/eur_w_ld_chr/%d.l2.ldscore.gz" % i for i in range(1,23)],
        ["resources/ldsc/eur_w_ld_chr/%d.l2.M_5_50" % i for i in range(1,23)],
        sum_stats_A = "results/ldsc/munged_sum_stats/ukbb/{snp_set}/{trait_A}.tsv.sumstats.gz",
        sum_stats_B = "results/ldsc/munged_sum_stats/ukbb/{snp_set}/{trait_B}.tsv.sumstats.gz"
    output:
        "results/ldsc/rg/ukbb/{snp_set}/{h2_intercept,fixed|free}_h2_{rg_intercept,fixed|free}_rg_intercept/{trait_A}-{trait_B}.log"
    log:
        log = "results/ldsc/rg/ukbb/{snp_set}/{h2_intercept,fixed|free}_h2_{rg_intercept,fixed|free}_rg_intercept/{trait_A}-{trait_B}.actual_log"
    params:
        log_file_par = "results/ldsc/rg/ukbb/{snp_set}/{h2_intercept,fixed|free}_h2_{rg_intercept,fixed|free}_rg_intercept/{trait_A}-{trait_B}",
        # NB: Trailing '/' is needed in ld_score_root
        ld_score_root = "resources/ldsc/eur_w_ld_chr/",
        h2_intercept = lambda wildcards: "--intercept-h2 1,1" if wildcards.h2_intercept == "fixed" else "",
        rg_intercept = lambda wildcards: "--intercept-gencov 0,0" if wildcards.rg_intercept == "fixed" else ""
    resources:
        runtime = 2
    group: "ukbb_ldsc"
    conda:
        "envs/ldsc.yaml"
    shell:
        # Hacky fix to retain 'log' file (the output is regrettably so-called), which we need whether or not the estimation process failed
        "python $ldsc/ldsc.py --rg {input.sum_stats_A},{input.sum_stats_B} --ref-ld-chr {params.ld_score_root} --w-ld-chr {params.ld_score_root} --out {params.log_file_par} {params.h2_intercept} {params.rg_intercept} >{log.log} || true"
