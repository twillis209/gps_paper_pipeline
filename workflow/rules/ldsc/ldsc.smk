from string import ascii_lowercase

tags = list(ascii_lowercase[:20])
tag_pvalue_dict = dict(zip(tags, [f"p.{x}" for x in range(1,21)]))
tag_signed_sumstats_dict = dict(zip(tags, [f"betasim.{x},{x+17}" for x in range(1,21)]))

rule download_ld_scores:
    output:
        "resources/ldsc/eur_w_ld_chr.tar.bz2"
    shell:
        "wget -O {output} https://storage.googleapis.com/broad-alkesgroup-public/LDSCORE/eur_w_ld_chr.tar.bz2"

rule extract_ld_scores:
    input:
        "resources/ldsc/eur_w_ld_chr.tar.bz2"
    output:
        ["resources/ldsc/eur_w_ld_chr/%d.l2.ldscore.gz" % i for i in range(1,23)],
        ["resources/ldsc/eur_w_ld_chr/%d.l2.M_5_50" % i for i in range(1,23)],
        "resources/ldsc/eur_w_ld_chr/w_hm3.snplist",
        "resources/ldsc/eur_w_ld_chr/README"
    params:
        output_root = "resources/ldsc"
    shell:
        "tar -xjf {input} -C {params.output_root}"

rule write_out_ld_score_snps:
    input:
        ldsc_files = ["resources/ldsc/eur_w_ld_chr/%d.l2.ldscore.gz" % i for i in range(1,23)],
        snplist_file = "resources/ldsc/eur_w_ld_chr/w_hm3.snplist",
    output:
        "resources/ldsc/ldsc_snps.tsv.gz"
    threads: 4
    script: "../../scripts/ldsc/write_out_ld_score_snps.R"

# TODO Fix like I fixed theoretical rg
rule calculate_theoretical_h2:
    input:
        "results/simgwas/combined_causal_variants.tsv"
    output:
        "results/ldsc/h2/whole_genome/{effect_blocks}_theo_h2.tsv"
    params:
        population_prevalence = 0.02,
        sample_prevalence = 0.5,
        odds_ratio = lambda wildcards: odds_ratio_dict[re.search("[smlvhin]", wildcards.effect_blocks).group()],
    threads: 4
    shell:
        "Rscript workflow/scripts/ldsc/calculate_theoretical_h2.R --cv_file {input} --effect_blocks {wildcards.effect_blocks} --odds_ratio {params.odds_ratio} -P {params.sample_prevalence} -K {params.population_prevalence} -o {output} -nt {threads}"

rule preprocess_ukbb_sum_stats_trait:
    input:
        flag_file = "resources/ukbb_sum_stats/{trait}.done",
        sum_stats_file = "resources/ukbb_sum_stats/{join}/merged_ukbb_sum_stats.tsv.gz",
        snplist_file = "resources/ldsc/ldsc_snps.tsv.gz"
    params:
        pval_col = lambda wildcards: f"pval.{wildcards.trait}",
        tstat_col = lambda wildcards: f"tstat.{wildcards.trait}",
        n_col = lambda wildcards: f"n_complete_samples.{wildcards.trait}",
    output:
        temp("results/ldsc/munged_sum_stats/ukbb/{join}/{trait}_preprocessed_sum_stats.tsv.gz")
    threads: 4
    resources:
        tmpdir = 'tmp'
    script: "../../scripts/ldsc/preprocess_ukbb_sum_stats_for_ldsc.R"

rule munge_ukbb_sum_stats:
    input:
        "results/ldsc/munged_sum_stats/ukbb/{join}/{trait}_preprocessed_sum_stats.tsv.gz"
    output:
        "results/ldsc/munged_sum_stats/ukbb/{join}/{trait}.tsv.sumstats.gz"
    params:
        output_filename = "results/ldsc/munged_sum_stats/ukbb/{join}/{trait}.tsv",
        # NB: The '0' below gives the null value for beta
        signed_sumstats_col = lambda wildcards: f"tstat.{wildcards.trait},0",
        pvalue_col = lambda wildcards: f"pval.{wildcards.trait}",
        n_col = lambda wildcards: f"n_complete_samples.{wildcards.trait}",
        id_col = "SNP"
    log:
        log = "results/ldsc/munged_sum_stats/ukbb/{join}/{trait}.log"
    threads: 1
    resources:
        runtime = 10
    conda:
        "envs/ldsc.yaml"
    shell:
        """
        python $ldsc/munge_sumstats.py --sumstats {input} --N-col {params.n_col} --snp {params.id_col} --out {params.output_filename} --signed-sumstats {params.signed_sumstats_col} --p {params.pvalue_col} --a1 a1 --a2 a2 --frq EUR;
        """

rule estimate_rg_for_ukbb_sum_stats:
    input:
        ["resources/ldsc/eur_w_ld_chr/%d.l2.ldscore.gz" % i for i in range(1,23)],
        ["resources/ldsc/eur_w_ld_chr/%d.l2.M_5_50" % i for i in range(1,23)],
        sum_stats_A = "results/ldsc/munged_sum_stats/ukbb/{join}/{trait_A}.tsv.sumstats.gz",
        sum_stats_B = "results/ldsc/munged_sum_stats/ukbb/{join}/{trait_B}.tsv.sumstats.gz"
    output:
        "results/ldsc/rg/ukbb/{join}/{h2_intercept,fixed|free}_h2_{rg_intercept,fixed|free}_rg_intercept/{trait_A}-{trait_B}.log"
    log:
        log = "results/ldsc/rg/ukbb/{join}/{h2_intercept,fixed|free}_h2_{rg_intercept,fixed|free}_rg_intercept/{trait_A}-{trait_B}.actual_log"
    params:
        log_file_par = "results/ldsc/rg/ukbb/{join}/{h2_intercept,fixed|free}_h2_{rg_intercept,fixed|free}_rg_intercept/{trait_A}-{trait_B}",
        # NB: Trailing '/' is needed in ld_score_root
        ld_score_root = "resources/ldsc/eur_w_ld_chr/",
        h2_intercept = lambda wildcards: "--intercept-h2 1,1" if wildcards.h2_intercept == "fixed" else "",
        rg_intercept = lambda wildcards: "--intercept-gencov 0,0" if wildcards.rg_intercept == "fixed" else ""
    resources:
        runtime = 2
    group: "ldsc_hoeffding_sumher_gps_sans_permutation"
    conda:
        "envs/ldsc.yaml"
    shell:
        # Hacky fix to retain 'log' file (the output is regrettably so-called), which we need whether or not the estimation process failed
        "python $ldsc/ldsc.py --rg {input.sum_stats_A},{input.sum_stats_B} --ref-ld-chr {params.ld_score_root} --w-ld-chr {params.ld_score_root} --out {params.log_file_par} {params.h2_intercept} {params.rg_intercept} >{log.log} || true"

        # TODO some of the variants specified seems to be missing from the cv_dat, 1_100, 5_108, 6_6, 6_25, 6_38, 8_22
#rule calculate_theoretical_rg:
#    input:
#        cv_file = "results/simgwas/combined_causal_variants.tsv",
#        available_blocks_file = "resources/simgwas/available_blocks.tsv"
#    output:
#        "results/ldsc/rg/whole_genome/{effect_blocks_A}_{effect_blocks_B}_theo_rg.tsv"
#    params:
#        population_prevalence_A = 0.02,
#        sample_prevalence_A = 0.5,
#        population_prevalence_B = 0.02,
#        sample_prevalence_B = 0.5,
#        # TODO handle null and non-uniform effect cases
#        odds_ratio_a = lambda wildcards: odds_ratio_dict[re.search("[smlvhin]", wildcards.effect_blocks_A).group()],
#        odds_ratio_b = lambda wildcards: odds_ratio_dict[re.search("[smlvhin]", wildcards.effect_blocks_B).group()]
#    threads: 4
#    shell:
#        "Rscript workflow/scripts/ldsc/calculate_theoretical_rg.R --cv_file {input.cv_file} --blocks_file {input.available_blocks_file} --effect_blocks_a {wildcards.effect_blocks_A} --effect_blocks_b {wildcards.effect_blocks_B} --odds_ratio_a {params.odds_ratio_a} --odds_ratio_b {params.odds_ratio_b} --P_a {params.sample_prevalence_A} --P_b {params.sample_prevalence_B} --K_a {params.population_prevalence_A} --K_b {params.population_prevalence_B} -o {output} -nt {threads}"
