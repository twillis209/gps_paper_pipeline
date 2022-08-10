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
        "results/ldsc/munged_sum_stats/ukbb/{join}/{trait}_preprocessed_sum_stats.tsv.gz",
    threads: 4
    resources:
        tmpdir = 'tmp'
    script: "../../scripts/ldsc/preprocess_ukbb_sum_stats_for_ldsc.R"

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
