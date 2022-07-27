rule munge_randomised_sum_stats:
    input:
        "results/simgwas/simulated_sum_stats/whole_genome_sum_stats/randomised/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/{effect_blocks}_seed_{seed}_sum_stats_{pair_label}_tag_{tag}_of_{tags}.tsv.gz"
    output:
        temp("results/ldsc/munged_sum_stats/whole_genome_sum_stats/randomised/{ncases_A,\d+}_{ncontrols_A,\d+}_{ncases_B,\d+}_{ncontrols_B,\d+}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/{effect_blocks}_seed_{seed,\d+}_{pair_label,[AB]}_{tag,\w}_of_{tags,\w+}.tsv.sumstats.gz")
    params:
        output_filename = "results/ldsc/munged_sum_stats/whole_genome_sum_stats/randomised/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/{effect_blocks}_seed_{seed}_{pair_label}_{tag}_of_{tags}.tsv",
         signed_sumstats_col = lambda wildcards: tag_signed_sumstats_dict[wildcards.tag],
         pvalue_col = lambda wildcards: tag_pvalue_dict[wildcards.tag]
    threads: 2
    resources:
        runtime = 10
    group: "ldsc_hoeffding_and_gps_sans_permutation"
    conda:
        "envs/ldsc.yaml"
    shell:
        """
        python $ldsc/munge_sumstats.py --sumstats {input} --N-con-col ncontrols --N-cas-col ncases --snp id --out {params.output_filename} --signed-sumstats {params.signed_sumstats_col} --p {params.pvalue_col} --a1 a0 --a2 a1 --frq EUR --a1-inc;
        """

rule estimate_rg_for_randomised_sum_stats:
    input:
        ["resources/ldsc/eur_w_ld_chr/%d.l2.ldscore.gz" % i for i in range(1,23)],
        ["resources/ldsc/eur_w_ld_chr/%d.l2.M_5_50" % i for i in range(1,23)],
        sum_stats_A = "results/ldsc/munged_sum_stats/whole_genome_sum_stats/randomised/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/{effect_blocks_A}_seed_{seed}_A_{tag_A}_of_{tag_A}{tag_B}.tsv.sumstats.gz",
        sum_stats_B = "results/ldsc/munged_sum_stats/whole_genome_sum_stats/randomised/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/{effect_blocks_B}_seed_{seed}_B_{tag_B}_of_{tag_A}{tag_B}.tsv.sumstats.gz"
    output:
        "results/ldsc/rg/{h2_intercept,fixed|free}_h2_{rg_intercept,fixed|free}_rg_intercept/whole_genome/randomised/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}_seed_{seed}_{tag_A}{tag_B}.log"
    log:
        log = "results/ldsc/rg/{h2_intercept}_h2_{rg_intercept}_rg_intercept/whole_genome/randomised/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}_seed_{seed}_{tag_A}{tag_B}.actual_log"
    params:
        log_file_par = "results/ldsc/rg/{h2_intercept}_h2_{rg_intercept}_rg_intercept/whole_genome/randomised/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}_seed_{seed}_{tag_A}{tag_B}",
        # NB: Trailing '/' is needed in ld_score_root
        ld_score_root = "resources/ldsc/eur_w_ld_chr/",
        population_prevalence = 0.02,
        sample_prevalence = 0.5,
        h2_intercept = lambda wildcards: "--intercept-h2 1,1" if wildcards.h2_intercept == "fixed" else "",
        rg_intercept = lambda wildcards: "--intercept-gencov 0,0" if wildcards.rg_intercept == "fixed" else ""
    resources:
        runtime = 2
    group: "ldsc_hoeffding_and_gps_sans_permutation"
    conda:
        "envs/ldsc.yaml"
    shell:
        # Hacky fix to retain 'log' file (the output is regrettably so-called), which we need whether or not the estimation process failed
        "python $ldsc/ldsc.py --rg {input.sum_stats_A},{input.sum_stats_B} --ref-ld-chr {params.ld_score_root} --w-ld-chr {params.ld_score_root} --out {params.log_file_par} --samp-prev {params.sample_prevalence},{params.sample_prevalence} --pop-prev {params.population_prevalence},{params.population_prevalence} {params.h2_intercept} {params.rg_intercept} >{log.log} || true"

rule write_out_randomised_blocks_for_pair:
    output:
        a_chrom_blocks_file = "results/ldsc/rg/whole_genome/randomised/theoretical_rg/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/block_files/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}_seed_{seed}_{tag_A}_{tag_A}{tag_B}.tsv",
        b_chrom_blocks_file = "results/ldsc/rg/whole_genome/randomised/theoretical_rg/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/block_files/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}_seed_{seed}_{tag_B}_{tag_A}{tag_B}.tsv",
        shared_chrom_blocks_file = "results/ldsc/rg/whole_genome/randomised/theoretical_rg/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/block_files/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}_seed_{seed}_shared_{tag_A}{tag_B}.tsv"
    params:
        chrom_block_tuples = lambda wildcards: get_randomised_chrom_block_tuples_for_pair(wildcards),
    resources:
        runtime = 2
    group: "calculate_theoretical_rg"
    run:
        _, _, _, shared_chrom_block_dict, a_chrom_block_dict, b_chrom_block_dict = params.chrom_block_tuples

        block_a_daf = chrom_block_dict_to_dataframe(a_chrom_block_dict)
        block_a_daf.sort_values(by = ['chr', 'block', 'effect'], inplace = True)

        block_b_daf = chrom_block_dict_to_dataframe(b_chrom_block_dict)
        block_b_daf.sort_values(by = ['chr', 'block', 'effect'], inplace = True)

        shared_block_daf = chrom_block_dict_to_dataframe(shared_chrom_block_dict)
        shared_block_daf.sort_values(by = ['chr', 'block', 'effect'], inplace = True)

        pd.concat([block_a_daf, shared_block_daf]).to_csv(output.a_chrom_blocks_file, index = False, sep = '\t')

        pd.concat([block_b_daf, shared_block_daf]).to_csv(output.b_chrom_blocks_file, index = False, sep = '\t')

        shared_block_daf.to_csv(output.shared_chrom_blocks_file, index = False, sep = '\t')

rule calculate_theoretical_rg_for_randomised_sum_stats:
    input:
        combined_causal_variants_file = "results/simgwas/combined_causal_variants.tsv",
        a_chrom_blocks_file = "results/ldsc/rg/whole_genome/randomised/theoretical_rg/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/block_files/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}_seed_{seed}_{tag_A}_{tag_A}{tag_B}.tsv",
        b_chrom_blocks_file = "results/ldsc/rg/whole_genome/randomised/theoretical_rg/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/block_files/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}_seed_{seed}_{tag_B}_{tag_A}{tag_B}.tsv"
    output:
        theo_rg_file = "results/ldsc/rg/whole_genome/randomised/theoretical_rg/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A,[smlh0-9\-]+}_{effect_blocks_B,[smlh0-9\-]+}_{shared_effect_blocks,[smlh0-9\-]+}_seed_{seed,\d+}_{tag_A,\w}{tag_B,\w}_theo_rg.tsv",
    params:
        population_prevalence_A = 0.02,
        sample_prevalence_A = 0.5,
        population_prevalence_B = 0.02,
        sample_prevalence_B = 0.5,
    resources:
        runtime = 2
    group: "calculate_theoretical_rg"
    shell:
        "Rscript workflow/scripts/ldsc/calculate_theoretical_rg_randomised_blocks.R --cv_file {input.combined_causal_variants_file} --a_blocks_file {input.a_chrom_blocks_file} --b_blocks_file {input.b_chrom_blocks_file} --P_a {params.sample_prevalence_A} --P_b {params.sample_prevalence_B} --K_a {params.population_prevalence_A} --K_b {params.population_prevalence_B} -o {output.theo_rg_file} -nt {threads}"
