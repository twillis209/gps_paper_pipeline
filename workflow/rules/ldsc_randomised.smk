rule munge_randomised_sum_stats:
    input:
        "results/simgwas/simulated_sum_stats/whole_genome_sum_stats/randomised/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/{effect_blocks}_seed_{seed}_sum_stats_{pair_label}_tag_{tag}_of_{tags}.tsv.gz"
    output:
        temp("results/ldsc/munged_sum_stats/whole_genome_sum_stats/randomised/{ncases_A,\d+}_{ncontrols_A,\d+}_{ncases_B,\d+}_{ncontrols_B,\d+}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/{effect_blocks}_seed_{seed,\d+}_{pair_label,[AB]}_{tag,\w}_of_{tags,\w+}.tsv.sumstats.gz")
    params:
        output_filename = "results/ldsc/munged_sum_stats/whole_genome_sum_stats/randomised/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/{effect_blocks}_seed_{seed}_{pair_label}_{tag}_of_{tags}.tsv",
         signed_sumstats_col = lambda wildcards: tag_signed_sumstats_dict[wildcards.tag],
         pvalue_col = lambda wildcards: tag_pvalue_dict[wildcards.tag]
    conda:
        "envs/ldsc.yaml"
    shell:
        """
        python $ldsc/munge_sumstats.py --sumstats {input} --N-con-col ncontrols --N-cas-col ncases --snp id --out {params.output_filename} --signed-sumstats {params.signed_sumstats_col} --p {params.pvalue_col} --a1 a0 --a2 a1 --frq EUR --a1-inc;
        """

# TODO revise simgwas_randomised to add in single-trait scheme
rule estimate_randomised_h2_B0_1:
    input:
        ["resources/ldsc/eur_w_ld_chr/%d.l2.ldscore.gz" % i for i in range(1,23)],
        ["resources/ldsc/eur_w_ld_chr/%d.l2.M_5_50" % i for i in range(1,23)],
        sum_stats = "results/ldsc/munged_sum_stats/whole_genome_sum_stats/{ncases}_{ncontrols}/{effect_blocks}_{tag}.tsv.sumstats.gz"
    output:
        "results/ldsc/h2/whole_genome/{ncases,\d+}_{ncontrols,\d+}/B0_1/{effect_blocks,[^/]+}_{tag,[abcde]}.log"
    params:
        log_file_par = "results/ldsc/h2/whole_genome/{ncases}_{ncontrols}/B0_1/{effect_blocks}_{tag,[abcde]}",
        # NB: Trailing '/' is needed in ld_score_root
        ld_score_root = "resources/ldsc/eur_w_ld_chr/",
        population_prevalence = 0.02,
        sample_prevalence = 0.5,
        intercept = 1
    resources:
        disk_mb = 2048
    conda:
        "envs/ldsc.yaml"
    shell:
        "python $ldsc/ldsc.py --h2 {input.sum_stats} --ref-ld-chr {params.ld_score_root} --w-ld-chr {params.ld_score_root} --out {params.log_file_par} --samp-prev {params.sample_prevalence} --pop-prev {params.population_prevalence} --intercept-h2 {params.intercept}"

rule estimate_rg_for_randomised_sum_stats:
    input:
        ["resources/ldsc/eur_w_ld_chr/%d.l2.ldscore.gz" % i for i in range(1,23)],
        ["resources/ldsc/eur_w_ld_chr/%d.l2.M_5_50" % i for i in range(1,23)],
        sum_stats_A = "results/ldsc/munged_sum_stats/whole_genome_sum_stats/randomised/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/{effect_blocks_A}_seed_{seed}_A_{tag_A}_of_{tag_A}{tag_B}.tsv.sumstats.gz",
        sum_stats_B = "results/ldsc/munged_sum_stats/whole_genome_sum_stats/randomised/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/{effect_blocks_B}_seed_{seed}_B_{tag_B}_of_{tag_A}{tag_B}.tsv.sumstats.gz"
    output:
        "results/ldsc/rg/whole_genome/randomised/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}_seed_{seed}_{tag_A}{tag_B}.log"
    params:
        log_file_par = "results/ldsc/rg/whole_genome/randomised/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}_seed_{seed}_{tag_A}{tag_B}",
        # NB: Trailing '/' is needed in ld_score_root
        ld_score_root = "resources/ldsc/eur_w_ld_chr/",
        population_prevalence = 0.02,
        sample_prevalence = 0.5
    conda:
        "envs/ldsc.yaml"
    shell:
        "python $ldsc/ldsc.py --rg {input.sum_stats_A},{input.sum_stats_B} --ref-ld-chr {params.ld_score_root} --w-ld-chr {params.ld_score_root} --out {params.log_file_par} --samp-prev {params.sample_prevalence},{params.sample_prevalence} --pop-prev {params.population_prevalence},{params.population_prevalence} --intercept-h2 1,1"

rule write_out_randomised_blocks_for_pair:
    output:
        a_chrom_blocks_file = temp("results/ldsc/rg/whole_genome/randomised/theoretical_rg/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/block_files/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}_seed_{seed}_{tag_A}_{tag_A}{tag_B}.tsv"),
        b_chrom_blocks_file = temp("results/ldsc/rg/whole_genome/randomised/theoretical_rg/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/block_files/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}_seed_{seed}_{tag_B}_{tag_A}{tag_B}.tsv"),
        shared_chrom_blocks_file = temp("results/ldsc/rg/whole_genome/randomised/theoretical_rg/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/block_files/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}_seed_{seed}_shared_{tag_A}{tag_B}.tsv")
    params:
        chrom_block_tuples = lambda wildcards: get_randomised_chrom_block_tuples_for_pair(wildcards),
    run:
        block_a_daf = pd.DataFrame(params.chrom_block_tuples[0]+params.chrom_block_tuples[1], columns = ['chr', 'block'])
        block_a_daf.sort_values(by = ['chr', 'block'], inplace = True)
        block_a_daf.to_csv(output.a_chrom_blocks_file, index = False, sep = '\t')

        block_b_daf = pd.DataFrame(params.chrom_block_tuples[0]+params.chrom_block_tuples[2], columns = ['chr', 'block'])
        block_b_daf.sort_values(by = ['chr', 'block'], inplace = True)
        block_b_daf.to_csv(output.b_chrom_blocks_file, index = False, sep = '\t')

        shared_block_daf = pd.DataFrame(params.chrom_block_tuples[0], columns = ['chr', 'block'])
        shared_block_daf.sort_values(by = ['chr', 'block'], inplace = True)
        shared_block_daf.to_csv(output.shared_chrom_blocks_file, index = False, sep = '\t')

rule calculate_theoretical_rg_for_randomised_sum_stats:
    input:
        combined_causal_variants_file = "results/simgwas/combined_causal_variants.tsv",
        a_chrom_blocks_file = "results/ldsc/rg/whole_genome/randomised/theoretical_rg/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/block_files/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}_seed_{seed}_{tag_A}_{tag_A}{tag_B}.tsv",
        b_chrom_blocks_file = "results/ldsc/rg/whole_genome/randomised/theoretical_rg/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/block_files/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}_seed_{seed}_{tag_B}_{tag_A}{tag_B}.tsv"
    output:
        theo_rg_file = "results/ldsc/rg/whole_genome/randomised/theoretical_rg/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A,[smlhv]\d+}_{effect_blocks_B,[smlhv]\d+}_{shared_effect_blocks,[smlhv]\d+}_seed_{seed,\d+}_{tag_A}{tag_B}_theo_rg.tsv",
    params:
        population_prevalence_A = 0.02,
        sample_prevalence_A = 0.5,
        population_prevalence_B = 0.02,
        sample_prevalence_B = 0.5,
        odds_ratio_a = lambda wildcards: odds_ratio_dict[re.search("[smlvhin]", wildcards.effect_blocks_A).group()],
        odds_ratio_b = lambda wildcards: odds_ratio_dict[re.search("[smlvhin]", wildcards.effect_blocks_B).group()]
    threads: 1
    shell:
        "Rscript workflow/scripts/ldsc/calculate_theoretical_rg_randomised_blocks.R --cv_file {input.combined_causal_variants_file} --a_blocks_file {input.a_chrom_blocks_file} --b_blocks_file {input.b_chrom_blocks_file} --odds_ratio_a {params.odds_ratio_a} --odds_ratio_b {params.odds_ratio_b} --P_a {params.sample_prevalence_A} --P_b {params.sample_prevalence_B} --K_a {params.population_prevalence_A} --K_b {params.population_prevalence_B} -o {output.theo_rg_file} -nt {threads}"
