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

rule munge_randomised_sum_stats:
    input:
        "results/simgwas/simulated_sum_stats/whole_genome_sum_stats/randomised/seed_{seed}/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/{effect_blocks}_sum_stats_{pair_label}.tsv.gz"
    output:
        temp("results/ldsc/munged_sum_stats/whole_genome_sum_stats/randomised/seed_{seed}/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/{effect_blocks}_{pair_label,[AB]}_{tag}.tsv.sumstats.gz")
    params:
        output_filename = "results/ldsc/munged_sum_stats/whole_genome_sum_stats/randomised/seed_{seed}/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/{effect_blocks}_sum_stats_{pair_label}_{tag}.tsv",
         signed_sumstats_col = lambda wildcards: tag_signed_sumstats_dict[wildcards.tag],
         pvalue_col = lambda wildcards: tag_pvalue_dict[wildcards.tag]
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
        sum_stats_A = "results/ldsc/munged_sum_stats/whole_genome_sum_stats/randomised/seed_{seed}/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/{effect_blocks_A}_A_{tag_A}.tsv.sumstats.gz",
        sum_stats_B = "results/ldsc/munged_sum_stats/whole_genome_sum_stats/randomised/seed_{seed}/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/{effect_blocks_B}_B_{tag_B}.tsv.sumstats.gz"
    output:
        "results/ldsc/rg/whole_genome/randomised/seed_{seed}/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}_{tag_A}{tag_B}.log"
    params:
        log_file_par = "results/ldsc/rg/whole_genome/randomised/seed_{seed}/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}_{tag_A}{tag_B}",
        # NB: Trailing '/' is needed in ld_score_root
        ld_score_root = "resources/ldsc/eur_w_ld_chr/",
        population_prevalence = 0.02,
        sample_prevalence = 0.5
    conda:
        "envs/ldsc.yaml"
    shell:
        "python $ldsc/ldsc.py --rg {input.sum_stats_A},{input.sum_stats_B} --ref-ld-chr {params.ld_score_root} --w-ld-chr {params.ld_score_root} --out {params.log_file_par} --samp-prev {params.sample_prevalence},{params.sample_prevalence} --pop-prev {params.population_prevalence},{params.population_prevalence} --intercept-h2 1,1"

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

        # TODO some of the variants specified seems to be missing from the cv_dat, 1_100, 5_108, 6_6, 6_25, 6_38, 8_22
rule calculate_theoretical_rg:
    input:
        cv_file = "results/simgwas/combined_causal_variants.tsv",
        available_blocks_file = "resources/simgwas/available_blocks.tsv"
    output:
        "results/ldsc/rg/whole_genome/{effect_blocks_A}_{effect_blocks_B}_theo_rg.tsv"
    params:
        population_prevalence_A = 0.02,
        sample_prevalence_A = 0.5,
        population_prevalence_B = 0.02,
        sample_prevalence_B = 0.5,
        # TODO handle null and non-uniform effect cases
        odds_ratio_a = lambda wildcards: odds_ratio_dict[re.search("[smlvhin]", wildcards.effect_blocks_A).group()],
        odds_ratio_b = lambda wildcards: odds_ratio_dict[re.search("[smlvhin]", wildcards.effect_blocks_B).group()]
    threads: 4
    shell:
        "Rscript workflow/scripts/ldsc/calculate_theoretical_rg.R --cv_file {input.cv_file} --blocks_file {input.available_blocks_file} --effect_blocks_a {wildcards.effect_blocks_A} --effect_blocks_b {wildcards.effect_blocks_B} --odds_ratio_a {params.odds_ratio_a} --odds_ratio_b {params.odds_ratio_b} --P_a {params.sample_prevalence_A} --P_b {params.sample_prevalence_B} --K_a {params.population_prevalence_A} --K_b {params.population_prevalence_B} -o {output} -nt {threads}"

rule munge_whole_genome_sum_stats:
    input:
        "results/simgwas/simulated_sum_stats/whole_genome_sum_stats/{ncases}_{ncontrols}/{effect_blocks}_sum_stats.tsv.gz"
    output:
        # NB: Awkward output filename due to the way the LDSC script works
        temp("results/ldsc/munged_sum_stats/whole_genome_sum_stats/{ncases}_{ncontrols}/{effect_blocks}_{tag,[abcde]}.tsv.sumstats.gz")
    params:
         output_filename = "results/ldsc/munged_sum_stats/whole_genome_sum_stats/{ncases}_{ncontrols}/{effect_blocks}_{tag}.tsv",
         signed_sumstats_col = lambda wildcards: tag_signed_sumstats_dict[wildcards.tag],
         pvalue_col = lambda wildcards: tag_pvalue_dict[wildcards.tag]
    resources:
        disk_mb = 2048
    conda:
        "envs/ldsc.yaml"
    shell:
        """
        python $ldsc/munge_sumstats.py --sumstats {input} --N-con-col ncontrols --N-cas-col ncases --snp id --out {params.output_filename} --signed-sumstats {params.signed_sumstats_col} --p {params.pvalue_col} --a1 a0 --a2 a1 --frq EUR --a1-inc;
        """

rule estimate_single_chr_h2:
    input:
        ["resources/ldsc/eur_w_ld_chr/%d.l2.ldscore.gz" % i for i in range(1,23)],
        ["resources/ldsc/eur_w_ld_chr/%d.l2.M_5_50" % i for i in range(1,23)],
        sum_stats = "results/ldsc/munged_sum_stats/chr{ch}/{ncases}_{ncontrols}/{effect_blocks}_{tag}.tsv.sumstats.gz"
    output:
        "results/ldsc/h2/chr{ch}/{ncases}_{ncontrols}/{effect_blocks}_{tag,[abcde]}.log"
    params:
        log_file_par = "results/ldsc/h2/chr{ch}/{ncases}_{ncontrols}/{effect_blocks}_{tag,[abcde]}",
        # NB: Trailing '/' is needed in ld_score_root
        ld_score_root = "resources/ldsc/eur_w_ld_chr/",
        population_prevalence = 0.02,
        sample_prevalence = 0.5
    conda:
        "envs/ldsc.yaml"
    shell:
        "python $ldsc/ldsc.py --h2 {input.sum_stats} --ref-ld-chr {params.ld_score_root} --w-ld-chr {params.ld_score_root} --out {params.log_file_par} --samp-prev {params.sample_prevalence} --pop-prev {params.population_prevalence}"

rule estimate_single_chr_h2_B0_1:
    input:
        ["resources/ldsc/eur_w_ld_chr/%d.l2.ldscore.gz" % i for i in range(1,23)],
        ["resources/ldsc/eur_w_ld_chr/%d.l2.M_5_50" % i for i in range(1,23)],
        sum_stats = "results/ldsc/munged_sum_stats/chr{ch}/{ncases}_{ncontrols}/{effect_blocks}_{tag}.tsv.sumstats.gz"
    output:
        "results/ldsc/h2/chr{ch}/{ncases}_{ncontrols}/B0_1/{effect_blocks}_{tag,[abcde]}.log"
    params:
        log_file_par = "results/ldsc/h2/chr{ch}/{ncases}_{ncontrols}/B0_1/{effect_blocks}_{tag,[abcde]}",
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

rule estimate_whole_genome_h2:
    input:
        ["resources/ldsc/eur_w_ld_chr/%d.l2.ldscore.gz" % i for i in range(1,23)],
        ["resources/ldsc/eur_w_ld_chr/%d.l2.M_5_50" % i for i in range(1,23)],
        sum_stats = "results/ldsc/munged_sum_stats/whole_genome_sum_stats/{ncases}_{ncontrols}/{effect_blocks}_{tag}.tsv.sumstats.gz"
    output:
        "results/ldsc/h2/whole_genome/{ncases,\d+}_{ncontrols,\d+}/{effect_blocks,[^/]+}_{tag,[abcde]}.log"
    params:
        log_file_par = "results/ldsc/h2/whole_genome/{ncases}_{ncontrols}/{effect_blocks}_{tag,[abcde]}",
        # NB: Trailing '/' is needed in ld_score_root
        ld_score_root = "resources/ldsc/eur_w_ld_chr/",
        population_prevalence = 0.02,
        sample_prevalence = 0.5
    resources:
        disk_mb = 2048
    conda:
        "envs/ldsc.yaml"
    shell:
        "python $ldsc/ldsc.py --h2 {input.sum_stats} --ref-ld-chr {params.ld_score_root} --w-ld-chr {params.ld_score_root} --out {params.log_file_par} --samp-prev {params.sample_prevalence} --pop-prev {params.population_prevalence}"

rule estimate_whole_genome_h2_B0_1:
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

rule estimate_single_chr_rg:
    input:
        ["resources/ldsc/eur_w_ld_chr/%d.l2.ldscore.gz" % i for i in range(1,23)],
        ["resources/ldsc/eur_w_ld_chr/%d.l2.M_5_50" % i for i in range(1,23)],
        sum_stats_A = "results/ldsc/munged_sum_stats/chr{ch_A}/{ncases_A}_{ncontrols_A}/{effect_blocks_A}_a.tsv.sumstats.gz",
        sum_stats_B = "results/ldsc/munged_sum_stats/chr{ch_B}/{ncases_B}_{ncontrols_B}/{effect_blocks_B}_b.tsv.sumstats.gz",
    output:
        "results/ldsc/rg/chr{ch_A}/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}.log"
    params:
        log_file_par = "results/ldsc/rg/chr{ch_A}_{ch_B}/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}",
        # NB: Trailing '/' is needed in ld_score_root
        ld_score_root = "resources/ldsc/eur_w_ld_chr/",
        population_prevalence = 0.02,
        sample_prevalence = 0.5
    conda:
        "envs/ldsc.yaml"
    shell:
        "python $ldsc/ldsc.py --rg {input.sum_stats_A},{input.sum_stats_B} --ref-ld-chr {params.ld_score_root} --w-ld-chr {params.ld_score_root} --out {params.log_file_par} --samp-prev {params.sample_prevalence} --pop-prev {params.population_prevalence}"

rule estimate_whole_genome_rg:
    input:
        ["resources/ldsc/eur_w_ld_chr/%d.l2.ldscore.gz" % i for i in range(1,23)],
        ["resources/ldsc/eur_w_ld_chr/%d.l2.M_5_50" % i for i in range(1,23)],
        sum_stats_A = "results/ldsc/munged_sum_stats/whole_genome_sum_stats/{ncases_A}_{ncontrols_A}/{effect_blocks_A}_{tag_A}.tsv.sumstats.gz",
        sum_stats_B = "results/ldsc/munged_sum_stats/whole_genome_sum_stats/{ncases_B}_{ncontrols_B}/{effect_blocks_B}_{tag_B}.tsv.sumstats.gz",
    output:
        "results/ldsc/rg/whole_genome/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{tag_A,[abcde]}{tag_B,[abcde]}.log"
    params:
        log_file_par = "results/ldsc/rg/whole_genome/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{tag_A}{tag_B}",
        # NB: Trailing '/' is needed in ld_score_root
        ld_score_root = "resources/ldsc/eur_w_ld_chr/",
        population_prevalence = 0.02,
        sample_prevalence = 0.5
    conda:
        "envs/ldsc.yaml"
    shell:
        "python $ldsc/ldsc.py --rg {input.sum_stats_A},{input.sum_stats_B} --ref-ld-chr {params.ld_score_root} --w-ld-chr {params.ld_score_root} --out {params.log_file_par} --samp-prev {params.sample_prevalence},{params.sample_prevalence} --pop-prev {params.population_prevalence},{params.population_prevalence} --intercept-h2 1,1"
