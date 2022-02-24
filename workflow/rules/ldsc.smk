tag_pvalue_dict = {"a": "p.1", "b": "p.2", "c": "p.3", "d": "p.4", "e": "p.5"}
tag_signed_sumstats_dict = {"a": "betasim.1,18", "b": "betasim.2,19", "c": "betasim.3,20", "d": "betasim.4,21", "e": "betasim.5,22"}

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

rule munge_single_chr_sum_stats:
    input:
        "results/simgwas/simulated_sum_stats/chr{ch}/whole_chr_sum_stats/{ncases}_{ncontrols}/{effect_blocks}_sum_stats.tsv.gz"
    output:
        # NB: Awkward output filename due to the way the LDSC script works
         temp("results/ldsc/munged_sum_stats/chr{ch}/{ncases}_{ncontrols}/{effect_blocks}_{tag,[abcde]}.tsv.sumstats.gz")
    params:
         output_filename = "results/ldsc/munged_sum_stats/chr{ch}/{ncases}_{ncontrols}/{effect_blocks}_{tag}.tsv",
         signed_sumstats_col = lambda wildcards: tag_signed_sumstats_dict[wildcards.tag],
         pvalue_col = lambda wildcards: tag_pvalue_dict[wildcards.tag]
    conda:
        "envs/ldsc.yaml"
    shell:
        """
        python $ldsc/munge_sumstats.py --sumstats {input} --N-con-col ncontrols --N-cas-col ncases --snp id --out {params.output_filename} --signed-sumstats {params.signed_sumstats_col} --p {params.pvalue_col} --a1 a0 --a2 a1 --frq EUR --a1-inc;
        """

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

rule estimate_whole_genome_h2:
    input:
        ["resources/ldsc/eur_w_ld_chr/%d.l2.ldscore.gz" % i for i in range(1,23)],
        ["resources/ldsc/eur_w_ld_chr/%d.l2.M_5_50" % i for i in range(1,23)],
        sum_stats = "results/ldsc/munged_sum_stats/whole_genome_sum_stats/{ncases}_{ncontrols}/{effect_blocks}_{tag}.tsv.sumstats.gz"
    output:
        "results/ldsc/h2/whole_genome/{ncases}_{ncontrols}/{effect_blocks}_{tag,[abcde]}.log"
    params:
        log_file_par = "results/ldsc/h2/whole_genome/{ncases}_{ncontrols}/{effect_blocks}_{tag,[abcde]}",
        # NB: Trailing '/' is needed in ld_score_root
        ld_score_root = "resources/ldsc/eur_w_ld_chr/",
        population_prevalence = 0.02,
        sample_prevalence = 0.5
    conda:
        "envs/ldsc.yaml"
    shell:
        "python $ldsc/ldsc.py --h2 {input.sum_stats} --ref-ld-chr {params.ld_score_root} --w-ld-chr {params.ld_score_root} --out {params.log_file_par} --samp-prev {params.sample_prevalence} --pop-prev {params.population_prevalence}"

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
        sum_stats_A = "results/ldsc/munged_sum_stats/whole_genome_sum_stats/{ncases_A}_{ncontrols_A}/{effect_blocks_A}_a.tsv.sumstats.gz",
        sum_stats_B = "results/ldsc/munged_sum_stats/whole_genome_sum_stats/{ncases_B}_{ncontrols_B}/{effect_blocks_B}_b.tsv.sumstats.gz",
    output:
        "results/ldsc/rg/whole_genome/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}.log"
    params:
        log_file_par = "results/ldsc/rg/whole_genome/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}",
        # NB: Trailing '/' is needed in ld_score_root
        ld_score_root = "resources/ldsc/eur_w_ld_chr/",
        population_prevalence = 0.02,
        sample_prevalence = 0.5
    conda:
        "envs/ldsc.yaml"
    shell:
        "python $ldsc/ldsc.py --rg {input.sum_stats_A},{input.sum_stats_B} --ref-ld-chr {params.ld_score_root} --w-ld-chr {params.ld_score_root} --out {params.log_file_par} --samp-prev {params.sample_prevalence} --pop-prev {params.population_prevalence}"
