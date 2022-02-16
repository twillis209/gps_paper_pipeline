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

rule munge_sum_stats:
    input:
        "results/simgwas/simulated_sum_stats/chr{ch}/whole_chr_sum_stats/{ncases}_{ncontrols}/blocks_{first}_{last}_{effect_block}_sum_stats.tsv.gz"
    output:
        # NB: Awkward output filename due to the way the LDSC script works
         temp("results/ldsc/munged_sum_stats/chr{ch}/{ncases}_{ncontrols}/blocks_{first}_{last}_{effect_block}_{tag,[ab]}.tsv.sumstats.gz")
    params:
         output_filename = "results/ldsc/munged_sum_stats/chr{ch}/{ncases}_{ncontrols}/blocks_{first}_{last}_{effect_block}_{tag}.tsv",
         signed_sumstats_col = lambda wildcards: "betasim.1,12" if wildcards.tag == 'a' else "betasim.2,13",
         pvalue_col = lambda wildcards: "p.1" if wildcards.tag == 'a' else "p.2"
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
        temp("results/ldsc/munged_sum_stats/whole_genome_sum_stats/{ncases}_{ncontrols}/{effect_blocks}_{tag,[ab]}.tsv.sumstats.gz")
    params:
         output_filename = "results/ldsc/munged_sum_stats/whole_genome_sum_stats/{ncases}_{ncontrols}/{effect_blocks}_{tag}.tsv",
         signed_sumstats_col = lambda wildcards: "betasim.1,12" if wildcards.tag == 'a' else "betasim.2,13",
         pvalue_col = lambda wildcards: "p.1" if wildcards.tag == 'a' else "p.2"
    conda:
        "envs/ldsc.yaml"
    shell:
        """
        python $ldsc/munge_sumstats.py --sumstats {input} --N-con-col ncontrols --N-cas-col ncases --snp id --out {params.output_filename} --signed-sumstats {params.signed_sumstats_col} --p {params.pvalue_col} --a1 a0 --a2 a1 --frq EUR --a1-inc;
        """

rule estimate_h2:
    input:
        ["resources/ldsc/eur_w_ld_chr/%d.l2.ldscore.gz" % i for i in range(1,23)],
        ["resources/ldsc/eur_w_ld_chr/%d.l2.M_5_50" % i for i in range(1,23)],
        sum_stats = "results/ldsc/munged_sum_stats/chr{ch}/{ncases}_{ncontrols}/blocks_{first}_{last}_{effect_block}_{tag}.tsv.sumstats.gz"
    output:
        "results/ldsc/h2/chr{ch}/{ncases}_{ncontrols}/blocks_{first}_{last}_{effect_block}_{tag,[ab]}.log"
    params:
        log_file_par = "results/ldsc/h2/chr{ch}/{ncases}_{ncontrols}/blocks_{first}_{last}_{effect_block}_{tag,[ab]}",
        # NB: Trailing '/' is needed in ld_score_root
        ld_score_root = "resources/ldsc/eur_w_ld_chr/"
    conda:
        "envs/ldsc.yaml"
    shell:
        "python $ldsc/ldsc.py --h2 {input.sum_stats} --ref-ld-chr {params.ld_score_root} --w-ld-chr {params.ld_score_root} --out {params.log_file_par}"

rule estimate_whole_genome_h2:
    input:
        ["resources/ldsc/eur_w_ld_chr/%d.l2.ldscore.gz" % i for i in range(1,23)],
        ["resources/ldsc/eur_w_ld_chr/%d.l2.M_5_50" % i for i in range(1,23)],
        sum_stats = "results/ldsc/munged_sum_stats/whole_genome_sum_stats/{ncases}_{ncontrols}/{effect_blocks}_{tag}.tsv.sumstats.gz"
    output:
        "results/ldsc/h2/whole_genome/{ncases}_{ncontrols}/{effect_blocks}_{tag,[ab]}.log"
    params:
        log_file_par = "results/ldsc/h2/whole_genome/{ncases}_{ncontrols}/{effect_blocks}_{tag,[ab]}",
        # NB: Trailing '/' is needed in ld_score_root
        ld_score_root = "resources/ldsc/eur_w_ld_chr/"
    conda:
        "envs/ldsc.yaml"
    shell:
        "python $ldsc/ldsc.py --h2 {input.sum_stats} --ref-ld-chr {params.ld_score_root} --w-ld-chr {params.ld_score_root} --out {params.log_file_par}"

rule estimate_rg:
    input:
        ["resources/ldsc/eur_w_ld_chr/%d.l2.ldscore.gz" % i for i in range(1,23)],
        ["resources/ldsc/eur_w_ld_chr/%d.l2.M_5_50" % i for i in range(1,23)],
        sum_stats_A = "results/ldsc/munged_sum_stats/chr{ch_A}/{ncases_A}_{ncontrols_A}/blocks_{first_A}_{last_A}_{effect_block_A}_a.tsv.sumstats.gz",
        sum_stats_B = "results/ldsc/munged_sum_stats/chr{ch_B}/{ncases_B}_{ncontrols_B}/blocks_{first_B}_{last_B}_{effect_block_B}_b.tsv.sumstats.gz",
    output:
        "results/ldsc/rg/chr{ch_A}_{ch_B}/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/blocks_{first_A}_{last_A}_{effect_block_A}_{first_B}_{last_B}_{effect_block_B}.log"
    params:
        log_file_par = "results/ldsc/rg/chr{ch_A}_{ch_B}/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/blocks_{first_A}_{last_A}_{effect_block_A}_{first_B}_{last_B}_{effect_block_B}",
        # NB: Trailing '/' is needed in ld_score_root
        ld_score_root = "resources/ldsc/eur_w_ld_chr/"
    conda:
        "envs/ldsc.yaml"
    shell:
        "python $ldsc/ldsc.py --rg {input.sum_stats_A},{input.sum_stats_B} --ref-ld-chr {params.ld_score_root} --w-ld-chr {params.ld_score_root} --out {params.log_file_par}"
