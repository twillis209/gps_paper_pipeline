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

# TODO do something about the $ldsc environment variable; ldsc as a submodule?
rule munge_sum_stats:
    input:
        "results/simgwas/simulated_sum_stats/chr{ch}/whole_chr_sum_stats/{effect_size}/blocks_{first}_{last}_effect_{effect_block}_sum_stats.tsv.gz"
    output:
        # NB: Awkward output filename due to the way the LDSC script works
        "results/ldsc/munged_sum_stats/chr{ch}/{effect_size}/blocks_{first}_{last}_effect_{effect_block}.tsv.sumstats.gz"
    params:
        output_filename = "results/ldsc/munged_sum_stats/chr{ch}/{effect_size}/blocks_{first}_{last}_effect_{effect_block}.tsv"
    conda:
        "envs/ldsc.yaml"
    shell:
        "python $ldsc/munge_sumstats.py --sumstats {input} --N-con-col ncontrols --N-cas-col ncases --snp id --out {params.output_filename} --signed-sumstats betasim.1,9 --p p.1 --a1 a0 --a2 a1 --frq EUR --a1-inc"

rule estimate_h2:
    input:
        ["resources/ldsc/eur_w_ld_chr/%d.l2.ldscore.gz" % i for i in range(1,23)],
        ["resources/ldsc/eur_w_ld_chr/%d.l2.M_5_50" % i for i in range(1,23)],
        sum_stats = "results/ldsc/munged_sum_stats/chr{ch}/{effect_size}/blocks_{first}_{last}_effect_{effect_block}.tsv.sumstats.gz"
    output:
        sum_stats = "results/ldsc/h2/chr{ch}/{effect_size}/blocks_{first}_{last}_effect_{effect_block}.log"
    params:
        log_file_par = "results/ldsc/h2/chr{ch}/{effect_size}/blocks_{first}_{last}_effect_{effect_block}",
        # NB: Trailing '/' is needed in ld_score_root
        ld_score_root = "resources/ldsc/eur_w_ld_chr/"
    conda:
        "envs/ldsc.yaml"
    shell:
        "python $ldsc/ldsc.py --h2 {input.sum_stats} --ref-ld-chr {params.ld_score_root} --w-ld-chr {params.ld_score_root} --out {params.log_file_par}"
