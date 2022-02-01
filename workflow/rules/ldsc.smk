rule download_ld_scores:
    output:
        "resources/ldsc/eur_w_ld_chr.tar.bz2"
    shell:
        "wget -O {output} https://storage.googleapis.com/broad-alkesgroup-public/LDSCORE/eur_w_ld_chr.tar.bz2"

rule extract_ld_scores:
    input:
        "resources/ldsc/eur_w_ld_chr.tar.bz2"
    output:
        "resources/ldsc/eur_w_ld_chr"
    shell:
        "tar -xjf {input}"
"""
rule munge_sum_stats:
    pass

rule univariate_ldsc:
    input:
      "resources/ldsc/eur_w_ld_chr",
      sum_stats_file = 
    shell:
      "munge_sumstats.py"
"""
