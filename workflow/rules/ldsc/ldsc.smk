from string import ascii_lowercase
import pandas as pd

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
