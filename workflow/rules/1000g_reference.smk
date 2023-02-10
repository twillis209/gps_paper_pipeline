import pandas as pd

rule download_1000g_genotype_data:
  output:
      "resources/1000g/{chr}.vcf.gz"
  run:
    if wildcards.chr == 'chrX':
        shell("wget -O resources/1000g/{wildcards.chr}.vcf.gz       http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.{wildcards.chr}.phase3_shapeit2_mvncall_integrated_v1c.20130502.genotypes.vcf.gz")
    elif wildcards.chr == 'chrY':
        shell("wget -O resources/1000g/{wildcards.chr}.vcf.gz http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chrY.phase3_integrated_v2b.20130502.genotypes.vcf.gz")
    else:
      shell("wget -O resources/1000g/{wildcards.chr}.vcf.gz http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.{wildcards.chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz")

rule download_1000g_sample_metadata:
  output:
      "resources/1000g/integrated_call_samples_v3.20130502.ALL.panel.txt",
      "resources/1000g/integrated_call_samples_v3.20200731.ALL.ped.txt"
  shell:
      """
      wget -O resources/1000g/integrated_call_samples_v3.20130502.ALL.panel.txt http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel
      wget -O resources/1000g/integrated_call_samples_v3.20200731.ALL.ped.txt http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20200731.ALL.ped
      """

rule vcf_to_bed:
  input:
    "resources/1000g/{chr}.vcf.gz"
  output:
    temp("resources/1000g/{chr}.bed"),
    "resources/1000g/{chr}.bim",
    temp("resources/1000g/{chr}.fam")
  params:
    out = "resources/1000g/{chr}"
  threads: 8
  resources:
      mem_mb=get_mem_mb
  shell:
      "plink2 --memory {resources.mem_mb} --threads {threads} --vcf {input} --make-bed --out {params.out} --set-all-var-ids @:#:\$r:\$a --max-alleles 2 --new-id-max-allele-len 20 'truncate'"

# NB: Removes related individuals, so we have 498 of 503 Europeans left
rule make_euro_fam:
  input:
      panel = "resources/1000g/integrated_call_samples_v3.20130502.ALL.panel.txt",
      ped = "resources/1000g/integrated_call_samples_v3.20200731.ALL.ped.txt"
  output:
    "resources/1000g/euro.fam"
  script:
    "../scripts/1000g/get_euro_fam_file.R"

rule get_euro_samples:
  input:
    "resources/1000g/{chr}.bed",
    "resources/1000g/{chr}.bim",
    "resources/1000g/{chr}.fam",
    "resources/1000g/euro.fam"
  output:
    temp("resources/1000g/euro/{chr}.bed"),
    temp("resources/1000g/euro/{chr}.bim"),
    temp("resources/1000g/euro/{chr}.fam")
  threads: 8
  resources:
      mem_mb=get_mem_mb
  shell:
    "plink2 --memory {resources.mem_mb} --threads {threads} --bfile resources/1000g/{wildcards.chr} --keep resources/1000g/euro.fam --make-bed --silent --out resources/1000g/euro/{wildcards.chr}_euro"

rule qc:
    input:
        "resources/1000g/euro/{chr}.bed",
        "resources/1000g/euro/{chr}.bim",
        "resources/1000g/euro/{chr}.fam"
    output:
        "resources/1000g/euro/qc/{chr}.bed",
        "resources/1000g/euro/qc/{chr}.bim",
        "resources/1000g/euro/qc/{chr}.fam"
    threads: 8
    resources:
        mem_mb=get_mem_mb
    shell:
        "plink2 --memory {resources.mem_mb} --threads {threads} --bfile resources/1000g/euro/{wildcards.chr}_euro --geno 0.1 --mind 0.1 --maf 0.005 --hwe 1e-50 --rm-dup 'force-first' --make-bed --silent --out resources/1000g/euro/qc/{wildcards.chr}_qc"

rule make_subset_ranges:
    input:
        bim = "resources/1000g/euro/qc/{chr}.bim",
        ukbb = "resources/ukbb/ukbb_sum_stats/merged_ukbb_sum_stats.tsv.gz"
    output:
        "resources/1000g/euro/qc/{snp_set}/ranges/{chr}.txt"
    threads: 8
    resources:
        mem_mb=get_mem_mb
    script: "../scripts/1000g/make_subset_ranges.R"

rule subset_reference:
    input:
      "resources/1000g/euro/qc/{chr}.bed",
      "resources/1000g/euro/qc/{chr}.bim",
      "resources/1000g/euro/qc/{chr}.fam",
      range_file = "resources/1000g/euro/qc/{snp_set}/ranges/{chr}.txt"
    output:
      "resources/1000g/euro/qc/{snp_set}/all/{chr}.bed",
      "resources/1000g/euro/qc/{snp_set}/all/{chr}.bim",
      "resources/1000g/euro/qc/{snp_set}/all/{chr}.fam"
    params:
      bfile = "resources/1000g/euro/qc/{chr}",
      out = "resources/1000g/euro/qc/{snp_set}/all/{chr}"
    threads: 8
    resources:
        mem_mb=get_mem_mb
    shell:
      "plink2 --memory {resources.mem_mb} --threads {threads} --bfile {params.bfile} --extract {input.range_file} --make-bed --out {params.out}"

rule make_pruned_ranges:
    input:
        "resources/1000g/euro/qc/{snp_set}/{variant_set}/{chr}.bed",
        "resources/1000g/euro/qc/{snp_set}/{variant_set}/{chr}.bim",
        "resources/1000g/euro/qc/{snp_set}/{variant_set}/{chr}.fam"
    output:
        "resources/1000g/euro/qc/{snp_set}/{variant_set}/ranges/prune/window_{window}_step_{step}_r2_{r2}/{chr}.prune.in",
        "resources/1000g/euro/qc/{snp_set}/{variant_set}/ranges/prune/window_{window}_step_{step}_r2_{r2}/{chr}.prune.out"
    params:
      bfile = "resources/1000g/euro/qc/{snp_set}/{variant_set}/{chr}",
      prune_out = "resources/1000g/euro/qc/{snp_set}/{variant_set}/ranges/prune/window_{window}_step_{step}_r2_{r2}/{chr}",
      r2 = lambda wildcards: float(wildcards.r2.replace('_', '.'))
    threads: 8
    resources:
        mem_mb=get_mem_mb
    shell:
      "plink --memory {resources.mem_mb} --threads {threads} --bfile {params.bfile} --indep-pairwise {wildcards.window} {wildcards.step} {params.r2} --out {params.prune_out}"

rule cat_pruned_ranges:
    input:
      ("resources/1000g/euro/qc/{snp_set}/{variant_set}/ranges/prune/window_{window}_step_{step}_r2_{r2}/chr%d.prune.in" % x for x in range(1,23))
    output:
        "resources/1000g/euro/qc/{snp_set}/{variant_set}/ranges/prune/window_{window}_step_{step}_r2_{r2}/all.prune.in"
    shell:
      "for x in {input}; do cat $x >>{output}; done"

rule cat_bim_files:
    input:
        [f"resources/1000g/euro/qc/{{snp_set}}/{{variant_set}}/chr{x}.bim" for x in range(1,23)]
    output:
        "resources/1000g/euro/qc/{snp_set}/{{variant_set}}/all.bim"
    shell:
        "for x in {input}; do cat $x >>{output}; done"

rule subset_to_get_snp_variants_only:
    input:
        "resources/1000g/euro/qc/{snp_set}/all/{chr}.bed",
        "resources/1000g/euro/qc/{snp_set}/all/{chr}.bim",
        "resources/1000g/euro/qc/{snp_set}/all/{chr}.fam",
    output:
        "resources/1000g/euro/qc/{snp_set}/snps_only/{chr}.bed",
        "resources/1000g/euro/qc/{snp_set}/snps_only/{chr}.bim",
        "resources/1000g/euro/qc/{snp_set}/snps_only/{chr}.fam"
    params:
        input_stem = "resources/1000g/euro/qc/{snp_set}/all/{chr}",
        output_stem = "resources/1000g/euro/qc/{snp_set}/snps_only/{chr}"
    threads: 8
    resources:
        mem_mb=get_mem_mb
    shell:
        "plink2 --memory {resources.mem_mb} --threads {threads} --bfile {params.input_stem} --snps-only --make-bed --silent --out {params.output_stem}"
