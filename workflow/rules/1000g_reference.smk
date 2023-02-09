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
  shell:
    "Rscript workflow/scripts/get_euro_fam_file.R --panel_file {input.panel} --ped_file {input.ped} --output_file {output}"

rule get_euro_samples:
  input:
    "resources/1000g/{chr}.bed",
    "resources/1000g/{chr}.bim",
    "resources/1000g/{chr}.fam",
    "resources/1000g/euro.fam"
  output:
    temp("resources/1000g/euro/{chr}_euro.bed"),
    temp("resources/1000g/euro/{chr}_euro.bim"),
    temp("resources/1000g/euro/{chr}_euro.fam")
  threads: 8
  resources:
      mem_mb=get_mem_mb
  shell:
    "plink2 --memory {resources.mem_mb} --threads {threads} --bfile resources/1000g/{wildcards.chr} --keep resources/1000g/euro.fam --make-bed --silent --out resources/1000g/euro/{wildcards.chr}_euro"

rule qc:
    input:
        "resources/1000g/euro/{chr}_euro.bed",
        "resources/1000g/euro/{chr}_euro.bim",
        "resources/1000g/euro/{chr}_euro.fam"
    output:
        "resources/1000g/euro/qc/{chr}_qc.bed",
        "resources/1000g/euro/qc/{chr}_qc.bim",
        "resources/1000g/euro/qc/{chr}_qc.fam"
    threads: 8
    resources:
        mem_mb=get_mem_mb
    shell:
        "plink2 --memory {resources.mem_mb} --threads {threads} --bfile resources/1000g/euro/{wildcards.chr}_euro --geno 0.1 --mind 0.1 --maf 0.005 --hwe 1e-50 --make-bed --silent --out resources/1000g/euro/qc/{wildcards.chr}_qc"

rule subset_reference:
    input:
      "resources/1000g/euro/qc/{chr}_qc.bed",
      "resources/1000g/euro/qc/{chr}_qc.bim",
      "resources/1000g/euro/qc/{chr}_qc.fam",
      range_file = "resources/plink_ranges/{join}/{snp_set}/{chr}.txt"
    output:
      "resources/plink_subsets/{join}/{snp_set}/{chr}.bed",
      "resources/plink_subsets/{join}/{snp_set}/{chr}.bim",
      "resources/plink_subsets/{join}/{snp_set}/{chr}.fam"
    params:
      bfile = "resources/1000g/euro/qc/{chr}_qc",
      out = "resources/plink_subsets/{join}/{snp_set}/{chr}"
    threads: 8
    resources:
        mem_mb=get_mem_mb
    shell:
      "plink2 --memory {resources.mem_mb} --threads {threads} --bfile {params.bfile} --extract {input.range_file} --make-bed --out {params.out}"

rule make_pruned_ranges:
    input:
      "resources/plink_subsets/{join}/{snp_set}/{chr}.bed",
      "resources/plink_subsets/{join}/{snp_set}/{chr}.bim",
      "resources/plink_subsets/{join}/{snp_set}/{chr}.fam"
    output:
      "resources/plink_ranges/{join}/{snp_set}/pruned_ranges/window_{window}_step_{step}_r2_{r2}/{chr}.prune.in",
      "resources/plink_ranges/{join}/{snp_set}/pruned_ranges/window_{window}_step_{step}_r2_{r2}/{chr}.prune.out"
    params:
      bfile = "resources/plink_subsets/{join}/{snp_set}/{chr}",
      prune_out = "resources/plink_ranges/{join}/{snp_set}/pruned_ranges/window_{window}_step_{step}_r2_{r2}/{chr}",
      r2 = lambda wildcards: float(wildcards.r2.replace('_', '.'))
    threads: 8
    resources:
        mem_mb=get_mem_mb
    shell:
      "plink --memory {resources.mem_mb} --threads {threads} --bfile {params.bfile} --indep-pairwise {wildcards.window} {wildcards.step} {params.r2} --out {params.prune_out}"

rule cat_pruned_ranges:
    input:
      ("resources/plink_ranges/{join}/{snp_set}/pruned_ranges/window_{window}_step_{step}_r2_{r2}/chr%d.prune.in" % x for x in range(1,23))
    output:
        "resources/plink_ranges/{join}/{snp_set}/pruned_ranges/window_{window}_step_{step}_r2_{r2}/all.prune.in"
    shell:
      "for x in {input}; do cat $x >>{output}; done"

rule cat_bim_files:
    input:
        ("resources/plink_subsets/{join}/{snp_set}/chr%d.bim" % x for x in range(1,23))
    output:
        "resources/plink_subsets/{join}/{snp_set}/all.bim"
    shell:
        "for x in {input}; do cat $x >>{output}; done"

rule deduplicate_variants:
    input:
        "resources/1000g/euro/qc/{chr}_qc.bed",
        "resources/1000g/euro/qc/{chr}_qc.bim",
        "resources/1000g/euro/qc/{chr}_qc.fam"
    output:
        "resources/1000g/euro/qc/nodup/{chr}.bed",
        "resources/1000g/euro/qc/nodup/{chr}.bim",
        "resources/1000g/euro/qc/nodup/{chr}.fam"
    params:
        input_stem = "resources/1000g/euro/qc/{chr}_qc",
        output_stem = "resources/1000g/euro/qc/nodup/{chr}"
    threads: 8
    resources:
        mem_mb=get_mem_mb
    shell:
        "plink2 --memory {resources.mem_mb} --threads {threads} --bfile {params.input_stem} --rm-dup 'force-first' --make-bed --silent --out {params.output_stem}"

rule subset_snp_variants:
    input:
        "resources/1000g/euro/qc/nodup/{chr}.bed",
        "resources/1000g/euro/qc/nodup/{chr}.bim",
        "resources/1000g/euro/qc/nodup/{chr}.fam",
        range_file = "resources/plink_ranges/{join}/{snp_set}/{chr}.txt"
    output:
        "resources/1000g/euro/qc/nodup/snps_only/{join}/{snp_set}/{chr}.bed",
        "resources/1000g/euro/qc/nodup/snps_only/{join}/{snp_set}/{chr}.bim",
        "resources/1000g/euro/qc/nodup/snps_only/{join}/{snp_set}/{chr}.fam"
    params:
        input_stem = "resources/1000g/euro/qc/nodup/{chr}",
        output_stem = "resources/1000g/euro/qc/nodup/snps_only/{join}/{snp_set}{chr}",
    threads: 8
    resources:
        mem_mb=get_mem_mb
    shell:
        "plink2 --memory {resources.mem_mb} --threads {threads} --bfile {params.input_stem} --snps-only --make-bed --extract {input.range_file} --silent --out {params.output_stem}"

rule concatenate_qc_bim_files:
    input:
        ["resources/1000g/euro/qc/chr%d_qc.bim" % x for x in range(1,23)]
    output:
        "resources/1000g/euro/qc/chr1-22_qc.bim"
    run:
        for i,x in enumerate(input):
            shell("cat %s >> %s" % (x, output[0]))

# NB: 'make_pruned_ranges' handles a subset of 1kGP SNPs, this prunes all QC-passing 1kGP SNPs for use with simGWAS
rule make_1000g_pruned_ranges:
    input:
        "resources/1000g/euro/qc/{chr}_qc.bed",
        "resources/1000g/euro/qc/{chr}_qc.bim",
        "resources/1000g/euro/qc/{chr}_qc.fam"
    output:
        "resources/1000g/euro/qc/pruned_ranges/window_{window}_step_{step}_r2_{r2}/{chr}.prune.in",
        "resources/1000g/euro/qc/pruned_ranges/window_{window}_step_{step}_r2_{r2}/{chr}.prune.out"
    params:
        bfile = "resources/1000g/euro/qc/{chr}_qc",
        prune_out = "resources/1000g/euro/qc/pruned_ranges/window_{window}_step_{step}_r2_{r2}/{chr}",
        r2 = lambda wildcards: float(wildcards.r2.replace('_', '.'))
    threads: 8
    resources:
        mem_mb=get_mem_mb
    shell:
        "plink --memory {resources.mem_mb} --threads {threads} --bfile {params.bfile} --indep-pairwise {wildcards.window} {wildcards.step} {params.r2} --out {params.prune_out}"

rule prune_qced_1000g_snps:
    input:
     [("resources/1000g/euro/qc/pruned_ranges/window_1000kb_step_50_r2_0_2/chr%s.prune.in" % x for x in range(1,23))]+
     [("resources/1000g/euro/qc/pruned_ranges/window_50_step_5_r2_0_2/chr%s.prune.in" % x for x in range(1,23))]
