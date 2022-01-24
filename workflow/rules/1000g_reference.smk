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
    "resources/1000g/{chr}.bed",
    "resources/1000g/{chr}.bim",
    "resources/1000g/{chr}.fam"
  threads: 8
  resources:
      mem_mb=get_mem_mb
  shell:
    "plink --memory {resources.mem_mb} --threads {threads} --vcf resources/1000g/{wildcards.chr}.vcf.gz --make-bed --out resources/1000g/{wildcards.chr}"

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
    "resources/1000g/euro/{chr}_euro.bed",
    "resources/1000g/euro/{chr}_euro.bim",
    "resources/1000g/euro/{chr}_euro.fam"
  threads: 8
  resources:
      mem_mb=get_mem_mb
  shell:
    "plink --memory {resources.mem_mb} --threads {threads} --bfile resources/1000g/{wildcards.chr} --keep resources/1000g/euro.fam --make-bed --silent --out resources/1000g/euro/{wildcards.chr}_euro"

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
    "plink --memory {resources.mem_mb} --threads {threads} --bfile resources/1000g/euro/{wildcards.chr}_euro --geno 0.1 --mind 0.1 --maf 0.005 --hwe 1e-50 --make-bed --silent --out resources/1000g/euro/qc/{wildcards.chr}_qc"
