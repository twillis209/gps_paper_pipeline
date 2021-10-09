wildcard_constraints:
    chr = "chr[0-9XY]{1,2}"

rule download_1000g_genotype_data:
  # hg19
  # TODO Really not sure that these files work with Plink v1.9
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
  shell:
    "plink --memory 16000 --threads 8 --vcf resources/1000g/{wildcards.chr}.vcf.gz --make-bed --out resources/1000g/{wildcards.chr}"

rule make_euro_fam:
  input:
      panel = "resources/1000g/integrated_call_samples_v3.20130502.ALL.panel.txt",
      ped = "resources/1000g/integrated_call_samples_v3.20200731.ALL.ped.txt"
  output:
    "resources/1000g/euro.fam"
  shell:
    "Rscript scripts/get_euro_fam_file.R --panel_file {input.panel} --ped_file {input.ped} --output_file {output}"

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
  shell:
    "plink --memory 16000 --threads 8 --bfile resources/1000g/{wildcards.chr} --keep resources/1000g/euro.fam --make-bed --silent --out resources/1000g/euro/{wildcards.chr}_euro"

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
  shell:
    "plink --memory 16000 --threads 8 --bfile resources/1000g/euro/{wildcards.chr}_euro --geno 0.1 --mind 0.1 --maf 0.005 --hwe 1e-50 --make-bed --silent --out resources/1000g/euro/qc/{wildcards.chr}_qc"

rule download_ukbb_sum_stats:
    output:
        "resources/ukbb_sum_stats/20002_1111.gwas.imputed_v3.both_sexes.tsv.bgz",
        "resources/ukbb_sum_stats/20002_1113.gwas.imputed_v3.both_sexes.tsv.bgz",
        "resources/ukbb_sum_stats/20002_1154.gwas.imputed_v3.both_sexes.tsv.bgz",
        "resources/ukbb_sum_stats/20002_1220.gwas.imputed_v3.both_sexes.tsv.bgz",
        "resources/ukbb_sum_stats/20002_1226.gwas.imputed_v3.both_sexes.tsv.bgz",
        "resources/ukbb_sum_stats/20002_1286.gwas.imputed_v3.both_sexes.tsv.bgz",
        "resources/ukbb_sum_stats/20002_1452.gwas.imputed_v3.both_sexes.tsv.bgz",
        "resources/ukbb_sum_stats/20002_1462.gwas.imputed_v3.both_sexes.tsv.bgz",
        "resources/ukbb_sum_stats/20002_1464.gwas.imputed_v3.both_sexes.tsv.bgz",
        "resources/ukbb_sum_stats/20002_1465.gwas.imputed_v3.both_sexes.tsv.bgz",
        "resources/ukbb_sum_stats/20002_1473.gwas.imputed_v3.both_sexes.tsv.bgz",
        "resources/ukbb_sum_stats/22126.gwas.imputed_v3.both_sexes.tsv.bgz",
        "resources/ukbb_sum_stats/6148_2.gwas.imputed_v3.both_sexes.tsv.bgz",
        "resources/ukbb_sum_stats/6148_5.gwas.imputed_v3.both_sexes.tsv.bgz",
        "resources/ukbb_sum_stats/D25.gwas.imputed_v3.both_sexes.tsv.bgz",
        "resources/ukbb_sum_stats/I9_IHD.gwas.imputed_v3.both_sexes.tsv.bgz",
        "resources/ukbb_sum_stats/K51.gwas.imputed_v3.both_sexes.tsv.bgz",
        "resources/ukbb_sum_stats/K57.gwas.imputed_v3.both_sexes.tsv.bgz",
        "resources/ukbb_sum_stats/K80.gwas.imputed_v3.both_sexes.tsv.bgz"
    shell:
        """
        wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/20002_1111.gwas.imputed_v3.both_sexes.tsv.bgz -O resources/ukbb_sum_stats/20002_1111.gwas.imputed_v3.both_sexes.tsv.bgz	
        wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/20002_1113.gwas.imputed_v3.both_sexes.tsv.bgz -O resources/ukbb_sum_stats/20002_1113.gwas.imputed_v3.both_sexes.tsv.bgz	
        wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/20002_1154.gwas.imputed_v3.both_sexes.tsv.bgz -O resources/ukbb_sum_stats/20002_1154.gwas.imputed_v3.both_sexes.tsv.bgz	
        wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/20002_1220.gwas.imputed_v3.both_sexes.tsv.bgz -O resources/ukbb_sum_stats/20002_1220.gwas.imputed_v3.both_sexes.tsv.bgz	
        wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/20002_1226.gwas.imputed_v3.both_sexes.tsv.bgz -O resources/ukbb_sum_stats/20002_1226.gwas.imputed_v3.both_sexes.tsv.bgz	
        wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/20002_1286.gwas.imputed_v3.both_sexes.tsv.bgz -O resources/ukbb_sum_stats/20002_1286.gwas.imputed_v3.both_sexes.tsv.bgz	
        wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/20002_1452.gwas.imputed_v3.both_sexes.tsv.bgz -O resources/ukbb_sum_stats/20002_1452.gwas.imputed_v3.both_sexes.tsv.bgz	
        wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/20002_1462.gwas.imputed_v3.both_sexes.tsv.bgz -O resources/ukbb_sum_stats/20002_1462.gwas.imputed_v3.both_sexes.tsv.bgz	
        wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/20002_1464.gwas.imputed_v3.both_sexes.tsv.bgz -O resources/ukbb_sum_stats/20002_1464.gwas.imputed_v3.both_sexes.tsv.bgz	
        wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/20002_1465.gwas.imputed_v3.both_sexes.tsv.bgz -O resources/ukbb_sum_stats/20002_1465.gwas.imputed_v3.both_sexes.tsv.bgz	
        wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/20002_1473.gwas.imputed_v3.both_sexes.tsv.bgz -O resources/ukbb_sum_stats/20002_1473.gwas.imputed_v3.both_sexes.tsv.bgz	
        wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/22126.gwas.imputed_v3.both_sexes.tsv.bgz -O resources/ukbb_sum_stats/22126.gwas.imputed_v3.both_sexes.tsv.bgz	
        wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/6148_2.gwas.imputed_v3.both_sexes.tsv.bgz -O resources/ukbb_sum_stats/6148_2.gwas.imputed_v3.both_sexes.tsv.bgz	
        wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/6148_5.gwas.imputed_v3.both_sexes.tsv.bgz -O resources/ukbb_sum_stats/6148_5.gwas.imputed_v3.both_sexes.tsv.bgz	
        wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/D25.gwas.imputed_v3.both_sexes.tsv.bgz -O resources/ukbb_sum_stats/D25.gwas.imputed_v3.both_sexes.tsv.bgz	
        wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/I9_IHD.gwas.imputed_v3.both_sexes.tsv.bgz -O resources/ukbb_sum_stats/I9_IHD.gwas.imputed_v3.both_sexes.tsv.bgz	
        wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/K51.gwas.imputed_v3.both_sexes.tsv.bgz -O resources/ukbb_sum_stats/K51.gwas.imputed_v3.both_sexes.tsv.bgz	
        wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/K57.gwas.imputed_v3.both_sexes.tsv.bgz -O resources/ukbb_sum_stats/K57.gwas.imputed_v3.both_sexes.tsv.bgz	
        wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/K80.gwas.imputed_v3.both_sexes.tsv.bgz -O resources/ukbb_sum_stats/K80.gwas.imputed_v3.both_sexes.tsv.bgz
        """

"""
ceeabfc0db7de5400f9a5ab1ad623f92
58e844f62b9ed983240cd0ffbab27aec
e647700a08c45ccb43cd6f23162dd6a2
05ff0825233a942af4a2882af49095ef
1bd4e7f6df2b74a20c20d61cea4236f7
08ea662fb003101dafbd87aeedff155c
d4b0fd4d8b33e9dc095ad673116448eb
42608dee500ae5d5ffd494cdd7540147
590f210d10c7b07d16a662e04b0ce21f
1f5a5675725d2675632a2339a6b556b8
32a1fee9b94a2670eff787ce9958ee49
96530bc3fe31c217833bbbf20bb4992f
f2ba980abfbd5d8dcb5a77ba4a38df26
c4aa8e26d88b28f96418ac380a2a60f1
f21035b64cc567eaff5a6af886d32f14
9cc7e9b25cbbf74f13321cac42753bd5
8acddc62ae088aa3fb962ba3af485ba1
6668642bbacc7f59d91c1a5d6a3bb1f7
6a7160c17c71e4c95415249f76ad96c1
"""

rule decompress_ukbb_sum_stats:
    input:
        "resources/ukbb_sum_stats/20002_1111.gwas.imputed_v3.both_sexes.tsv.bgz",
        "resources/ukbb_sum_stats/20002_1113.gwas.imputed_v3.both_sexes.tsv.bgz",
        "resources/ukbb_sum_stats/20002_1154.gwas.imputed_v3.both_sexes.tsv.bgz",
        "resources/ukbb_sum_stats/20002_1220.gwas.imputed_v3.both_sexes.tsv.bgz",
        "resources/ukbb_sum_stats/20002_1226.gwas.imputed_v3.both_sexes.tsv.bgz",
        "resources/ukbb_sum_stats/20002_1286.gwas.imputed_v3.both_sexes.tsv.bgz",
        "resources/ukbb_sum_stats/20002_1452.gwas.imputed_v3.both_sexes.tsv.bgz",
        "resources/ukbb_sum_stats/20002_1462.gwas.imputed_v3.both_sexes.tsv.bgz",
        "resources/ukbb_sum_stats/20002_1464.gwas.imputed_v3.both_sexes.tsv.bgz",
        "resources/ukbb_sum_stats/20002_1465.gwas.imputed_v3.both_sexes.tsv.bgz",
        "resources/ukbb_sum_stats/20002_1473.gwas.imputed_v3.both_sexes.tsv.bgz",
        "resources/ukbb_sum_stats/22126.gwas.imputed_v3.both_sexes.tsv.bgz",
        "resources/ukbb_sum_stats/6148_2.gwas.imputed_v3.both_sexes.tsv.bgz",
        "resources/ukbb_sum_stats/6148_5.gwas.imputed_v3.both_sexes.tsv.bgz",
        "resources/ukbb_sum_stats/D25.gwas.imputed_v3.both_sexes.tsv.bgz",
        "resources/ukbb_sum_stats/I9_IHD.gwas.imputed_v3.both_sexes.tsv.bgz",
        "resources/ukbb_sum_stats/K51.gwas.imputed_v3.both_sexes.tsv.bgz",
        "resources/ukbb_sum_stats/K57.gwas.imputed_v3.both_sexes.tsv.bgz",
        "resources/ukbb_sum_stats/K80.gwas.imputed_v3.both_sexes.tsv.bgz"
    output:
        "resources/ukbb_sum_stats/20002_1111.gwas.imputed_v3.both_sexes.tsv",
        "resources/ukbb_sum_stats/20002_1113.gwas.imputed_v3.both_sexes.tsv",
        "resources/ukbb_sum_stats/20002_1154.gwas.imputed_v3.both_sexes.tsv",
        "resources/ukbb_sum_stats/20002_1220.gwas.imputed_v3.both_sexes.tsv",
        "resources/ukbb_sum_stats/20002_1226.gwas.imputed_v3.both_sexes.tsv",
        "resources/ukbb_sum_stats/20002_1286.gwas.imputed_v3.both_sexes.tsv",
        "resources/ukbb_sum_stats/20002_1452.gwas.imputed_v3.both_sexes.tsv",
        "resources/ukbb_sum_stats/20002_1462.gwas.imputed_v3.both_sexes.tsv",
        "resources/ukbb_sum_stats/20002_1464.gwas.imputed_v3.both_sexes.tsv",
        "resources/ukbb_sum_stats/20002_1465.gwas.imputed_v3.both_sexes.tsv",
        "resources/ukbb_sum_stats/20002_1473.gwas.imputed_v3.both_sexes.tsv",
        "resources/ukbb_sum_stats/22126.gwas.imputed_v3.both_sexes.tsv",
        "resources/ukbb_sum_stats/6148_2.gwas.imputed_v3.both_sexes.tsv",
        "resources/ukbb_sum_stats/6148_5.gwas.imputed_v3.both_sexes.tsv",
        "resources/ukbb_sum_stats/D25.gwas.imputed_v3.both_sexes.tsv",
        "resources/ukbb_sum_stats/I9_IHD.gwas.imputed_v3.both_sexes.tsv",
        "resources/ukbb_sum_stats/K51.gwas.imputed_v3.both_sexes.tsv",
        "resources/ukbb_sum_stats/K57.gwas.imputed_v3.both_sexes.tsv",
        "resources/ukbb_sum_stats/K80.gwas.imputed_v3.both_sexes.tsv"
    run:
        for i,x in enumerate(input):
            shell("gunzip -c %s >%s" % (x, output[i]))

rule merge_ukbb_sum_stats:
    input:
        "resources/ukbb_sum_stats/20002_1111.gwas.imputed_v3.both_sexes.tsv",
        "resources/ukbb_sum_stats/20002_1113.gwas.imputed_v3.both_sexes.tsv",
        "resources/ukbb_sum_stats/20002_1154.gwas.imputed_v3.both_sexes.tsv",
        "resources/ukbb_sum_stats/20002_1220.gwas.imputed_v3.both_sexes.tsv",
        "resources/ukbb_sum_stats/20002_1226.gwas.imputed_v3.both_sexes.tsv",
        "resources/ukbb_sum_stats/20002_1286.gwas.imputed_v3.both_sexes.tsv",
        "resources/ukbb_sum_stats/20002_1452.gwas.imputed_v3.both_sexes.tsv",
        "resources/ukbb_sum_stats/20002_1462.gwas.imputed_v3.both_sexes.tsv",
        "resources/ukbb_sum_stats/20002_1464.gwas.imputed_v3.both_sexes.tsv",
        "resources/ukbb_sum_stats/20002_1465.gwas.imputed_v3.both_sexes.tsv",
        "resources/ukbb_sum_stats/20002_1473.gwas.imputed_v3.both_sexes.tsv",
        "resources/ukbb_sum_stats/22126.gwas.imputed_v3.both_sexes.tsv",
        "resources/ukbb_sum_stats/6148_2.gwas.imputed_v3.both_sexes.tsv",
        "resources/ukbb_sum_stats/6148_5.gwas.imputed_v3.both_sexes.tsv",
        "resources/ukbb_sum_stats/D25.gwas.imputed_v3.both_sexes.tsv",
        "resources/ukbb_sum_stats/I9_IHD.gwas.imputed_v3.both_sexes.tsv",
        "resources/ukbb_sum_stats/K51.gwas.imputed_v3.both_sexes.tsv",
        "resources/ukbb_sum_stats/K57.gwas.imputed_v3.both_sexes.tsv",
        "resources/ukbb_sum_stats/K80.gwas.imputed_v3.both_sexes.tsv"
    output:
        "resources/ukbb_sum_stats/merged_sum_stats.tsv"
    threads: 8
    shell:
        "Rscript scripts/merge_sum_stats.R"

rule make_plink_ranges:
    input:
      "resources/ukbb_sum_stats/merged_sum_stats.tsv",
      ("resources/1000g/euro/qc/chr%d_qc.bim" % x for x in range(1,23)),
      "resources/1000g/euro/qc/chrX_qc.bim"
    output:
      ("resources/ukbb_sum_stats/plink_ranges/chr%d.txt" % x for x in range(1,23)),
      "resources/ukbb_sum_stats/plink_ranges/chrX.txt"
    threads: 8
    shell:
      "Rscript scripts/make_plink_ranges.R"

rule subset_reference:
    input:
      "resources/1000g/euro/qc/{chr}_qc.bed",
      "resources/1000g/euro/qc/{chr}_qc.bim",
      "resources/1000g/euro/qc/{chr}_qc.fam",
      range_file = "resources/ukbb_sum_stats/plink_ranges/{chr}.txt"
    output:
      "resources/ukbb_sum_stats/plink_subsets/{chr}.bed",
      "resources/ukbb_sum_stats/plink_subsets/{chr}.bim",
      "resources/ukbb_sum_stats/plink_subsets/{chr}.fam"
    params:
      bfile = "resources/1000g/euro/qc/{chr}_qc",
      out = "resources/ukbb_sum_stats/plink_subsets/{chr}"
    threads: 8
    shell:
      "plink --memory 16000 --threads 8 --bfile {params.bfile} --extract {input.range_file} --make-bed --out {params.out}"

rule make_prune_ranges:
    input:
      "resources/ukbb_sum_stats/plink_subsets/{chr}.bed",
      "resources/ukbb_sum_stats/plink_subsets/{chr}.bim",
      "resources/ukbb_sum_stats/plink_subsets/{chr}.fam"
    output:
      "resources/ukbb_sum_stats/plink_subsets/pruned_ranges/{chr}.prune.in",
      "resources/ukbb_sum_stats/plink_subsets/pruned_ranges/{chr}.prune.out"
    params:
      bfile = "resources/ukbb_sum_stats/plink_subsets/{chr}",
      prune_out = "resources/ukbb_sum_stats/plink_subsets/pruned_ranges/{chr}"
    threads: 8
    shell:
      "plink --memory 16000 --threads 8 --bfile {params.bfile} --indep-pairwise 1000kb 50 0.2 --out {params.prune_out}"

rule cat_prune_ranges:
    input:
      ("resources/ukbb_sum_stats/plink_subsets/pruned_ranges/chr%d.prune.in" % x for x in range(1,23))
    output:
      "resources/ukbb_sum_stats/plink_subsets/pruned_ranges/all.prune.in"
    shell:
      "for x in {input}; do cat $x >>{output}; done"

rule cat_bim_files:
    input:
        ("resources/ukbb_sum_stats/plink_subsets/{chr}.bim" % x for x in range(1,23))
    output:
        "resources/ukbb_sum_stats/plink_subsets/all.bim"
    shell:
        "for x in {input}; do cat $x >>{output}; done"

rule prune_merged_sum_stats:
    input:
      "resources/ukbb_sum_stats/plink_subsets/all.bim",
      "resources/ukbb_sum_stats/plink_subsets/pruned_ranges/all.prune.in"
    output:
      "resources/ukbb_sum_stats/pruned_merged_sum_stats.tsv"
    shell:
      "Rscript prune_merged_sum_stats.R"
