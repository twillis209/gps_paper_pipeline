wildcard_constraints:
    chr = "chr[0-9XY]{1,2}"

def get_mem_mb(wildcards, threads):
    return threads * 3420

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

rule download_ukbb_sum_stats:
    output:
        "resources/ukbb_sum_stats/{id}.gwas.imputed_v3.both_sexes.tsv.bgz",
    shell:
        """
        wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/{wildcards.id}.gwas.imputed_v3.both_sexes.tsv.bgz -O resources/ukbb_sum_stats/{wildcards.id}.gwas.imputed_v3.both_sexes.tsv.bgz
        """

rule decompress_ukbb_sum_stats:
    input:
        "resources/ukbb_sum_stats/{id}.gwas.imputed_v3.both_sexes.tsv.bgz"
    output:
        "resources/ukbb_sum_stats/{id}.gwas.imputed_v3.both_sexes.tsv"
    run:
        for i,x in enumerate(input):
            shell("gunzip -c %s >%s" % (x, output[i]))

# TODO add temp files to later rules
rule merge_pid_and_ukbb_sum_stats:
    input:
        ukbb_files = ["resources/ukbb_sum_stats/%s.gwas.imputed_v3.both_sexes.tsv" % x for x in glob_wildcards("resources/ukbb_sum_stats/{id}.gwas.imputed_v3.both_sexes.tsv").id],
        pid_file = "resources/pid.tsv.gz"
    output:
        merged_file = "resources/merged_pid_ukbb_sum_stats.tsv.gz",
        temp_files = [temp("resources/%s.temp" % x) for x in glob_wildcards("resources/ukbb_sum_stats/{id}.gwas.imputed_v3.both_sexes.tsv").id]+["resources/pid.temp"]
    threads: 8
    shell:
      """
      Rscript scripts/merge_pid_and_ukbb_sum_stats.R --ukbb_files {input.ukbb_files} --pid_file {input.pid_file} -o {output.merged_file} -nt {threads};
      for x in {output.temp_files}; do touch $x; done
      """

# SNP sets differ between UKBB pairs and PID-UKBB pairs
rule make_ukbb_plink_ranges:
    input:
      sum_stats_file = "resources/merged_pid_ukbb_sum_stats.tsv.gz",
      bim_files = [("resources/1000g/euro/qc/chr%d_qc.bim" % x for x in range(1,23)),
      "resources/1000g/euro/qc/chrX_qc.bim"]
    output:
      bim_files = [("resources/plink_ranges/ukbb/chr%d.txt" % x for x in range(1,23)),
      "resources/plink_ranges/ukbb/chrX.txt"]
    params:
      non_na_cols = "pval.%s" % glob_wildcards("resources/ukbb_sum_stats/{id}.gwas.imputed_v3.both_sexes.tsv").id[0]
    threads: 8
    shell:
      "Rscript scripts/make_plink_ranges.R --sum_stats_file {input.sum_stats_file} --input_bim_files {input.bim_files} --non_na_cols {params.non_na_cols} --output_bim_files {output.bim_files} -nt {threads}"

rule make_pid_ukbb_plink_ranges:
    input:
      sum_stats_file = "resources/merged_pid_ukbb_sum_stats.tsv.gz",
      bim_files = [("resources/1000g/euro/qc/chr%d_qc.bim" % x for x in range(1,23)), "resources/1000g/euro/qc/chrX_qc.bim"]
    output:
      bim_files = [("resources/plink_ranges/pid_ukbb/chr%d.txt" % x for x in range(1,23)), "resources/plink_ranges/pid_ukbb/chrX.txt"]
    params:
      non_na_cols = ["pval.pid", "pval.%s" % glob_wildcards("resources/ukbb_sum_stats/{id}.gwas.imputed_v3.both_sexes.tsv").id[0]]
    threads: 8
    shell:
        "Rscript scripts/make_plink_ranges.R --sum_stats_file {input.sum_stats_file} --input_bim_files {input.bim_files} --non_na_cols {params.non_na_cols} --output_bim_files {output.bim_files} -nt {threads}"

rule subset_reference:
    input:
      "resources/1000g/euro/qc/{chr}_qc.bed",
      "resources/1000g/euro/qc/{chr}_qc.bim",
      "resources/1000g/euro/qc/{chr}_qc.fam",
      range_file = "resources/plink_ranges/{join}/{chr}.txt"
    output:
      "resources/plink_subsets/{join}/{chr}.bed",
      "resources/plink_subsets/{join}/{chr}.bim",
      "resources/plink_subsets/{join}/{chr}.fam"
    params:
      bfile = "resources/1000g/euro/qc/{chr}_qc",
      out = "resources/plink_subsets/{join}/{chr}"
    threads: 8
    resources:
        mem_mb=get_mem_mb
    shell:
      "plink --memory {resources.mem_mb} --threads {threads} --bfile {params.bfile} --extract {input.range_file} --make-bed --out {params.out}"

rule make_pruned_ranges:
    input:
      "resources/plink_subsets/{join}/{chr}.bed",
      "resources/plink_subsets/{join}/{chr}.bim",
      "resources/plink_subsets/{join}/{chr}.fam"
    output:
      "resources/plink_ranges/{join}/pruned_ranges/{chr}.prune.in",
      "resources/plink_ranges/{join}/pruned_ranges/{chr}.prune.out"
    params:
      bfile = "resources/plink_subsets/{join}/{chr}",
      prune_out = "resources/plink_ranges/{join}/pruned_ranges/{chr}"
    threads: 8
    resources:
        mem_mb=get_mem_mb
    shell:
      "plink --memory {resources.mem_mb} --threads {threads} --bfile {params.bfile} --indep-pairwise 1000kb 50 0.2 --out {params.prune_out}"

rule cat_pruned_ranges:
    input:
      ("resources/plink_ranges/{join}/pruned_ranges/chr%d.prune.in" % x for x in range(1,23))
    output:
      "resources/plink_ranges/{join}/pruned_ranges/all.prune.in"
    shell:
      "for x in {input}; do cat $x >>{output}; done"

rule cat_bim_files:
    input:
        ("resources/plink_subsets/{join}/chr%d.bim" % x for x in range(1,23))
    output:
        "resources/plink_subsets/{join}/all.bim"
    shell:
        "for x in {input}; do cat $x >>{output}; done"

rule prune_merged_sum_stats:
    input:
      sum_stats_file = "resources/merged_pid_ukbb_sum_stats.tsv.gz",
      bim_file = "resources/plink_subsets/{join}/all.bim",
      pruned_range_file = "resources/plink_ranges/{join}/pruned_ranges/all.prune.in"
    output:
      "resources/pruned_sum_stats/{join}/pruned_merged_sum_stats.tsv"
    threads: 8
    shell:
      """
      Rscript scripts/prune_merged_sum_stats.R --sum_stats_file {input.sum_stats_file} --bim_file {input.bim_file} --range_file {input.pruned_range_file} -o {output} -nt {threads};
      sed -i 's/pval\.//g' {output}
      """

rule compute_gps_for_trait_pair:
    input:
      "resources/{trait_A}.temp",
      "resources/{trait_B}.temp",
      sum_stats_file = "resources/pruned_sum_stats/{join}/pruned_merged_sum_stats.tsv",
    output:
      temp("results/{join}/{trait_A}-{trait_B}_gps_value.tsv")
    shell:
      "scripts/gps_cpp/build/apps/computeGpsForTraitPairCLI -i {input.sum_stats_file} -a {wildcards.trait_A} -b {wildcards.trait_B} -o {output}"

rule permute_trait_pair:
    input:
      "resources/{trait_A}.temp",
      "resources/{trait_B}.temp",
      sum_stats_file = "resources/pruned_sum_stats/{join}/pruned_merged_sum_stats.tsv",
    output:
      "results/{join}/permutations/{trait_A}-{trait_B}.tsv"
    threads: 8
    resources:
        mem_mb=get_mem_mb
    shell:
      "scripts/gps_cpp/build/apps/permuteTraitsCLI -i {input.sum_stats_file} -o {output} -a {wildcards.trait_A} -b {wildcards.trait_B} -c 8 -n 375"

rule compute_gps_p_value_for_trait_pair:
    input:
      gps_file = "results/{join}/{trait_A}-{trait_B}_gps_value.tsv",
      perm_file = "results/{join}/permutations/{trait_A}-{trait_B}.tsv"
    output:
      "results/{join}/{trait_A}-{trait_B}_gps_pvalue.tsv"
    shell:
      "Rscript scripts/compute_gps_pvalue.R -g {input.gps_file} -p {input.perm_file} -a {wildcards.trait_A} -b {wildcards.trait_B} -o {output}"

# TODO make CLI argument
rule collate_gps_p_value_data:
    input:
#        "resources/ukbb_sum_stats/traits_codes_abbrv_cases.csv",
#        "resources/ukbb_sum_stats/traits_rg.tsv",
        "results/ukbb/20002_1111-pid_gps_pvalue.tsv"
        "results/ukbb/20002_1111-22126_gps_pvalue.tsv",
        "results/ukbb/20002_1473-I9_IHD_gps_pvalue.tsv",
        "results/ukbb/20002_1220-I9_IHD_gps_pvalue.tsv",
        "results/ukbb/20002_1111-20002_1113_gps_pvalue.tsv",
        "results/ukbb/20002_1220-20002_1465_gps_pvalue.tsv",
        "results/ukbb/20002_1220-20002_1473_gps_pvalue.tsv",
        "results/ukbb/20002_1220-K80_gps_pvalue.tsv",
        "results/ukbb/20002_1111-20002_1465_gps_pvalue.tsv",
        "results/ukbb/20002_1111-20002_1452_gps_pvalue.tsv",
        "results/ukbb/20002_1226-I9_IHD_gps_pvalue.tsv",
        "results/ukbb/20002_1465-I9_IHD_gps_pvalue.tsv",
        "results/ukbb/I9_IHD-K80_gps_pvalue.tsv",
        "results/ukbb/20002_1113-I9_IHD_gps_pvalue.tsv",
        "results/ukbb/20002_1465-K57_gps_pvalue.tsv",
        "results/ukbb/20002_1111-20002_1220_gps_pvalue.tsv",
        "results/ukbb/20002_1113-K57_gps_pvalue.tsv",
        "results/ukbb/20002_1286-20002_1465_gps_pvalue.tsv",
        "results/ukbb/20002_1154-20002_1286_gps_pvalue.tsv",
        "results/ukbb/20002_1113-20002_1220_gps_pvalue.tsv",
        "results/ukbb/20002_1286-K57_gps_pvalue.tsv",
        "results/ukbb/20002_1465-20002_1473_gps_pvalue.tsv",
        "results/ukbb/K57-K80_gps_pvalue.tsv",
        "results/ukbb/20002_1286-20002_1291_gps_pvalue.tsv",
        "results/ukbb/20002_1452-22126_gps_pvalue.tsv",
        "results/ukbb/20002_1154-20002_1465_gps_pvalue.tsv",
        "results/ukbb/20002_1113-20002_1465_gps_pvalue.tsv",
        "results/ukbb/20002_1220-K57_gps_pvalue.tsv",
        "results/ukbb/20002_1111-K57_gps_pvalue.tsv",
        "results/ukbb/20002_1226-20002_1473_gps_pvalue.tsv",
        "results/ukbb/20002_1220-20002_1226_gps_pvalue.tsv",
        "results/ukbb/20002_1113-20002_1473_gps_pvalue.tsv",
        "results/ukbb/20002_1111-20002_1154_gps_pvalue.tsv",
        "results/ukbb/20002_1154-K57_gps_pvalue.tsv",
        "results/ukbb/20002_1464-20002_1465_gps_pvalue.tsv",
        "results/ukbb/20002_1111-I9_IHD_gps_pvalue.tsv",
        "results/ukbb/20002_1154-I9_IHD_gps_pvalue.tsv",
        "results/ukbb/20002_1154-22126_gps_pvalue.tsv",
        "results/ukbb/20002_1113-20002_1286_gps_pvalue.tsv",
        "results/ukbb/20002_1154-20002_1220_gps_pvalue.tsv",
        "results/ukbb/20002_1226-20002_1465_gps_pvalue.tsv",
        "results/ukbb/20002_1465-K80_gps_pvalue.tsv",
        "results/ukbb/20002_1286-K80_gps_pvalue.tsv",
        "results/ukbb/20002_1473-K80_gps_pvalue.tsv",
        "results/ukbb/20002_1286-20002_1473_gps_pvalue.tsv",
        "results/ukbb/20002_1220-20002_1464_gps_pvalue.tsv",
        "results/ukbb/20002_1226-20002_1286_gps_pvalue.tsv",
        "results/ukbb/20002_1111-20002_1473_gps_pvalue.tsv",
        "results/ukbb/20002_1286-22126_gps_pvalue.tsv",
        "results/ukbb/20002_1462-K51_gps_pvalue.tsv",
        "results/ukbb/20002_1154-K80_gps_pvalue.tsv",
        "results/ukbb/20002_1111-20002_1286_gps_pvalue.tsv",
        "results/ukbb/20002_1465-22126_gps_pvalue.tsv",
        "results/ukbb/20002_1286-I9_IHD_gps_pvalue.tsv",
        "results/ukbb/20002_1220-20002_1286_gps_pvalue.tsv",
        "results/ukbb/I9_IHD-K57_gps_pvalue.tsv",
        "results/ukbb/20002_1111-K80_gps_pvalue.tsv",
        "results/ukbb/20002_1113-22126_gps_pvalue.tsv",
        "results/ukbb/20002_1154-20002_1226_gps_pvalue.tsv",
        "results/ukbb/20002_1226-K80_gps_pvalue.tsv",
        "results/ukbb/20002_1291-20002_1452_gps_pvalue.tsv",
        "results/ukbb/20002_1473-K57_gps_pvalue.tsv",
        "results/ukbb/20002_1462-I9_IHD_gps_pvalue.tsv",
        "results/ukbb/20002_1464-I9_IHD_gps_pvalue.tsv",
        "results/ukbb/20002_1154-20002_1473_gps_pvalue.tsv",
        "results/ukbb/22126-K57_gps_pvalue.tsv",
        "results/ukbb/20002_1113-K51_gps_pvalue.tsv",
        "results/ukbb/20002_1226-20002_1464_gps_pvalue.tsv",
        "results/ukbb/20002_1113-20002_1464_gps_pvalue.tsv",
        "results/ukbb/20002_1226-K57_gps_pvalue.tsv",
        "results/ukbb/20002_1226-D25_gps_pvalue.tsv",
        "results/ukbb/20002_1462-6148_2_gps_pvalue.tsv",
        "results/ukbb/20002_1464-20002_1473_gps_pvalue.tsv",
        "results/ukbb/20002_1111-20002_1464_gps_pvalue.tsv",
        "results/ukbb/20002_1286-D25_gps_pvalue.tsv",
        "results/ukbb/20002_1111-20002_1226_gps_pvalue.tsv",
        "results/ukbb/22126-6148_2_gps_pvalue.tsv",
        "results/ukbb/I9_IHD-K51_gps_pvalue.tsv",
        "results/ukbb/20002_1220-6148_2_gps_pvalue.tsv",
        "results/ukbb/20002_1113-K80_gps_pvalue.tsv",
        "results/ukbb/20002_1464-K51_gps_pvalue.tsv",
        "results/ukbb/20002_1113-20002_1226_gps_pvalue.tsv",
        "results/ukbb/6148_2-K51_gps_pvalue.tsv",
        "results/ukbb/20002_1452-20002_1473_gps_pvalue.tsv",
        "results/ukbb/20002_1462-20002_1473_gps_pvalue.tsv",
        "results/ukbb/20002_1452-I9_IHD_gps_pvalue.tsv",
        "results/ukbb/20002_1291-K51_gps_pvalue.tsv",
        "results/ukbb/22126-K51_gps_pvalue.tsv",
        "results/ukbb/6148_2-I9_IHD_gps_pvalue.tsv",
        "results/ukbb/20002_1462-22126_gps_pvalue.tsv",
        "results/ukbb/20002_1464-6148_2_gps_pvalue.tsv",
        "results/ukbb/20002_1111-20002_1462_gps_pvalue.tsv",
        "results/ukbb/20002_1464-K80_gps_pvalue.tsv",
        "results/ukbb/20002_1291-20002_1462_gps_pvalue.tsv",
        "results/ukbb/20002_1111-D25_gps_pvalue.tsv",
        "results/ukbb/20002_1113-20002_1154_gps_pvalue.tsv",
        "results/ukbb/20002_1452-6148_2_gps_pvalue.tsv",
        "results/ukbb/20002_1452-20002_1465_gps_pvalue.tsv",
        "results/ukbb/20002_1154-20002_1452_gps_pvalue.tsv",
        "results/ukbb/20002_1113-6148_5_gps_pvalue.tsv",
        "results/ukbb/20002_1286-20002_1452_gps_pvalue.tsv",
        "results/ukbb/20002_1464-22126_gps_pvalue.tsv",
        "results/ukbb/20002_1465-6148_2_gps_pvalue.tsv",
        "results/ukbb/20002_1291-K57_gps_pvalue.tsv",
        "results/ukbb/20002_1220-20002_1291_gps_pvalue.tsv",
        "results/ukbb/6148_5-K57_gps_pvalue.tsv",
        "results/ukbb/20002_1291-20002_1473_gps_pvalue.tsv",
        "results/ukbb/D25-K57_gps_pvalue.tsv",
        "results/ukbb/20002_1473-D25_gps_pvalue.tsv",
        "results/ukbb/20002_1465-D25_gps_pvalue.tsv",
        "results/ukbb/20002_1286-6148_5_gps_pvalue.tsv",
        "results/ukbb/20002_1464-K57_gps_pvalue.tsv",
        "results/ukbb/20002_1220-K51_gps_pvalue.tsv",
        "results/ukbb/20002_1291-K80_gps_pvalue.tsv",
        "results/ukbb/20002_1286-20002_1462_gps_pvalue.tsv",
        "results/ukbb/20002_1291-6148_2_gps_pvalue.tsv",
        "results/ukbb/20002_1291-22126_gps_pvalue.tsv",
        "results/ukbb/20002_1291-20002_1465_gps_pvalue.tsv",
        "results/ukbb/20002_1452-K80_gps_pvalue.tsv",
        "results/ukbb/K51-K57_gps_pvalue.tsv",
        "results/ukbb/20002_1291-6148_5_gps_pvalue.tsv",
        "results/ukbb/20002_1462-6148_5_gps_pvalue.tsv",
        "results/ukbb/20002_1226-K51_gps_pvalue.tsv",
        "results/ukbb/20002_1111-6148_2_gps_pvalue.tsv",
        "results/ukbb/20002_1464-D25_gps_pvalue.tsv",
        "results/ukbb/22126-D25_gps_pvalue.tsv",
        "results/ukbb/6148_2-D25_gps_pvalue.tsv",
        "results/ukbb/20002_1113-20002_1462_gps_pvalue.tsv",
        "results/ukbb/20002_1462-K57_gps_pvalue.tsv",
        "results/ukbb/20002_1473-6148_5_gps_pvalue.tsv",
        "results/ukbb/20002_1226-20002_1452_gps_pvalue.tsv",
        "results/ukbb/20002_1220-20002_1462_gps_pvalue.tsv",
        "results/ukbb/20002_1111-K51_gps_pvalue.tsv",
        "results/ukbb/20002_1226-6148_5_gps_pvalue.tsv",
        "results/ukbb/6148_2-6148_5_gps_pvalue.tsv",
        "results/ukbb/20002_1291-D25_gps_pvalue.tsv",
        "results/ukbb/20002_1226-20002_1291_gps_pvalue.tsv",
        "results/ukbb/20002_1286-20002_1464_gps_pvalue.tsv",
        "results/ukbb/20002_1462-D25_gps_pvalue.tsv",
        "results/ukbb/20002_1462-K80_gps_pvalue.tsv",
        "results/ukbb/D25-I9_IHD_gps_pvalue.tsv",
        "results/ukbb/20002_1111-20002_1291_gps_pvalue.tsv",
        "results/ukbb/20002_1113-D25_gps_pvalue.tsv",
        "results/ukbb/20002_1154-20002_1291_gps_pvalue.tsv",
        "results/ukbb/6148_5-D25_gps_pvalue.tsv",
        "results/ukbb/20002_1154-K51_gps_pvalue.tsv",
        "results/ukbb/20002_1291-I9_IHD_gps_pvalue.tsv",
        "results/ukbb/20002_1462-20002_1465_gps_pvalue.tsv",
        "results/ukbb/20002_1452-K51_gps_pvalue.tsv",
        "results/ukbb/20002_1111-6148_5_gps_pvalue.tsv",
        "results/ukbb/20002_1452-6148_5_gps_pvalue.tsv",
        "results/ukbb/20002_1226-6148_2_gps_pvalue.tsv",
        "results/ukbb/20002_1291-20002_1464_gps_pvalue.tsv",
        "results/ukbb/6148_2-K80_gps_pvalue.tsv",
        "results/ukbb/6148_5-K51_gps_pvalue.tsv",
        "results/ukbb/20002_1154-6148_5_gps_pvalue.tsv",
        "results/ukbb/20002_1154-D25_gps_pvalue.tsv",
        "results/ukbb/22126-K80_gps_pvalue.tsv",
        "results/ukbb/20002_1220-D25_gps_pvalue.tsv",
        "results/ukbb/20002_1462-20002_1464_gps_pvalue.tsv",
        "results/ukbb/D25-K80_gps_pvalue.tsv",
        "results/ukbb/20002_1113-20002_1452_gps_pvalue.tsv",
        "results/ukbb/20002_1154-6148_2_gps_pvalue.tsv",
        "results/ukbb/22126-I9_IHD_gps_pvalue.tsv",
        "results/ukbb/20002_1113-20002_1291_gps_pvalue.tsv",
        "results/ukbb/20002_1113-6148_2_gps_pvalue.tsv",
        "results/ukbb/20002_1220-6148_5_gps_pvalue.tsv",
        "results/ukbb/D25-K51_gps_pvalue.tsv",
        "results/ukbb/20002_1220-20002_1452_gps_pvalue.tsv",
        "results/ukbb/6148_2-K57_gps_pvalue.tsv",
        "results/ukbb/22126-6148_5_gps_pvalue.tsv",
        "results/ukbb/20002_1464-6148_5_gps_pvalue.tsv",
        "results/ukbb/20002_1154-20002_1464_gps_pvalue.tsv",
        "results/ukbb/20002_1473-6148_2_gps_pvalue.tsv",
        "results/ukbb/20002_1465-K51_gps_pvalue.tsv",
        "results/ukbb/20002_1452-D25_gps_pvalue.tsv",
        "results/ukbb/20002_1226-22126_gps_pvalue.tsv",
        "results/ukbb/20002_1473-22126_gps_pvalue.tsv",
        "results/ukbb/20002_1226-20002_1462_gps_pvalue.tsv",
        "results/ukbb/6148_5-I9_IHD_gps_pvalue.tsv",
        "results/ukbb/K51-K80_gps_pvalue.tsv",
        "results/ukbb/20002_1286-6148_2_gps_pvalue.tsv",
        "results/ukbb/20002_1473-K51_gps_pvalue.tsv",
        "results/ukbb/6148_5-K80_gps_pvalue.tsv",
        "results/ukbb/20002_1154-20002_1462_gps_pvalue.tsv",
        "results/ukbb/20002_1465-6148_5_gps_pvalue.tsv",
        "results/ukbb/20002_1452-20002_1462_gps_pvalue.tsv",
        "results/ukbb/20002_1220-22126_gps_pvalue.tsv",
        "results/ukbb/20002_1286-K51_gps_pvalue.tsv",
        "results/ukbb/20002_1452-20002_1464_gps_pvalue.tsv",
        "results/ukbb/20002_1452-K57_gps_pvalue.tsv",
        "results/ukbb/20002_1111-20002_1289_gps_pvalue.tsv",
        "results/ukbb/20002_1111-20002_1381_gps_pvalue.tsv",
        "results/ukbb/20002_1113-20002_1289_gps_pvalue.tsv",
        "results/ukbb/20002_1113-20002_1381_gps_pvalue.tsv",
        "results/ukbb/20002_1154-20002_1289_gps_pvalue.tsv",
        "results/ukbb/20002_1154-20002_1381_gps_pvalue.tsv",
        "results/ukbb/20002_1220-20002_1289_gps_pvalue.tsv",
        "results/ukbb/20002_1220-20002_1381_gps_pvalue.tsv",
        "results/ukbb/20002_1226-20002_1289_gps_pvalue.tsv",
        "results/ukbb/20002_1226-20002_1381_gps_pvalue.tsv",
        "results/ukbb/20002_1286-20002_1289_gps_pvalue.tsv",
        "results/ukbb/20002_1286-20002_1381_gps_pvalue.tsv",
        "results/ukbb/20002_1289-20002_1291_gps_pvalue.tsv",
        "results/ukbb/20002_1289-20002_1381_gps_pvalue.tsv",
        "results/ukbb/20002_1289-20002_1452_gps_pvalue.tsv",
        "results/ukbb/20002_1289-20002_1462_gps_pvalue.tsv",
        "results/ukbb/20002_1289-20002_1464_gps_pvalue.tsv",
        "results/ukbb/20002_1289-20002_1465_gps_pvalue.tsv",
        "results/ukbb/20002_1289-20002_1473_gps_pvalue.tsv",
        "results/ukbb/20002_1289-22126_gps_pvalue.tsv",
        "results/ukbb/20002_1289-6148_2_gps_pvalue.tsv",
        "results/ukbb/20002_1289-6148_5_gps_pvalue.tsv",
        "results/ukbb/20002_1289-D25_gps_pvalue.tsv",
        "results/ukbb/20002_1289-I9_IHD_gps_pvalue.tsv",
        "results/ukbb/20002_1289-K51_gps_pvalue.tsv",
        "results/ukbb/20002_1289-K57_gps_pvalue.tsv",
        "results/ukbb/20002_1289-K80_gps_pvalue.tsv",
        "results/ukbb/20002_1291-20002_1381_gps_pvalue.tsv",
        "results/ukbb/20002_1381-20002_1452_gps_pvalue.tsv",
        "results/ukbb/20002_1381-20002_1462_gps_pvalue.tsv",
        "results/ukbb/20002_1381-20002_1464_gps_pvalue.tsv",
        "results/ukbb/20002_1381-20002_1465_gps_pvalue.tsv",
        "results/ukbb/20002_1381-20002_1473_gps_pvalue.tsv",
        "results/ukbb/20002_1381-22126_gps_pvalue.tsv",
        "results/ukbb/20002_1381-6148_2_gps_pvalue.tsv",
        "results/ukbb/20002_1381-6148_5_gps_pvalue.tsv",
        "results/ukbb/20002_1381-D25_gps_pvalue.tsv",
        "results/ukbb/20002_1381-I9_IHD_gps_pvalue.tsv",
        "results/ukbb/20002_1381-K51_gps_pvalue.tsv",
        "results/ukbb/20002_1381-K57_gps_pvalue.tsv",
        "results/ukbb/20002_1381-K80_gps_pvalue.tsv"
#    output:
#        "results/collated_gps_pvalues.tsv"
#    shell:
#        "Rscript scripts/collate_gps_pvalues.R"

