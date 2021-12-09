import re
from itertools import chain

wildcard_constraints:
    chr = "chr[0-9XY]{1,2}"

def get_mem_mb(wildcards, threads):
    return threads * 3420

# TODO implement for production version
ukbb_trait_codes = ["20002_1220","20002_1289","20002_1286","20002_1291","I9_IHD","20002_1464","20002_1381","20002_1111","22126","20002_1462","K51","20002_1465","20002_1473","K57","K80","20002_1452","20002_1154","D25","6148_2","20002_1113","6148_5","20002_1226","E4_DM1"]

#ukbb_trait_pairs = list(chain(*[[ukbb_trait_codes[i]+'-'+ukbb_trait_codes[j] for j in range(i+1,len(ukbb_trait_codes))] for i in range(len(ukbb_trait_codes))]))

ukbb_trait_pairs = ["20002_1111-22126",
"20002_1473-I9_IHD",
"20002_1220-I9_IHD",
"20002_1111-20002_1113",
"20002_1220-20002_1465",
"20002_1220-20002_1473",
"20002_1220-K80",
"20002_1111-20002_1465",
"20002_1111-20002_1452",
"20002_1226-I9_IHD",
"20002_1465-I9_IHD",
"I9_IHD-K80",
"20002_1113-I9_IHD",
"20002_1465-K57",
"20002_1111-20002_1220",
"20002_1113-K57",
"20002_1286-20002_1465",
"20002_1154-20002_1286",
"20002_1113-20002_1220",
"20002_1286-K57",
"20002_1465-20002_1473",
"K57-K80",
"20002_1286-20002_1291",
"20002_1452-22126",
"20002_1154-20002_1465",
"20002_1113-20002_1465",
"20002_1220-K57",
"20002_1111-K57",
"20002_1226-20002_1473",
"20002_1220-20002_1226",
"20002_1113-20002_1473",
"20002_1111-20002_1154",
"20002_1154-K57",
"20002_1464-20002_1465",
"20002_1111-I9_IHD",
"20002_1154-I9_IHD",
"20002_1154-22126",
"20002_1113-20002_1286",
"20002_1154-20002_1220",
"20002_1226-20002_1465",
"20002_1465-K80",
"20002_1286-K80",
"20002_1473-K80",
"20002_1286-20002_1473",
"20002_1220-20002_1464",
"20002_1226-20002_1286",
"20002_1111-20002_1473",
"20002_1286-22126",
"20002_1462-K51",
"20002_1154-K80",
"20002_1111-20002_1286",
"20002_1465-22126",
"20002_1286-I9_IHD",
"20002_1220-20002_1286",
"I9_IHD-K57",
"20002_1111-K80",
"20002_1113-22126",
"20002_1154-20002_1226",
"20002_1226-K80",
"20002_1291-20002_1452",
"20002_1473-K57",
"20002_1462-I9_IHD",
"20002_1464-I9_IHD",
"20002_1154-20002_1473",
"22126-K57",
"20002_1113-K51",
"20002_1226-20002_1464",
"20002_1113-20002_1464",
"20002_1226-K57",
"20002_1226-D25",
"20002_1462-6148_2",
"20002_1464-20002_1473",
"20002_1111-20002_1464",
"20002_1286-D25",
"20002_1111-20002_1226",
"22126-6148_2",
"I9_IHD-K51",
"20002_1220-6148_2",
"20002_1113-K80",
"20002_1464-K51",
"20002_1113-20002_1226",
"6148_2-K51",
"20002_1452-20002_1473",
"20002_1462-20002_1473",
"20002_1452-I9_IHD",
"20002_1291-K51",
"22126-K51",
"6148_2-I9_IHD",
"20002_1462-22126",
"20002_1464-6148_2",
"20002_1111-20002_1462",
"20002_1464-K80",
"20002_1291-20002_1462",
"20002_1111-D25",
"20002_1113-20002_1154",
"20002_1452-6148_2",
"20002_1452-20002_1465",
"20002_1154-20002_1452",
"20002_1113-6148_5",
"20002_1286-20002_1452",
"20002_1464-22126",
"20002_1465-6148_2",
"20002_1291-K57",
"20002_1220-20002_1291",
"6148_5-K57",
"20002_1291-20002_1473",
"D25-K57",
"20002_1473-D25",
"20002_1465-D25",
"20002_1286-6148_5",
"20002_1464-K57",
"20002_1220-K51",
"20002_1291-K80",
"20002_1286-20002_1462",
"20002_1291-6148_2",
"20002_1291-22126",
"20002_1291-20002_1465",
"20002_1452-K80",
"K51-K57",
"20002_1291-6148_5",
"20002_1462-6148_5",
"20002_1226-K51",
"20002_1111-6148_2",
"20002_1464-D25",
"22126-D25",
"6148_2-D25",
"20002_1113-20002_1462",
"20002_1462-K57",
"20002_1473-6148_5",
"20002_1226-20002_1452",
"20002_1220-20002_1462",
"20002_1111-K51",
"20002_1226-6148_5",
"6148_2-6148_5",
"20002_1291-D25",
"20002_1226-20002_1291",
"20002_1286-20002_1464",
"20002_1462-D25",
"20002_1462-K80",
"D25-I9_IHD",
"20002_1111-20002_1291",
"20002_1113-D25",
"20002_1154-20002_1291",
"6148_5-D25",
"20002_1154-K51",
"20002_1291-I9_IHD",
"20002_1462-20002_1465",
"20002_1452-K51",
"20002_1111-6148_5",
"20002_1452-6148_5",
"20002_1226-6148_2",
"20002_1291-20002_1464",
"6148_2-K80",
"6148_5-K51",
"20002_1154-6148_5",
"20002_1154-D25",
"22126-K80",
"20002_1220-D25",
"20002_1462-20002_1464",
"D25-K80",
"20002_1113-20002_1452",
"20002_1154-6148_2",
"22126-I9_IHD",
"20002_1113-20002_1291",
"20002_1113-6148_2",
"20002_1220-6148_5",
"D25-K51",
"20002_1220-20002_1452",
"6148_2-K57",
"22126-6148_5",
"20002_1464-6148_5",
"20002_1154-20002_1464",
"20002_1473-6148_2",
"20002_1465-K51",
"20002_1452-D25",
"20002_1226-22126",
"20002_1473-22126",
"20002_1226-20002_1462",
"6148_5-I9_IHD",
"K51-K80",
"20002_1286-6148_2",
"20002_1473-K51",
"6148_5-K80",
"20002_1154-20002_1462",
"20002_1465-6148_5",
"20002_1452-20002_1462",
"20002_1220-22126",
"20002_1286-K51",
"20002_1452-20002_1464",
"20002_1452-K57",
"20002_1111-20002_1289",
"20002_1111-20002_1381",
"20002_1113-20002_1289",
"20002_1113-20002_1381",
"20002_1154-20002_1289",
"20002_1154-20002_1381",
"20002_1220-20002_1289",
"20002_1220-20002_1381",
"20002_1226-20002_1289",
"20002_1226-20002_1381",
"20002_1286-20002_1289",
"20002_1286-20002_1381",
"20002_1289-20002_1291",
"20002_1289-20002_1381",
"20002_1289-20002_1452",
"20002_1289-20002_1462",
"20002_1289-20002_1464",
"20002_1289-20002_1465",
"20002_1289-20002_1473",
"20002_1289-22126",
"20002_1289-6148_2",
"20002_1289-6148_5",
"20002_1289-D25",
"20002_1289-I9_IHD",
"20002_1289-K51",
"20002_1289-K57",
"20002_1289-K80",
"20002_1291-20002_1381",
"20002_1381-20002_1452",
"20002_1381-20002_1462",
"20002_1381-20002_1464",
"20002_1381-20002_1465",
"20002_1381-20002_1473",
"20002_1381-22126",
"20002_1381-6148_2",
"20002_1381-6148_5",
"20002_1381-D25",
"20002_1381-I9_IHD",
"20002_1381-K51",
"20002_1381-K57",
"20002_1381-K80",
"20002_1220-E4_DM1",
"20002_1289-E4_DM1",
"20002_1286-E4_DM1",
"20002_1291-E4_DM1",
"I9_IHD-E4_DM1",
"20002_1464-E4_DM1",
"20002_1381-E4_DM1",
"20002_1111-E4_DM1",
"22126-E4_DM1",
"20002_1462-E4_DM1",
"K51-E4_DM1",
"20002_1465-E4_DM1",
"20002_1473-E4_DM1",
"K57-E4_DM1",
"K80-E4_DM1",
"20002_1452-E4_DM1",
"20002_1154-E4_DM1",
"D25-E4_DM1",
"6148_2-E4_DM1",
"20002_1113-E4_DM1",
"6148_5-E4_DM1",
"20002_1226-E4_DM1"]

pid_ukbb_trait_pairs = ["pid-20002_1220",
"pid-20002_1286",
"pid-20002_1289",
"pid-20002_1291",
"pid-20002_1381",
"pid-I9_IHD",
"pid-20002_1464",
"pid-20002_1111",
"pid-22126",
"pid-20002_1462",
"pid-K51",
"pid-20002_1465",
"pid-20002_1473",
"pid-K57",
"pid-K80",
"pid-20002_1452",
"pid-20002_1154",
"pid-D25",
"pid-6148_2",
"pid-20002_1113",
"pid-6148_5",
"pid-20002_1226",
"pid-E4_DM1"]

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

rule merge_pid_and_ukbb_sum_stats:
    input:
        ukbb_files = ["resources/ukbb_sum_stats/%s.gwas.imputed_v3.both_sexes.tsv" % x for x in ukbb_trait_codes],
        pid_file = "resources/pid.tsv.gz"
    output:
        merged_file = "resources/merged_pid_ukbb_sum_stats.tsv.gz",
        temp_files = ["resources/%s.temp" % x for x in glob_wildcards("resources/ukbb_sum_stats/{id}.gwas.imputed_v3.both_sexes.tsv").id]+["resources/pid.temp"]
    threads: 8
    shell:
      """
      Rscript scripts/merge_pid_and_ukbb_sum_stats.R --ukbb_files {input.ukbb_files} --pid_file {input.pid_file} -o {output.merged_file} -nt {threads};
      for x in {output.temp_files}; do touch $x; done
      """

# SNP sets differ between UKBB pairs and PID-UKBB pairs
rule make_ukbb_plink_ranges:
    input:
      sum_stats_file = ancient("resources/merged_pid_ukbb_sum_stats.tsv.gz"),
      bim_files = ancient([("resources/1000g/euro/qc/chr%d_qc.bim" % x for x in range(1,23)),
      "resources/1000g/euro/qc/chrX_qc.bim"])
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
      sum_stats_file = ancient("resources/merged_pid_ukbb_sum_stats.tsv.gz"),
      bim_files = ancient([("resources/1000g/euro/qc/chr%d_qc.bim" % x for x in range(1,23)), "resources/1000g/euro/qc/chrX_qc.bim"])
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
      sum_stats_file = ancient("resources/merged_pid_ukbb_sum_stats.tsv.gz"),
      bim_file = ancient("resources/plink_subsets/{join}/all.bim"),
      pruned_range_file = ancient("resources/plink_ranges/{join}/pruned_ranges/all.prune.in")
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
      ancient("resources/{trait_A}.temp"),
      ancient("resources/{trait_B}.temp"),
      sum_stats_file = ancient("resources/pruned_sum_stats/{join}/pruned_merged_sum_stats.tsv"),
    output:
      temp("results/{join}/{trait_A}-{trait_B}_gps_value.tsv")
    shell:
      "scripts/gps_cpp/build/apps/computeGpsCLI -i {input.sum_stats_file} -a {wildcards.trait_A} -b {wildcards.trait_B} -c {wildcards.trait_A} -d {wildcards.trait_B} -o {output}"

rule permute_trait_pair:
    input:
      ancient("resources/{trait_A}.temp"),
      ancient("resources/{trait_B}.temp"),
      sum_stats_file = ancient("resources/pruned_sum_stats/{join}/pruned_merged_sum_stats.tsv"),
    output:
      "results/{join}/permutations/{trait_A}-{trait_B}.tsv"
    threads: 8
    resources:
        mem_mb=get_mem_mb
    shell:
      "scripts/gps_cpp/build/apps/permuteTraitsCLI -i {input.sum_stats_file} -o {output} -a {wildcards.trait_A} -b {wildcards.trait_B} -c 8 -n 375"

rule compute_gps_pvalue_for_trait_pair:
    input:
      gps_file = "results/{join}/{trait_A}-{trait_B}_gps_value.tsv",
      perm_file = "results/{join}/permutations/{trait_A}-{trait_B}.tsv"
    output:
      "results/{join}/{trait_A}-{trait_B}_gps_pvalue.tsv"
    shell:
      "Rscript scripts/compute_gps_pvalue.R -g {input.gps_file} -p {input.perm_file} -a {wildcards.trait_A} -b {wildcards.trait_B} -o {output}"

rule collate_gps_pvalue_data:
    input:
        ["results/ukbb/%s_gps_pvalue.tsv" % x for x in ukbb_trait_pairs]+
        ["results/pid_ukbb/%s_gps_pvalue.tsv" % x for x in pid_ukbb_trait_pairs]
    output:
        "results/combined_pvalues.tsv"
    run:
        with open(output[0], 'w') as outfile:
            outfile.write(("\t".join(["trait_A", "trait_B", "gps", "n", "loc", "loc.sd", "scale", "scale.sd", "shape", "shape.sd", "pval"]))+"\n")
            for i,x in enumerate(input):
                with open(x, 'r') as infile:
                    data_line = infile.readlines()[1]

                m = re.match("results/(pid_){0,1}ukbb/(\w+)-(\w+)_gps_pvalue.tsv", x)

                outfile.write(("\t".join([m[2], m[3], data_line])))

rule add_trait_labels_to_gps_results:
    input:
      pvalue_file = "results/combined_pvalues.tsv",
      lookup_file = "resources/ukbb_sum_stats/trait_metadata.tsv"
    output:
      "results/combined_pvalues_with_labels.tsv"
    shell:
      "Rscript scripts/add_trait_labels_to_gps_results.R -p {input.pvalue_file} -l {input.lookup_file} -o {output}"

rule plot_gps_heatmaps:
    input:
        gps_file = "results/combined_pvalues_with_labels.tsv"
    params:
        # traits argument allows ordering of traits in matrix
        traits = "bipolar disorder,glaucoma,MD,depression,schizophrenia,leiomyoma,emphysema/chronic bronchitis,hypercholesterolaemia,osteoarthritis,IHD,diverticulosis,cholelithiasis,IBS,CD,UC,hayfever,asthma,eczema/derm,rheumatoid arthritis,lupus,hypothyroidism,T1D"
    output:
        "results/plots/gps_pval_heatmap.png",
    shell:
        "Rscript scripts/gps_heatmaps.R -i {input.gps_file} -t '{params.traits}' -o {output}"

rule plot_exemplar_rg_gps_heatmap:
    input:
      pvalue_file = "results/combined_pvalues_with_labels.tsv",
      rg_file = "resources/exemplar_rg_values.tsv",
    params:
      # traits argument allows ordering of traits in matrix
      traits = "bipolar disorder,glaucoma,MD,depression,schizophrenia,leiomyoma,emphysema/chronic bronchitis,hypercholesterolaemia,osteoarthritis,IHD,diverticulosis,cholelithiasis,IBS,CD,UC,hayfever,asthma,eczema/derm,rheumatoid arthritis,lupus,hypothyroidism"
    output:
        rg_pvalue_file = "results/plots/exemplar_rg_pval_vs_gps.png",
        rg_estimate_file = "results/plots/exemplar_rg_estimate_vs_gps.png"
    shell:
        "Rscript scripts/rg_vs_gps_heatmap.R -p {input.pvalue_file} -r {input.rg_file} -t '{params.traits}' --rg_pvalue_output {output.rg_pvalue_file} --rg_estimate_output {output.rg_estimate_file}"

rule plot_ukbb_rg_gps_heatmap:
    input:
        pvalue_file = "results/combined_pvalues_with_labels.tsv",
        rg_file = "resources/ukbb_rg_values.tsv",
    params:
        # traits argument allows ordering of traits in matrix
        traits = "bipolar disorder,glaucoma,MD,depression,schizophrenia,leiomyoma,emphysema/chronic bronchitis,hypercholesterolaemia,osteoarthritis,IHD,diverticulosis,cholelithiasis,IBS,CD,UC,hayfever,asthma,eczema/derm,rheumatoid arthritis,lupus,hypothyroidism"
    output:
        rg_pvalue_file = "results/plots/ukbb_rg_pval_vs_gps.png",
        rg_estimate_file = "results/plots/ukbb_rg_estimate_vs_gps.png"
    shell:
        "Rscript scripts/rg_vs_gps_heatmap.R -p {input.pvalue_file} -r {input.rg_file} -t '{params.traits}' --rg_pvalue_output {output.rg_pvalue_file} --rg_estimate_output {output.rg_estimate_file}"

rule plot_rg_heatmaps:
    input:
        rg_file = "resources/rg_values.tsv"
    params:
        # traits argument allows ordering of traits in matrix
        traits = "bipolar disorder,glaucoma,MD,depression,schizophrenia,leiomyoma,emphysema/chronic bronchitis,hypercholesterolaemia,osteoarthritis,IHD,diverticulosis,cholelithiasis,IBS,CD,UC,hayfever,asthma,eczema/derm,rheumatoid arthritis,lupus,hypothyroidism"
    output:
        rg_pvalue_file = "results/plots/rg_pval.png",
        rg_estimate_file = "results/plots/rg_estimate.png",
        rg_pvalue_estimate_file = "results/plots/rg_pvalue_estimate.png"
    shell:
        "Rscript scripts/rg_heatmaps.R -r {input.rg_file} -t '{params.traits}' --rg_pvalue_output {output.rg_pvalue_file} --rg_estimate_output {output.rg_estimate_file} --rg_pvalue_estimate_output {output.rg_pvalue_estimate_file}"

rule compute_hoeffdings_for_trait_pair:
    input:
        sum_stats_file = ancient("resources/pruned_sum_stats/{join}/pruned_merged_sum_stats.tsv"),
    output:
        "results/{join}/{trait_A}-{trait_B}_hoeffdings.tsv"
    shell:
        "Rscript scripts/compute_hoeffdings.R -i {input.sum_stats_file} -a {wildcards.trait_A} -b {wildcards.trait_B} -o {output} -nt 1"

rule collate_hoeffdings_results:
    input:
      ["results/ukbb/%s_hoeffdings.tsv" % x for x in ukbb_trait_pairs]+
      ["results/pid_ukbb/%s_hoeffdings.tsv" % x for x in pid_ukbb_trait_pairs]
    output:
      "results/hoeffdings_results.tsv"
    shell:
       """
       echo -e 'trait_A\ttrait_B\tn\tDn\tscaled\tp.value' >> {output}
       for x in {input}; do
          cat <(tail -n 1 $x) >> {output}
       done
       """

rule add_trait_labels_to_hoeffdings_results:
    input:
        results_file = "results/hoeffdings_results.tsv",
        lookup_file = "resources/ukbb_sum_stats/trait_metadata.tsv"
    output:
        "results/hoeffdings_results_with_labels.tsv"
    shell:
        "Rscript scripts/add_trait_labels_to_hoeffdings_results.R -r {input.results_file} -l {input.lookup_file} -o {output}"

rule plot_exemplar_rg_hoeffdings_heatmaps:
    input:
        hoeffdings_file = "results/hoeffdings_results_with_labels.tsv",
        rg_file = "resources/exemplar_rg_values.tsv",
    params:
        # traits argument allows me to set order of traits in matrix
        traits = "bipolar disorder,glaucoma,MD,depression,schizophrenia,leiomyoma,emphysema/chronic bronchitis,hypercholesterolaemia,osteoarthritis,IHD,diverticulosis,cholelithiasis,IBS,CD,UC,hayfever,asthma,eczema/derm,rheumatoid arthritis,lupus,hypothyroidism"
    output:
        rg_pvalue_file = "results/plots/exemplar_rg_pval_vs_hoeffdings.png",
        rg_estimate_file = "results/plots/exemplar_rg_estimate_vs_hoeffdings.png"
    shell:
        "Rscript scripts/rg_vs_hoeffdings_heatmap.R -f {input.hoeffdings_file} -r {input.rg_file} -t '{params.traits}' --rg_pvalue_output {output.rg_pvalue_file} --rg_estimate_output {output.rg_estimate_file}"

rule compute_correlation_measures_for_trait_pair:
    input:
        ancient("resources/{trait_A}.temp"),
        ancient("resources/{trait_B}.temp"),
        sum_stats_file = ancient("resources/pruned_sum_stats/{join}/pruned_merged_sum_stats.tsv"),
    output:
        "results/{join}/{trait_A}-{trait_B}_correlation_measures.tsv"
    shell:
        "Rscript scripts/compute_correlation_measures.R -i {input.sum_stats_file} -a {wildcards.trait_A} -b {wildcards.trait_B} -o {output}"

rule collate_correlation_measures:
    input:
        ["results/ukbb/%s_correlation_measures.tsv" % x for x in ukbb_trait_pairs]+
        ["results/pid_ukbb/%s_correlation_measures.tsv" % x for x in pid_ukbb_trait_pairs]
    output:
        "results/correlation_measures.tsv"
    shell:
        """
        echo -e 'trait_A\ttrait_B\tn\tpearson\tspearman' >> {output}
        for x in {input}; do
        cat <(tail -n 1 $x) >> {output}
        done
        """

rule add_labels_to_collated_correlation_measures:
    input:
      corr_file = "results/correlation_measures.tsv",
      lookup_file = "resources/ukbb_sum_stats/trait_metadata.tsv"
    output:
      "results/correlation_measures_with_labels.tsv"
    shell:
      "Rscript scripts/add_trait_labels_to_correlation_measures.R -p {input.corr_file} -l {input.lookup_file} -o {output}"

rule plot_null_dists_to_compare:
    input:
      perm_file = "results/{join}/permutations/{trait_A}-{trait_B}.tsv",
      fit_file = "results/{join}/{trait_A}-{trait_B}_gps_pvalue.tsv"
    output:
      exp1_null = "results/plots/{join}/{trait_A}-{trait_B}_exp1_null_dist.png",
      gev_null = "results/plots/{join}/{trait_A}-{trait_B}_gev_null_dist.png",
      exp1_gev_combined = "results/plots/{join}/{trait_A}-{trait_B}_null_dists.png"
    shell:
      "Rscript scripts/plot_null_dists_to_compare.R -f {input.fit_file} -p {input.perm_file} --exp1_null {output.exp1_null} --gev_null {output.gev_null} --exp1_gev_combined {output.exp1_gev_combined}"

rule plot_all_null_dists_to_compare:
    input:
        ["results/plots/ukbb/%s_null_dists.png" % x for x in ukbb_trait_pairs]+
        ["results/plots/pid_ukbb/%s_null_dists.png" % x for x in pid_ukbb_trait_pairs]

rule plot_gof_plots:
    input:
      perm_file = "results/{join}/permutations/{trait_A}-{trait_B}.tsv",
      fit_file = "results/{join}/{trait_A}-{trait_B}_gps_pvalue.tsv"
    output:
      "results/plots/{join}/{trait_A}-{trait_B}_gof_plots.png"
    shell:
      "Rscript scripts/plot_gev_gof_plots.R -f {input.fit_file} -p {input.perm_file} -l {wildcards.trait_A}-{wildcards.trait_B} -o {output}"

rule plot_all_gof_plots:
    input:
        ["results/plots/ukbb/%s_gof_plots.png" % x for x in ukbb_trait_pairs]+
        ["results/plots/pid_ukbb/%s_gof_plots.png" % x for x in pid_ukbb_trait_pairs]

# TODO remove me
rule write_freq_map:
    input:
        ancient("resources/{trait}.temp"),
        sum_stats_file = ancient("resources/pruned_sum_stats/{join}/pruned_merged_sum_stats.tsv"),
    output:
        "results/{join}/freq_map/add_epsilon/{trait}.tsv"
    shell:
        "scripts/gps_cpp/build/apps/writeFreqMapCLI -i {input.sum_stats_file} -a {wildcards.trait} -b {output}"

# TODO remove me
rule freq_maps:
    input:
        ["results/ukbb/freq_map/add_epsilon/%s.tsv" % x for x in ukbb_trait_codes]
