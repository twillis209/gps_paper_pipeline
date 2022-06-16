ukbb_trait_codes = ["20002_1220","20002_1289","20002_1286","20002_1291","I9_IHD","20002_1464","20002_1381","20002_1111","22126","20002_1462","K51","20002_1465","20002_1473","K57","K80","20002_1452","20002_1154","D25","6148_2","20002_1113","6148_5","20002_1226","E4_DM1","I42","N80"]

# TODO implement for production version
ukbb_trait_pairs_new = list(chain(*[[ukbb_trait_codes[i]+'-'+ukbb_trait_codes[j] for j in range(i+1,len(ukbb_trait_codes))] for i in range(len(ukbb_trait_codes))]))

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
"20002_1226-E4_DM1",
"20002_1220-I42",
"20002_1289-I42",
"20002_1286-I42",
"20002_1291-I42",
"I9_IHD-I42",
"20002_1464-I42",
"20002_1381-I42",
"20002_1111-I42",
"22126-I42",
"20002_1462-I42",
"K51-I42",
"20002_1465-I42",
"20002_1473-I42",
"K57-I42",
"K80-I42",
"20002_1452-I42",
"20002_1154-I42",
"D25-I42",
"6148_2-I42",
"20002_1113-I42",
"6148_5-I42",
"20002_1226-I42",
"E4_DM1-I42",
"20002_1220-N80",
"20002_1289-N80",
"20002_1286-N80",
"20002_1291-N80",
"I9_IHD-N80",
"20002_1464-N80",
"20002_1381-N80",
"20002_1111-N80",
"22126-N80",
"20002_1462-N80",
"K51-N80",
"20002_1465-N80",
"20002_1473-N80",
"K57-N80",
"K80-N80",
"20002_1452-N80",
"20002_1154-N80",
"D25-N80",
"6148_2-N80",
"20002_1113-N80",
"6148_5-N80",
"20002_1226-N80",
"E4_DM1-N80",
"I42-N80"]

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

trait_pairs_for_increasing_perm_fits = ["20002_1111-20002_1113", "20002_1473-I9_IHD", "20002_1220-K80", "20002_1465-K57", "20002_1452-22126", "20002_1465-K80",  "20002_1286-D25", "6148_2-K51", "20002_1291-K80"]

rule download_ukbb_sum_stats:
    output:
        temp("resources/ukbb_sum_stats/{id}.gwas.imputed_v3.both_sexes.tsv.bgz"),
    shell:
        """
        wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/{wildcards.id}.gwas.imputed_v3.both_sexes.tsv.bgz -O resources/ukbb_sum_stats/{wildcards.id}.gwas.imputed_v3.both_sexes.tsv.bgz
        """

rule decompress_ukbb_sum_stats:
    input:
        "resources/ukbb_sum_stats/{id}.gwas.imputed_v3.both_sexes.tsv.bgz"
    output:
        temp("resources/ukbb_sum_stats/{id}.gwas.imputed_v3.both_sexes.tsv")
    run:
        for i,x in enumerate(input):
            shell("gunzip -c %s >%s" % (x, output[i]))

rule merge_ukbb_sum_stats:
    input:
        ukbb_files = ["resources/ukbb_sum_stats/%s.gwas.imputed_v3.both_sexes.tsv" % x for x in ukbb_trait_codes]
    output:
        merged_file = "resources/ukbb_sum_stats/merged_ukbb_sum_stats.tsv.gz"
    params:
        ukbb_trait_codes = ukbb_trait_codes
    threads: 8
    script: "merge_ukbb_sum_stats.R"

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
      Rscript workflow/scripts/merge_pid_and_ukbb_sum_stats.R --ukbb_files {input.ukbb_files} --pid_file {input.pid_file} -o {output.merged_file} -nt {threads};
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
      non_na_cols = [f"pval.{x}" for x in ukbb_trait_codes]
    threads: 8
    shell:
      "Rscript workflow/scripts/make_plink_ranges.R --sum_stats_file {input.sum_stats_file} --input_bim_files {input.bim_files} --non_na_cols {params.non_na_cols} --output_bim_files {output.bim_files} -nt {threads}"

rule make_pid_ukbb_plink_ranges:
    input:
      sum_stats_file = ancient("resources/merged_pid_ukbb_sum_stats.tsv.gz"),
      bim_files = ancient([("resources/1000g/euro/qc/chr%d_qc.bim" % x for x in range(1,23)), "resources/1000g/euro/qc/chrX_qc.bim"])
    output:
      bim_files = [("resources/plink_ranges/pid_ukbb/chr%d.txt" % x for x in range(1,23)), "resources/plink_ranges/pid_ukbb/chrX.txt"]
    params:
        non_na_cols = ['pval.pid']+[f"pval.{x}" for x in ukbb_trait_codes]
    threads: 8
    shell:
        "Rscript workflow/scripts/make_plink_ranges.R --sum_stats_file {input.sum_stats_file} --input_bim_files {input.bim_files} --non_na_cols {params.non_na_cols} --output_bim_files {output.bim_files} -nt {threads}"

rule prune_merged_sum_stats:
    input:
      sum_stats_file = ancient("resources/merged_pid_ukbb_sum_stats.tsv.gz"),
      bim_file = ancient("resources/plink_subsets/{join}/all.bim"),
      pruned_range_file = ancient("resources/plink_ranges/{join}/pruned_ranges/window_{window}_step_{step}/all.prune.in")
    output:
        "resources/pruned_sum_stats/{join}/all_pruned_snps/window_{window}_step_{step}/pruned_merged_sum_stats.tsv"
    threads: 8
    shell:
      """
      Rscript workflow/scripts/prune_merged_sum_stats.R --sum_stats_file {input.sum_stats_file} --bim_file {input.bim_file} --range_file {input.pruned_range_file} -o {output} -nt {threads};
      sed -i 's/pval\.//g' {output}
      """

rule downsample_pruned_merged_sum_stats:
    input:
        ancient("resources/pruned_sum_stats/{join}/all_pruned_snps/window_{window}_step_{step}/pruned_merged_sum_stats.tsv")
    output:
        temp("resources/pruned_sum_stats/{join}/{no_snps}_snps/window_{window}_step_{step}/pruned_merged_sum_stats.tsv")
    threads: 8
    shell:
     "Rscript workflow/scripts/downsample_sum_stats.R -i {input} -n {wildcards.no_snps} -o {output} -nt {threads}"

rule excise_mhc_from_pruned_merged_sum_stats:
    input:
        ancient("resources/pruned_sum_stats/{join}/all_pruned_snps/window_{window}_step_{step}/pruned_merged_sum_stats.tsv")
    output:
        "resources/pruned_sum_stats/{join}/sans_mhc/window_{window}_step_{step}/pruned_merged_sum_stats.tsv"
    threads: 8
    shell:
        "Rscript workflow/scripts/excise_mhc_from_sum_stats.R -i {input} -o {output} -nt {threads}"
