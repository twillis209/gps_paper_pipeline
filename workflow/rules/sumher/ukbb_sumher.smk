from scipy.stats import chi2

rule thin_predictors_for_ukbb:
    input:
        "resources/ukbb_sum_stats/{join}/nodup/snps_only/{chr}.bed",
        "resources/ukbb_sum_stats/{join}/nodup/snps_only/{chr}.bim",
        "resources/ukbb_sum_stats/{join}/nodup/snps_only/{chr}.fam"
    output:
        thin_file = "results/ldak/ldak-thin/weights/ukbb/{join}/{chr}/thin.in",
        weights_file = "results/ldak/ldak-thin/weights/ukbb/{join}/{chr}/weights.thin"
    log:
        log_file = "results/ldak/ldak-thin/weights/ukbb/{join}/{chr}/thin.log"
    params:
        input_stem = "resources/ukbb_sum_stats/{join}/nodup/snps_only/{chr}",
        output_stem = "results/ldak/ldak-thin/weights/ukbb/{join}/{chr}/thin"
    group: "ldsc_hoeffding_sumher_gps_sans_permutation"
    shell:
        """
        $ldakRoot/ldak --thin {params.output_stem} --bfile {params.input_stem} --window-prune .98 --window-kb 100 > {log.log_file};
        awk < {output.thin_file} '{{print $1, 1}}' > {output.weights_file}
        """

rule calculate_ldak_thin_taggings_for_chromosome_for_ukbb:
    input:
        "resources/ukbb_sum_stats/{join}/nodup/snps_only/{chr}.bed",
        "resources/ukbb_sum_stats/{join}/nodup/snps_only/{chr}.bim",
        "resources/ukbb_sum_stats/{join}/nodup/snps_only/{chr}.fam",
        weights_file = "results/ldak/ldak-thin/weights/ukbb/{join}/{chr}/weights.thin"
    output:
        tagging_file = temp("results/ldak/ldak-thin/taggings/ukbb/{join}/{chr}.tagging"),
    log:
        log_file = "results/ldak/ldak-thin/taggings/ukbb/{join}/{chr}.tagging.log"
    params:
        input_stem = "resources/ukbb_sum_stats/{join}/nodup/snps_only/{chr}",
        output_stem = "results/ldak/ldak-thin/taggings/ukbb/{join}/{chr}"
    group: "ldsc_hoeffding_sumher_gps_sans_permutation"
    shell:
        "$ldakRoot/ldak --calc-tagging {params.output_stem} --bfile {params.input_stem} --weights {input.weights_file} --chr {wildcards.chr} --window-kb 1000 --power -.25 > {log.log_file}"

rule join_ldak_thin_taggings_for_ukbb:
    input:
        [f"results/ldak/ldak-thin/taggings/ukbb/{{join}}/chr{x}.tagging" for x in range(1, 23)]
    output:
        wg_tagging_file = "results/ldak/ldak-thin/ukbb/{join}/whole_genome.tagging",
        chrom_taggings_file = temp("results/ldak/ldak-thin/ukbb/{join}/taggings.txt")
    log:
        log_file = "results/ldak/ldak-thin/ukbb/{join}/whole_genome.tagging.log"
    params:
        output_stem = "results/ldak/ldak-thin/ukbb/{join}/whole_genome"
    group: "ldsc_hoeffding_sumher_gps_sans_permutation"
    shell:
        """
        for x in {input}; do
            echo $x >> {output.chrom_taggings_file}
        done;

        $ldakRoot/ldak --join-tagging {params.output_stem} --taglist {output.chrom_taggings_file} > {log.log_file}
        """

rule process_ukbb_sum_stats:
    input:
        ancient("resources/ukbb_sum_stats/{join}/merged_ukbb_sum_stats.tsv.gz")
    output:
        "resources/ukbb_sum_stats/{join}/{trait_code}.assoc"
    params:
        z_colname = lambda wildcards: '.'.join(['tstat', wildcards.trait_code]),
        sample_size_colname = lambda wildcards: '.'.join(['n_complete_samples', wildcards.trait_code])
    threads: 8
    resources:
        runtime = 10
    group: "ldsc_hoeffding_sumher_gps_sans_permutation"
    script:
        "../../scripts/process_ukbb_sum_stats.R"

rule estimate_rg_with_ldak_thin_for_ukbb:
    input:
        wg_tagging_file = "results/ldak/ldak-thin/ukbb/{join}/whole_genome.tagging",
        sum_stats_file_A = "resources/ukbb_sum_stats/{join}/{trait_A}.assoc",
        sum_stats_file_B = "resources/ukbb_sum_stats/{join}/{trait_B}.assoc"
    output:
        progress_file = "results/ldak/ldak-thin/ukbb/{join}/rg/{trait_A}-{trait_B}.progress",
        cors_file = "results/ldak/ldak-thin/ukbb/{join}/rg/{trait_A}-{trait_B}.cors",
        cors_full_file = "results/ldak/ldak-thin/ukbb/{join}/rg/{trait_A}-{trait_B}.cors.full",
        labels_file = "results/ldak/ldak-thin/ukbb/{join}/rg/{trait_A}-{trait_B}.cors.labels",
        overlap_file = "results/ldak/ldak-thin/ukbb/{join}/rg/{trait_A}-{trait_B}.cors.overlap"
    log:
        log_file = "results/ldak/ldak-thin/ukbb/{join}/rg/{trait_A}-{trait_B}.log"
    params:
        output_stem = "results/ldak/ldak-thin/ukbb/{join}/rg/{trait_A}-{trait_B}"
    resources:
        runtime = 5
    group: "ldsc_hoeffding_sumher_gps_sans_permutation"
    shell:
        """
        $ldakRoot/ldak --sum-cors {params.output_stem} --tagfile {input.wg_tagging_file} --summary {input.sum_stats_file_A} --summary2 {input.sum_stats_file_B} --allow-ambiguous YES --check-sums NO --cutoff 0.01 > {log.log_file}
        """

rule ukbb_sumher:
    input:
       sumher_files = [f"results/ldak/ldak-thin/ukbb/{{join}}/rg/{trait_pair}.cors.full" for trait_pair in ukbb_trait_pairs],
       metadata_file = "resources/ukbb_sum_stats/trait_metadata.tsv"
    output:
        "results/ldak/ldak-thin/ukbb/{join}/rg/compiled_ukbb_sumher_results.tsv"
    run:
        meta_daf = pd.read_csv(input.metadata_file, sep = '\t')
        with open(output[0], 'w') as outfile:
            outfile.write("trait_A_code\ttrait_B_code\ttrait_A_name\ttrait_B_name\th2.A\th2.A.se\th2.B\th2.B.se\tgcov\tgcov.se\trg\trg.se\trg.z\trg.p\n")
            for x in input.sumher_files:
                head, tail = os.path.split(x)

                trait_A, trait_B = re.match("(\w+)-(\w+).cors.full", tail).groups()

                trait_A_name = meta_daf.loc[meta_daf.code == trait_A, 'long_abbrv'].values[0]
                trait_B_name = meta_daf.loc[meta_daf.code == trait_B, 'long_abbrv'].values[0]

                with open(x, 'r') as infile:
                    line = infile.readline()
                    line = infile.readline()

                    # Category Trait1_Her SD Trait2_Her SD Both_Coher SD Correlation SD
                    _, h2_A, h2_A_se, h2_B, h2_B_se, cov, cov_se, rg, rg_se = line.split()
                    rg_z = float(rg)/float(rg_se)

                    rg_p = chi2.sf(rg_z**2, df = 1, loc = 0, scale = 1)
                    outfile.write(f"{trait_A}\t{trait_B}\t{trait_A_name}\t{trait_B_name}\t{float(h2_A)}\t{float(h2_A_se)}\t{float(h2_B)}\t{float(h2_B_se)}\t{float(cov)}\t{float(cov_se)}\t{float(rg)}\t{float(rg_se)}\t{rg_z}\t{rg_p}\n")
