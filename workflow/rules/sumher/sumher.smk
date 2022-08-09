from scipy.stats import chi2

rule thin_predictors:
    input:
        "resources/1000g/euro/qc/nodup/snps_only/{join}/{chr}.bed",
        "resources/1000g/euro/qc/nodup/snps_only/{join}/{chr}.bim",
        "resources/1000g/euro/qc/nodup/snps_only/{join}/{chr}.fam"
    output:
        thin_file = "results/ldak/ldak-thin/weights/{join}/{chr}/thin.in",
        weights_file = "results/ldak/ldak-thin/weights/{join}/{chr}/weights.thin"
    log:
        log_file = "results/ldak/ldak-thin/weights/{join}/{chr}/thin.log"
    params:
        input_stem = "resources/1000g/euro/qc/nodup/snps_only/{join}/{chr}",
        output_stem = "results/ldak/ldak-thin/weights/{join}/{chr}/thin"
    group: "ldsc_hoeffding_sumher_gps_sans_permutation"
    shell:
        """
        $ldakRoot/ldak --thin {params.output_stem} --bfile {params.input_stem} --window-prune .98 --window-kb 100 > {log.log_file};
        awk < {output.thin_file} '{{print $1, 1}}' > {output.weights_file}
        """

rule calculate_ldak_thin_taggings_for_chromosome:
    input:
        "resources/1000g/euro/qc/nodup/snps_only/{join}/{chr}.bed",
        "resources/1000g/euro/qc/nodup/snps_only/{join}/{chr}.bim",
        "resources/1000g/euro/qc/nodup/snps_only/{join}/{chr}.fam",
        weights_file = "results/ldak/ldak-thin/weights/{join}/{chr}/weights.thin"
    output:
        tagging_file = temp("results/ldak/ldak-thin/taggings/{join}/{chr}.tagging"),
    log:
        log_file = "results/ldak/ldak-thin/taggings/{join}/{chr}.tagging.log"
    params:
        input_stem = "resources/1000g/euro/qc/nodup/snps_only/{join}/{chr}",
        output_stem = "results/ldak/ldak-thin/taggings/{join}/{chr}"
    group: "ldsc_hoeffding_sumher_gps_sans_permutation"
    shell:
        "$ldakRoot/ldak --calc-tagging {params.output_stem} --bfile {params.input_stem} --weights {input.weights_file} --chr {wildcards.chr} --window-kb 1000 --power -.25 > {log.log_file}"

rule join_ldak_thin_taggings:
    input:
        [f"results/ldak/ldak-thin/taggings/{{join}}/chr{x}.tagging" for x in range(1, 23)]
    output:
        wg_tagging_file = "results/ldak/ldak-thin/{join}/whole_genome.tagging",
        chrom_taggings_file = temp("results/ldak/ldak-thin/{join}/taggings.txt")
    log:
        log_file = "results/ldak/ldak-thin/{join}/whole_genome.tagging.log"
    params:
        output_stem = "results/ldak/ldak-thin/{join}/whole_genome"
    group: "ldsc_hoeffding_sumher_gps_sans_permutation"
    shell:
        """
        for x in {input}; do
            echo $x >> {output.chrom_taggings_file}
        done;

        $ldakRoot/ldak --join-tagging {params.output_stem} --taglist {output.chrom_taggings_file} > {log.log_file}
        """

# TODO: Currently assuming in the script that ncases_A == ncases_B and ncontrols_A == ncontrols_B
rule process_combined_simgwas_sum_stats:
    input:
        "results/simgwas/simulated_sum_stats/whole_genome_sum_stats/{no_reps}_reps/randomised/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/seed_{seed}_sum_stats_{pair_label}_tag_{tag}_of_{tag_A}-{tag_B}.tsv.gz"
    output:
        temp("results/simgwas/simulated_sum_stats/whole_genome_sum_stats/{no_reps}_reps/randomised/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/seed_{seed}_sum_stats_{pair_label}_tag_{tag}_of_{tag_A}-{tag_B}.assoc")
    params:
        z_colname = lambda wildcards: f'zsim.{wildcards.tag}',
        chr_colname = 'chr',
        bp_colname = 'position',
        a1_colname = 'a0',
        a2_colname = 'a1',
        sample_size = lambda wildcards: int(wildcards.ncases_A)+int(wildcards.ncontrols_A)
    threads: 2
    resources:
        runtime = 10
    group: "ldsc_hoeffding_sumher_gps_sans_permutation"
    script:
        "../../process_combined_simgwas_sum_stats_for_sumher.R"

# NB: Currently assuming in the script that ncases_A == ncases_B and ncontrols_A == ncontrols_B
rule process_ukbb_sum_stats:
    input:
        ancient("resources/ukbb_sum_stats/merged_ukbb_sum_stats.tsv.gz")
    output:
        "resources/ukbb_sum_stats/{trait_code}.assoc"
    params:
        z_colname = lambda wildcards: '.'.join(['tstat', wildcards.trait_code]),
        sample_size_colname = lambda wildcards: '.'.join(['n_complete_samples', wildcards.trait_code])
    threads: 8
    resources:
        runtime = 10
    group: "ldsc_hoeffding_sumher_gps_sans_permutation"
    script:
        "../../process_ukbb_sum_stats.R"

rule estimate_rg_with_ldak_thin_for_simgwas:
    input:
        wg_tagging_file = "results/ldak/ldak-thin/1000g/whole_genome.tagging",
        sum_stats_file_A = "results/simgwas/simulated_sum_stats/whole_genome_sum_stats/{no_reps}_reps/randomised/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/seed_{seed}_sum_stats_A_tag_{tag_A}_of_{tag_A}-{tag_B}.assoc",
        sum_stats_file_B = "results/simgwas/simulated_sum_stats/whole_genome_sum_stats/{no_reps}_reps/randomised/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/seed_{seed}_sum_stats_B_tag_{tag_B}_of_{tag_A}-{tag_B}.assoc"
    output:
        progress_file = "results/ldak/ldak-thin/simgwas/{no_reps}_reps/randomised/rg/{ncases_A,\d+}_{ncontrols_A,\d+}_{ncases_B,\d+}_{ncontrols_B,\d+}/{effect_blocks_A,[smlvh\d-]+}_{effect_blocks_B,[smlvh\d-]+}_{shared_effect_blocks,[smlvh\d-]+}/seed_{seed,\d+}_tags_{tag_A,\d+}-{tag_B,\d+}.progress",
        cors_file = "results/ldak/ldak-thin/simgwas/{no_reps}_reps/randomised/rg/{ncases_A,\d+}_{ncontrols_A,\d+}_{ncases_B,\d+}_{ncontrols_B,\d+}/{effect_blocks_A,[smlvh\d-]+}_{effect_blocks_B,[smlvh\d-]+}_{shared_effect_blocks,[smlvh\d-]+}/seed_{seed,\d+}_tags_{tag_A,\d+}-{tag_B,\d+}.cors",
        cors_full_file = "results/ldak/ldak-thin/simgwas/{no_reps}_reps/randomised/rg/{ncases_A,\d+}_{ncontrols_A,\d+}_{ncases_B,\d+}_{ncontrols_B,\d+}/{effect_blocks_A,[smlvh\d-]+}_{effect_blocks_B,[smlvh\d-]+}_{shared_effect_blocks,[smlvh\d-]+}/seed_{seed,\d+}_tags_{tag_A,\d+}-{tag_B,\d+}.cors.full",
        labels_file = "results/ldak/ldak-thin/simgwas/{no_reps}_reps/randomised/rg/{ncases_A,\d+}_{ncontrols_A,\d+}_{ncases_B,\d+}_{ncontrols_B,\d+}/{effect_blocks_A,[smlvh\d-]+}_{effect_blocks_B,[smlvh\d-]+}_{shared_effect_blocks,[smlvh\d-]+}/seed_{seed,\d+}_tags_{tag_A,\d+}-{tag_B,\d+}.labels",
        overlap_file = "results/ldak/ldak-thin/simgwas/{no_reps}_reps/randomised/rg/{ncases_A,\d+}_{ncontrols_A,\d+}_{ncases_B,\d+}_{ncontrols_B,\d+}/{effect_blocks_A,[smlvh\d-]+}_{effect_blocks_B,[smlvh\d-]+}_{shared_effect_blocks,[smlvh\d-]+}/seed_{seed,\d+}_tags_{tag_A,\d+}-{tag_B,\d+}.overlap"
    log:
        log_file = "results/ldak/ldak-thin/simgwas/{no_reps}_reps/randomised/rg/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/seed_{seed}_tags_{tag_A}-{tag_B}.log"
    params:
        output_stem = "results/ldak/ldak-thin/simgwas/{no_reps}_reps/randomised/rg/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/seed_{seed}_tags_{tag_A}-{tag_B}"
    resources:
        runtime = 5
    group: "ldsc_hoeffding_sumher_gps_sans_permutation"
    shell:
        """
        $ldakRoot/ldak --sum-cors {params.output_stem} --tagfile {input.wg_tagging_file} --summary {input.sum_stats_file_A} --summary2 {input.sum_stats_file_B} --allow-ambiguous YES --check-sums NO --cutoff 0.01 > {log.log_file}
        """

rule estimate_rg_with_ldak_thin_for_ukbb:
    input:
        wg_tagging_file = "results/ldak/ldak-thin/ukbb/whole_genome.tagging",
        sum_stats_file_A = "resources/ukbb_sum_stats/{trait_A}.assoc",
        sum_stats_file_B = "resources/ukbb_sum_stats/{trait_B}.assoc"
    output:
        progress_file = "results/ldak/ldak-thin/ukbb/rg/{trait_A}-{trait_B}.progress",
        cors_file = "results/ldak/ldak-thin/ukbb/rg/{trait_A}-{trait_B}.cors",
        cors_full_file = "results/ldak/ldak-thin/ukbb/rg/{trait_A}-{trait_B}.cors.full",
#        labels_file = "results/ldak/ldak-thin/ukbb/rg/{trait_A}-{trait_B}.cors.labels",
#        overlap_file = "results/ldak/ldak-thin/ukbb/rg/{trait_A}-{trait_B}.cors.overlap"
    log:
        log_file = "results/ldak/ldak-thin/ukbb/rg/{trait_A}-{trait_B}.log"
    params:
        output_stem = "results/ldak/ldak-thin/ukbb/rg/{trait_A}-{trait_B}"
    resources:
        runtime = 5
    group: "ldsc_hoeffding_sumher_gps_sans_permutation"
    shell:
        """
        $ldakRoot/ldak --sum-cors {params.output_stem} --tagfile {input.wg_tagging_file} --summary {input.sum_stats_file_A} --summary2 {input.sum_stats_file_B} --allow-ambiguous YES --check-sums NO --cutoff 0.01 > {log.log_file}
        """

rule ukbb_sumher:
    input:
       sumher_files = [f"results/ldak/ldak-thin/ukbb/rg/{trait_pair}.cors.full" for trait_pair in ukbb_trait_pairs],
       metadata_file = "resources/ukbb_sum_stats/trait_metadata.tsv"
    output:
        "results/ldak/ldak-thin/ukbb/rg/compiled_ukbb_sumher_results.tsv"
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
