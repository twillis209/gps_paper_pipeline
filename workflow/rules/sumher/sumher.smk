# TODO experiment with LDAK-Thin and BLD-LDAK and see if it makes a difference

rule thin_predictors:
    input:
        "resources/1000g/euro/qc/nodup/snps_only/{chr}.bed",
        "resources/1000g/euro/qc/nodup/snps_only/{chr}.bim",
        "resources/1000g/euro/qc/nodup/snps_only/{chr}.fam"
    output:
        thin_file = "results/simgwas/ldak/ldak-thin/weights/{chr}/thin.in",
        weights_file = "results/simgwas/ldak/ldak-thin/weights/{chr}/weights.thin"
    log:
        log_file = "results/simgwas/ldak/ldak-thin/weights/{chr}/thin.log"
    params:
        input_stem = "resources/1000g/euro/qc/nodup/snps_only/{chr}",
        output_stem = "results/simgwas/ldak/ldak-thin/weights/{chr}/thin"
    shell:
        """
        $ldakRoot/ldak --thin {params.output_stem} --bfile {params.input_stem} --window-prune .98 --window-kb 100 > {log.log_file};
        awk < {output.thin_file} '{{print $1, 1}}' > {output.weights_file}
        """

# TODO need to decide whether assuming LDSC/GCTA model for which we ignore weights and set --power 1
# TODO the LDAK-Thin model uses --power -.25
rule calculate_ldak_thin_taggings_for_chromosome:
    input:
        "resources/1000g/euro/qc/nodup/snps_only/{chr}.bed",
        "resources/1000g/euro/qc/nodup/snps_only/{chr}.bim",
        "resources/1000g/euro/qc/nodup/snps_only/{chr}.fam",
        weights_file = "results/simgwas/ldak/ldak-thin/weights/{chr}/weights.thin"
    output:
        tagging_file = temp("results/simgwas/ldak/ldak-thin/taggings/{chr}.tagging"),
    log:
        log_file = "results/simgwas/ldak/ldak-thin/taggings/{chr}.tagging.log"
    params:
        input_stem = "resources/1000g/euro/qc/nodup/snps_only/{chr}",
        output_stem = "results/simgwas/ldak/ldak-thin/taggings/{chr}"
    shell:
        "$ldakRoot/ldak --calc-tagging {params.output_stem} --bfile {params.input_stem} --weights {input.weights_file} --chr {wildcards.chr} --window-kb 1000 --power -.25 > {log.log_file}"

rule join_ldak_thin_taggings:
    input:
        [f"results/simgwas/ldak/ldak-thin/taggings/chr{x}.tagging" for x in range(1, 23)]
    output:
        wg_tagging_file = "results/simgwas/ldak/ldak-thin/whole_genome.tagging",
        chrom_taggings_file = temp("results/simgwas/ldak/ldak-thin/taggings.txt")
    log:
        log_file = "results/simgwas/ldak/ldak-thin/whole_genome.tagging.log"
    params:
        output_stem = "results/simgwas/ldak/ldak-thin/whole_genome"
    shell:
        """
        for x in {input}; do
            echo $x >> {output.chrom_taggings_file}
        done;

        $ldakRoot/ldak --join-tagging {params.output_stem} --taglist {output.chrom_taggings_file} > {log.log_file}
        """

# NB: Currently assuming in the script that ncases_A == ncases_B and ncontrols_A == ncontrols_B
rule process_combined_sum_stats_for_sumher:
    input:
        "results/simgwas/simulated_sum_stats/whole_genome_sum_stats/randomised/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/{effect_blocks}_seed_{seed}_sum_stats_{pair_label}_tag_{tag}_of_{tags}.tsv.gz"
    output:
        temp("results/simgwas/simulated_sum_stats/whole_genome_sum_stats/randomised/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/{effect_blocks}_seed_{seed}_sum_stats_{pair_label}_tag_{tag}_of_{tags}.assoc")
    params:
        z_colname = lambda wildcards: f'zsim.{tags.index(wildcards.tag)+1}',
        chr_colname = 'chr',
        bp_colname = 'position',
        a1_colname = 'a0',
        a2_colname = 'a1',
        sample_size = lambda wildcards: int(wildcards.ncases_A)+int(wildcards.ncontrols_A)
    threads: 2
    group: "ldsc_hoeffding_and_gps_sans_permutation"
    script:
        "process_combined_sum_stats_for_sumher.R"

# TODO removing predictors explaining more than 1% of phenotypic variance?
# TODO use --genomic control or --intercept arguments?
# https://dougspeed.com/summary-statistics/
rule estimate_rg_with_ldak_thin:
    input:
        wg_tagging_file = "results/simgwas/ldak/ldak-thin/whole_genome.tagging",
        sum_stats_file_A = "results/simgwas/simulated_sum_stats/whole_genome_sum_stats/randomised/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/{effect_blocks_A}_seed_{seed}_sum_stats_A_tag_{tag_A}_of_{tag_A}{tag_B}.assoc",
        sum_stats_file_B = "results/simgwas/simulated_sum_stats/whole_genome_sum_stats/randomised/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/{effect_blocks_B}_seed_{seed}_sum_stats_B_tag_{tag_B}_of_{tag_A}{tag_B}.assoc"
    output:
        cors_file = "results/simgwas/ldak/ldak-thin/rg/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}_seed_{seed}_{tag_A}{tag_B}.cors",
        cors_full_file = "results/simgwas/ldak/ldak-thin/rg/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}_seed_{seed}_{tag_A}{tag_B}.cors.full",
        labels_file = "results/simgwas/ldak/ldak-thin/rg/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}_seed_{seed}_{tag_A}{tag_B}.cors.labels",
        overlap_file = "results/simgwas/ldak/ldak-thin/rg/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}_seed_{seed}_{tag_A}{tag_B}.cors.overlap",
        progress_file = "results/simgwas/ldak/ldak-thin/rg/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}_seed_{seed}_{tag_A}{tag_B}.cors.progress"
    log:
        log_file = "results/simgwas/ldak/ldak-thin/rg/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}_seed_{seed}_{tag_A}{tag_B}.log"
    params:
        output_stem = "results/simgwas/ldak/ldak-thin/rg/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}_seed_{seed}_{tag_A}{tag_B}"
    group: "ldsc_hoeffding_and_gps_sans_permutation"
    shell:
        """
        $ldakRoot/ldak --sum-cors {params.output_stem} --tagfile {input.wg_tagging_file} --summary {input.sum_stats_file_A} --summary2 {input.sum_stats_file_B} --allow-ambiguous YES --check-sums NO --cutoff 0.01 > {log.log_file}
        """