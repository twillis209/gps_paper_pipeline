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
    threads: 4
    priority: 1
    group: "ldsc_hoeffding_sumher_gps_sans_permutation"
    shell:
        """
        $ldakRoot/ldak --thin {params.output_stem} --bfile {params.input_stem} --window-prune .98 --window-kb 100 --max-threads {threads} > {log.log_file};
        awk < {output.thin_file} '{{print $1, 1}}' > {output.weights_file}
        """

rule calculate_ldak_thin_taggings_for_chromosome:
    input:
        "resources/1000g/euro/qc/nodup/snps_only/{join}/{chr}.bed",
        "resources/1000g/euro/qc/nodup/snps_only/{join}/{chr}.bim",
        "resources/1000g/euro/qc/nodup/snps_only/{join}/{chr}.fam",
        weights_file = "results/ldak/ldak-thin/weights/{join}/{chr}/weights.thin"
    output:
        tagging_file = "results/ldak/ldak-thin/taggings/{join}/{chr}.tagging",
    log:
        log_file = "results/ldak/ldak-thin/taggings/{join}/{chr}.tagging.log"
    params:
        input_stem = "resources/1000g/euro/qc/nodup/snps_only/{join}/{chr}",
        output_stem = "results/ldak/ldak-thin/taggings/{join}/{chr}"
    threads: 4
    priority: 1
    group: "ldsc_hoeffding_sumher_gps_sans_permutation"
    shell:
        "$ldakRoot/ldak --calc-tagging {params.output_stem} --bfile {params.input_stem} --weights {input.weights_file} --chr {wildcards.chr} --window-kb 1000 --power -.25 --max-threads {threads} > {log.log_file}"

rule join_ldak_thin_taggings:
    input:
        [f"results/ldak/ldak-thin/taggings/{{join}}/chr{x}.tagging" for x in range(1, 23)]
    output:
        wg_tagging_file = "results/ldak/ldak-thin/{join}/whole_genome.tagging",
        chrom_taggings_file = "results/ldak/ldak-thin/{join}/taggings.txt"
    log:
        log_file = "results/ldak/ldak-thin/{join}/whole_genome.tagging.log"
    params:
        output_stem = "results/ldak/ldak-thin/{join}/whole_genome"
    priority: 1
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
        "results/simgwas/simulated_sum_stats/whole_genome_sum_stats/{no_reps}_reps/randomised/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/seed_{seed}_sum_stats_{pair_label}_tag_{tag}_of_{tag_A}-{tag_B}.assoc"
    params:
        z_colname = lambda wildcards: f'zsim.{wildcards.tag}',
        chr_colname = 'chr',
        bp_colname = 'position',
        a1_colname = 'a1',
        a2_colname = 'a0',
        sample_size = lambda wildcards: int(wildcards.ncases_A)+int(wildcards.ncontrols_A)
    threads: 2
    resources:
        runtime = 5
    priority: 1
    group: "ldsc_hoeffding_sumher_gps_sans_permutation"
    script:
        "../../scripts/process_combined_simgwas_sum_stats_for_sumher.R"

rule estimate_rg_with_ldak_thin_for_simgwas:
    input:
        wg_tagging_file = "results/ldak/ldak-thin/1000g/whole_genome.tagging",
        sum_stats_file_A = "results/simgwas/simulated_sum_stats/whole_genome_sum_stats/{no_reps}_reps/randomised/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/seed_{seed}_sum_stats_A_tag_{tag_A}_of_{tag_A}-{tag_B}.assoc",
        sum_stats_file_B = "results/simgwas/simulated_sum_stats/whole_genome_sum_stats/{no_reps}_reps/randomised/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/seed_{seed}_sum_stats_B_tag_{tag_B}_of_{tag_A}-{tag_B}.assoc"
    output:
#        progress_file = "results/ldak/ldak-thin/simgwas/{no_reps}_reps/randomised/rg/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/seed_{seed}_tags_{tag_A}-{tag_B}.progress",
        cors_file = "results/ldak/ldak-thin/simgwas/{no_reps}_reps/randomised/rg/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/seed_{seed}_tags_{tag_A}-{tag_B}.cors",
#        cors_full_file = "results/ldak/ldak-thin/simgwas/{no_reps}_reps/randomised/rg/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/seed_{seed}_tags_{tag_A}-{tag_B}.cors.full",
#        labels_file = "results/ldak/ldak-thin/simgwas/{no_reps}_reps/randomised/rg/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/seed_{seed}_tags_{tag_A}-{tag_B}.labels",
#        overlap_file = "results/ldak/ldak-thin/simgwas/{no_reps}_reps/randomised/rg/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/seed_{seed}_tags_{tag_A}-{tag_B}.overlap"
    log:
        log_file = "results/ldak/ldak-thin/simgwas/{no_reps}_reps/randomised/rg/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/seed_{seed}_tags_{tag_A}-{tag_B}.log"
    params:
        output_stem = "results/ldak/ldak-thin/simgwas/{no_reps}_reps/randomised/rg/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/seed_{seed}_tags_{tag_A}-{tag_B}"
    resources:
        runtime = 5
    priority: 1
    group: "ldsc_hoeffding_sumher_gps_sans_permutation"
    shell:
        """
        $ldakRoot/ldak --sum-cors {params.output_stem} --tagfile {input.wg_tagging_file} --summary {input.sum_stats_file_A} --summary2 {input.sum_stats_file_B} --allow-ambiguous YES --check-sums NO --cutoff 0.01 > {log.log_file}
        """
