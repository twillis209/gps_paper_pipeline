# TODO experiment with LDAK-Thin and BLD-LDAK and see if it makes a difference

rule thin_predictors:
    input:
        "resources/1000g/euro/qc/{chr}_qc.bed",
        "resources/1000g/euro/qc/{chr}_qc.bim",
        "resources/1000g/euro/qc/{chr}_qc.fam"
    output:
        thin_file = "results/simgwas/ldak/ldak-thin/weights/{chr}/thin.in",
        weights_file = "results/simgwas/ldak/ldak-thin/weights/{chr}/weights.thin"
    params:
        input_stem = "resources/1000g/euro/qc/{chr}_qc",
        output_stem = "results/simgwas/ldak/ldak-thin/weights/{chr}/thin"
    shell:
        """
        ldak --thin {params.output_stem} --bfile {params.input_stem} --window-prune .98 --window-kb 100;
        awk < {output.thin_file} '{{print $1, 1}}' > {output.weights_file}
        """

# TODO need to decide whether assuming LDSC/GCTA model for which we ignore weights and set --power 1
# TODO the LDAK-Thin model uses --power -.25
rule calculate_ldak_thin_taggings_for_chromosome:
    input:
        "resources/1000g/euro/qc/{chr}_qc.bed",
        "resources/1000g/euro/qc/{chr}_qc.bim",
        "resources/1000g/euro/qc/{chr}_qc.fam",
        weights_file = "results/simgwas/ldak/ldak-thin/weights/{chr}/weights.thin"
    output:
        tagging_file = "results/simgwas/ldak/ldak-thin/taggings/{chr}.tagging"
    params:
        input_stem = "resources/1000g/euro/qc/{chr}_qc",
        output_stem = "results/simgwas/ldak/ldak-thin/taggings/{chr}"
    shell:
        "ldak --calc-tagging {params.output_stem} --bfile {params.input_stem} --weights {input.weights_file} --chr {wildcards.chr} --window-kb 1000 --power -.25"

rule join_ldak_thin_taggings:
    input:
        [f"results/simgwas/ldak/ldak-thin/taggings/{chr}.tagging" for chr in range(1, 23)]
    output:
        chrom_taggings_file = temp("results/simgwas/ldak/ldak-thin/taggings.txt"),
        wg_tagging_file = "results/simgwas/ldak/ldak-thin/whole_genome.tagging"
    shell:
        """
        for x in {input}; do
            echo $x >> {output.chrom_taggings_file}
        done;

        ldak --join-tagging {output.wg_tagging_file} --taglist {output.chrom_taggings_file}
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
    script:
        "process_combined_sum_stats_for_sumher.R"

# TODO removing predictors explaining more than 1% of phenotypic variance?
# TODO use --genomic control or --intercept arguments?
# https://dougspeed.com/summary-statistics/
rule estimate_rg_with_ldak_thin:
    input:
        wg_tagging_file = "results/simgwas/ldak/ldak-thin/whole_genome.tagging",
        sum_stats_file_A = "results/simgwas/simulated_sum_stats/whole_genome_sum_stats/randomised/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/{effect_blocks_A}_seed_{seed}_sum_stats_A_tag_{tag_A}_of_{tag_A}{tag_B}.tsv.gz",
        sum_stats_file_B = "results/simgwas/simulated_sum_stats/whole_genome_sum_stats/randomised/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/{effect_blocks_B}_seed_{seed}_sum_stats_B_tag_{tag_B}_of_{tag_A}{tag_B}.tsv.gz"
    output:
        "results/simgwas/ldak/ldak-thin/rg/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}_seed_{seed}_{tag_A}{tag_B}.cors"
    params:
        output_stem = "results/simgwas/ldak/ldak-thin/rg/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}_seed_{seed}_{tag_A}{tag_B}"
    shell:
        """
        ldak --sum-cors {params.output_stem} --tagfile {input.wg_tagging_file} --summary {input.sum_stats_file_A} --summary2 {input.sum_stats_file_B} --allow-ambiguous YES
        """

#        rule estimate_sumher_rg:
#            input:
#            output:
#            shell:
#                "ldak --sum-cors gencor --summary {input.sum_stats_A} --summary2 {input.sum_stats_B} --tagfile LDAK-Thin.tagging --allow-ambiguous YES"
#
#rule cut_weights:
#    input:
#        "resources/1000g/euro/qc/{chr}_qc.bed",
#        "resources/1000g/euro/qc/{chr}_qc.bim",
#        "resources/1000g/euro/qc/{chr}_qc.fam"
#    output:
#        predictors = "results/simgwas/ldak/weights/{chr}/weights.predictors",
#        log = "resources/simgwas/ldak/weights/{chr}/cut_weights.log"
#    params:
#        input_stem = "resources/1000g/euro/qc/{chr}_qc",
#        output_dir = "results/simgwas/ldak/weights/{chr}"
#    shell:
#        "ldak --cut-weights {params.output_dir} --bfile {params.input_stem} > {output.log}"
#
#rule calc_weights_all:
#     input:
#        "results/simgwas/ldak/{chr}/weights.predictors",
#        "resources/1000g/euro/qc/{chr}_qc.bed",
#        "resources/1000g/euro/qc/{chr}_qc.bim",
#        "resources/1000g/euro/qc/{chr}_qc.fam"
#     output:
#         weights_file = "results/simgwas/ldak/weights/{chr}/weights.all",
#         log = "results/simgwas/ldak/weights/{chr}/calc_weights_all.log"
#     params:
#        input_stem = "resources/1000g/euro/qc/{chr}_qc",
#        output_dir = "results/simgwas/ldak/weights/{chr}"
#     shell:
#      "ldak --calc-weights-all {params.output_dir} --bfile {params.input_stem} > {output.log}"
#
## NB: Actually the GCTA model, but recommended if one wishes to copy LDSC
#
#rule calculate_ldsc_taggings_for_chromosome:
#    input:
#        weights_file
#    output:
#        output_stem
#    shell:
#        "ldak --calc-tagging {params.output_stem} --bfile {params.input_dir} --ignore-weights YES --chr {wildcards.chr} --power 1 --window-kb 1000"
#
## TODO should I add --cutoff 0.01? Read about it
