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
        [f"results/simgwas/ldak/ldak-thin/taggings/{chr}.tagging" for chr in chroms]
    output:
        chrom_taggings_file = temp("results/simgwas/ldak/ldak-thin/taggings.txt")
        wg_tagging_file = "results/simgwas/ldak/ldak-thin/whole_genome.tagging"
    shell:
        """
        for x in {input}; do
            echo $x >> {output.chrom_taggings_file}
        done;

        ldak --join-tagging {output.wg_tagging_file} --taglist {output.chrom_taggings_file}
        """

rule estimate_rg_with_ldak_thin:

# TODO removing predictors explaining more than 1% of phenotypic variance?
# TODO use --genomic control or --intercept arguments?
# https://dougspeed.com/summary-statistics/
rule estimate_rg_with_ldak_thin:
    input:
        wg_tagging_file = "results/simgwas/ldak/ldak-thin/whole_genome.tagging",
        sum_stats_file_A = ,
        sum_stats_file_B = ,
    output:
        "results/simgwas/ldak/ldak-thin/rg/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}_seed_{seed}_{tag_A}{tag_B}.cors"
    params:
        output_stem = "results/simgwas/ldak/ldak-thin/rg/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}_seed_{seed}_{tag_A}{tag_B}"
    shell:
        """
        ldak --sum-cors {params.output_stem} --tagfile {input.wg_tagging_file} --summary {input.sum_stats_file_A} --summary2 {input.sum_stats_file_B} --allow-ambiguous YES
        """

rule cut_weights:
    input:
        "resources/1000g/euro/qc/{chr}_qc.bed",
        "resources/1000g/euro/qc/{chr}_qc.bim",
        "resources/1000g/euro/qc/{chr}_qc.fam"
    output:
        predictors = "results/simgwas/ldak/weights/{chr}/weights.predictors",
        log = "resources/simgwas/ldak/weights/{chr}/cut_weights.log"
    params:
        input_stem = "resources/1000g/euro/qc/{chr}_qc",
        output_dir = "results/simgwas/ldak/weights/{chr}"
    shell:
        "ldak --cut-weights {params.output_dir} --bfile {params.input_stem} > {output.log}"

rule calc_weights_all:
     input:
        "results/simgwas/ldak/{chr}/weights.predictors",
        "resources/1000g/euro/qc/{chr}_qc.bed",
        "resources/1000g/euro/qc/{chr}_qc.bim",
        "resources/1000g/euro/qc/{chr}_qc.fam"
     output:
         weights_file = "results/simgwas/ldak/weights/{chr}/weights.all",
         log = "results/simgwas/ldak/weights/{chr}/calc_weights_all.log"
     params:
        input_stem = "resources/1000g/euro/qc/{chr}_qc",
        output_dir = "results/simgwas/ldak/weights/{chr}"
     shell:
      "ldak --calc-weights-all {params.output_dir} --bfile {params.input_stem} > {output.log}"

# NB: Actually the GCTA model, but recommended if one wishes to copy LDSC
rule calculate_ldsc_taggings_for_chromosome:
    input:
        weights_file
    output:
        output_stem
    shell:
        "ldak --calc-tagging {params.output_stem} --bfile {params.input_dir} --ignore-weights YES --chr {wildcards.chr} --power 1 --window-kb 1000"

rule compile_taggings:
    input:
    output:
    params:
        output_stem
    shell:
        "ldak --join-tagging {params.output_stem} --taglist {params.tag_files}"

# TODO should I add --cutoff 0.01? Read about it
rule estimate_sumher_rg:
    input:
    output:
    shell:
        "ldak --sum-cors gencor --summary {input.sum_stats_A} --summary2 {input.sum_stats_B} --tagfile LDAK-Thin.tagging --allow-ambiguous YES"
