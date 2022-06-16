import os

rule compute_gps_for_sim_pair:
    input:
        sum_stats_file = "results/simgwas/simulated_sum_stats/whole_genome_sum_stats/randomised/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/window_{window}_step_{step}/seed_{seed}_pruned_sum_stats_tags_{tag_A}{tag_B}.tsv"
    output:
        temp("results/gps/simgwas/randomised/window_{window}_step_{step}/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}_seed_{seed}_tags_{tag_A}{tag_B}_gps_value.tsv")
    params:
        a_colname = lambda wildcards: "p.%d" % (tags.index(wildcards.tag_A)+1),
        b_colname = lambda wildcards: "p.%d" % (tags.index(wildcards.tag_B)+1)
    resources:
        time = 2
    group: "ldsc_hoeffding_and_gps_sans_permutation"
    shell:
        "workflow/scripts/gps_cpp/build/apps/computeGpsCLI -i {input.sum_stats_file} -a {params.a_colname} -b {params.b_colname} -c {wildcards.effect_blocks_A} -d {wildcards.effect_blocks_B} -l -o {output}"

rule permute_sim_pair:
    input:
        sum_stats_file = "results/simgwas/simulated_sum_stats/whole_genome_sum_stats/randomised/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/window_{window}_step_{step}/seed_{seed}_pruned_sum_stats_tags_{tag_A}{tag_B}.tsv"
    output:
        "results/gps/simgwas/randomised/window_{window}_step_{step}/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{draws,\d+}_permutations/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}_seed_{seed}_tags_{tag_A,[a-z]}{tag_B,[a-z]}.tsv"
    params:
        a_colname = lambda wildcards: "p.%d" % (tags.index(wildcards.tag_A)+1),
        b_colname = lambda wildcards: "p.%d" % (tags.index(wildcards.tag_B)+1)
    threads: 8
    resources:
        mem_mb = get_mem_mb,
        time = get_permute_time,
    group: "permutation"
    shell:
        "workflow/scripts/gps_cpp/build/apps/permuteTraitsCLI -i {input.sum_stats_file} -o {output} -a {params.a_colname} -b {params.b_colname} -c {threads} -n {wildcards.draws}"

rule fit_gev_and_compute_gps_pvalue_for_sim_pair:
    input:
        gps_file = "results/gps/simgwas/randomised/window_{window}_step_{step}/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}_seed_{seed}_tags_{tag_A}{tag_B}_gps_value.tsv",
        perm_file = "results/gps/simgwas/randomised/window_{window}_step_{step}/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{draws}_permutations/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}_seed_{seed}_tags_{tag_A}{tag_B}.tsv"
    output:
        "results/gps/simgwas/randomised/window_{window}_step_{step}/{ncases_A,\d+}_{ncontrols_A,\d+}_{ncases_B,\d+}_{ncontrols_B,\d+}/{draws,\d+}_permutations/{effect_blocks_A,[smlh\d-]+}_{effect_blocks_B,[smlh\d-]+}_{shared_effect_blocks,[smlh\d-]+}_seed_{seed,\d+}_tags_{tag_A,[a-z]}{tag_B,[a-z]}_gps_pvalue.tsv"
    resources:
        time = 2
    group: "permutation"
    shell:
      "Rscript workflow/scripts/fit_gev_and_compute_gps_pvalue.R -g {input.gps_file} -p {input.perm_file} -a {wildcards.effect_blocks_A} -b {wildcards.effect_blocks_B} -o {output}"

rule compute_hoeffdings_for_sim_pair:
    input:
        sum_stats_file = "results/simgwas/simulated_sum_stats/whole_genome_sum_stats/randomised/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/window_{window}_step_{step}/seed_{seed}_pruned_sum_stats_tags_{tag_A}{tag_B}.tsv"
    output:
        "results/hoeffdings/simgwas/randomised/window_{window}_step_{step}/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}_seed_{seed,\d+}_tags_{tag_A,[a-z]}{tag_B,[a-z]}_hoeffdings.tsv"
    params:
        a_colname = lambda wildcards: "p.%d" % (tags.index(wildcards.tag_A)+1),
        b_colname = lambda wildcards: "p.%d" % (tags.index(wildcards.tag_B)+1)
    group: "ldsc_hoeffding_and_gps_sans_permutation"
    resources:
        time = 2
    shell:
        """
        Rscript workflow/scripts/compute_hoeffdings.R -i {input.sum_stats_file} -a {params.a_colname} -b {params.b_colname} -o {output} -nt 1
        sed -i 's/{params.a_colname}/{wildcards.effect_blocks_A}_{wildcards.shared_effect_blocks}/' {output}
        sed -i 's/{params.b_colname}/{wildcards.effect_blocks_B}_{wildcards.shared_effect_blocks}/' {output}
        """

rule generate_ecdf_values_for_simgwas_trait_pair:
    input:
        sum_stats_file = "results/simgwas/simulated_sum_stats/whole_genome_sum_stats/randomised/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/window_{window}_step_{step}/seed_{seed}_pruned_sum_stats_tags_{tag_A,[a-z]}{tag_B,[a-z]}.tsv"
    output:
        temp("results/gps/simgwas/window_{window}_step_{step}/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}_seed_{seed}_{tag_A}{tag_B}_ecdf.tsv")
    params:
        pvalue_col_one = lambda wildcards: tag_pvalue_dict[wildcards.tag_A],
        pvalue_col_two = lambda wildcards: tag_pvalue_dict[wildcards.tag_B]
    shell:
        # TODO fix this: need to number p.x and p.y with x and y being the 
        "workflow/scripts/gps_cpp/build/apps/fitAndEvaluateEcdfsCLI -i {input.sum_stats_file} -a {params.pvalue_col_one} -b {params.pvalue_col_two} -o {output}"
