# Try groups next
rule test_random:
    input:
        block_files = lambda wildcards: get_randomised_block_files_for_pair(wildcards)[0]
    output:
        "results/simgwas/simulated_sum_stats/whole_genome_sum_stats/randomised/{ncases_A,\d+}_{ncontrols_A,\d+}_{ncases_B,\d+}_{ncontrols_B,\d+}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/{effect_blocks_A}_seed_{seed}_sum_stats_A_tag_{tag_A}_of_{tag_A}{tag_B}"
    resources:
        runtime = 1
    shell:
        "touch {output}"
    #"results/simgwas/simulated_sum_stats/whole_genome_sum_stats/randomised/{ncases_A,\d+}_{ncontrols_A,\d+}_{ncases_B,\d+}_{ncontrols_B,\d+}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/{effect_blocks_A}_seed_{seed}_sum_stats_A_tag_{tag_A}_of_{tag_A}{tag_B}.tsv.gz"

rule test_one:
    input:
        "results/simgwas/simulated_sum_stats/chr1/block_sum_stats/null/500_10000/block_0_sum_stats.tsv.gz"
    output:
        "test_one_{number}.txt"
    group: "test"
    resources:
        runtime = 1
    shell:
        "touch {output}"

rule test_two:
    input:
        "test_one_{number}.txt"
    output:
        "test_two_{number}.txt"
    group: "test"
    resources:
        runtime = 1
    shell:
        "touch {output}"
