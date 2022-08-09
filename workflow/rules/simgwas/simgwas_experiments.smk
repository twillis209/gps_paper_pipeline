import re
import os
from numpy import nan

include: 'simgwas_experiments_functions.py'

odds_ratio_dict = {"s": 1.05, "m": 1.2, 'l': 1.4, 'v': 2, 'r': 'random', 'n' : 1, 'i' : 1.1}

sample_sizes = [(500, 10000),
                (1000, 10000),
                (2000, 10000),
                (5000, 10000),
                (10000, 10000),
                (50000, 50000),
                (100000, 100000),
                (250000, 250000)]

localrules: simulation_result_quartet

rule run_all_block_simulations:
    input:
        block_files = get_all_block_files(sample_sizes)

rule write_out_simulation_parameters_file:
    output:
        "results/simgwas/simulation_parameters.tsv"
    params:
        sample_sizes = sample_sizes
    script: "../../scripts/simgwas/write_out_simulation_parameters.py"

rule simulation_result_quartet_with_values:
    input:
        "results/ldsc/simgwas/400_reps/randomised/500_10000_500_10000/s400_s400_s0/rg/fixed_h2_free_rg_intercept/seed_1_tags_1-2.log",
        "results/gps/simgwas/400_reps/randomised/500_10000_500_10000/s400_s400_s0/window_1000kb_step_50/3000_permutations/seed_1_tags_1-2_gps_pvalue.tsv",
        "results/hoeffdings/simgwas/400_reps/randomised/500_10000_500_10000/s400_s400_s0/window_1000kb_step_50/seed_1_tags_1-2_hoeffdings.tsv",
        "results/ldak/ldak-thin/simgwas/400_reps/randomised/rg/500_10000_500_10000/s400_s400_s0/seed_1_tags_1-2.cors"

rule simulation_result_quartet:
    input:
        "results/ldsc/simgwas/{no_reps}_reps/randomised/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/rg/fixed_h2_free_rg_intercept/seed_{seed}_tags_{tag_A}-{tag_B}.log",
        "results/gps/simgwas/{no_reps}_reps/randomised/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/window_1000kb_step_50/3000_permutations/seed_{seed}_tags_{tag_A}-{tag_B}_gps_pvalue.tsv",
        "results/hoeffdings/simgwas/{no_reps}_reps/randomised/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/window_1000kb_step_50/seed_{seed}_tags_{tag_A}-{tag_B}_hoeffdings.tsv",
        "results/ldak/ldak-thin/simgwas/{no_reps}_reps/randomised/rg/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/seed_{seed}_tags_{tag_A}-{tag_B}.cors"
    output:
        "results/simgwas/done/{no_reps}_reps/randomised/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/seed_{seed}_tags_{tag_A}-{tag_B}.done"
    shell: "touch {output}"

# TODO simulation parameters probably needs to be a resource as we can't induce its production with this rule
rule run_all_simulations:
    input:
        input_files = get_all_simulation_done_files("results/simgwas/simulation_parameters.tsv", reps = 400)