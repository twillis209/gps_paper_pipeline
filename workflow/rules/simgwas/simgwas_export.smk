import re
import os
from numpy import nan

include: 'simgwas_export_functions.py'

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

rule simulation_result_quartet:
    input:
        "results/ldsc/simgwas/{no_reps}_reps/randomised/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/rg/fixed_h2_free_rg_intercept/seed_{seed}_tags_{tag_A}-{tag_B}.log",
        "results/gps/simgwas/{no_reps}_reps/randomised/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/window_1000kb_step_50/3000_permutations/seed_{seed}_tags_{tag_A}-{tag_B}_gps_pvalue.tsv",
        "results/hoeffdings/simgwas/{no_reps}_reps/randomised/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/window_1000kb_step_50/seed_{seed}_tags_{tag_A}-{tag_B}_hoeffdings.tsv",
        "results/ldak/ldak-thin/simgwas/{no_reps}_reps/randomised/rg/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/seed_{seed}_tags_{tag_A}-{tag_B}.cors"
    output:
        "results/simgwas/done/{no_reps}_reps/randomised/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/seed_{seed}_tags_{tag_A}-{tag_B}.done"
    shell: "touch {output}"

rule run_m25_simulations_and_analyses:
    input:
        input_files = lambda wildcards: get_all_simulation_done_files("results/simgwas/simulation_parameters.tsv", reps = 400, subset = 'm25')
    output:
        "results/simgwas/simulated_sum_stats/whole_genome_sum_stats/400_reps/randomised/m25.done"
    shell: "touch {output}"

rule run_m50_simulations_and_analyses:
    input:
        input_files = lambda wildcards: get_all_simulation_done_files("results/simgwas/simulation_parameters.tsv", reps = 400, subset = 'm50')
    output:
        "results/simgwas/simulated_sum_stats/whole_genome_sum_stats/400_reps/randomised/m50.done"
    shell: "touch {output}"

rule run_s400_simulations_and_analyses:
    input:
        input_files = lambda wildcards: get_all_simulation_done_files("results/simgwas/simulation_parameters.tsv", reps = 400, subset = 's400')
    output:
        "results/simgwas/simulated_sum_stats/whole_genome_sum_stats/400_reps/randomised/s400.done"
    shell: "touch {output}"

rule run_s200_m25_simulations_and_analyses:
    input:
        input_files = lambda wildcards: get_all_simulation_done_files("results/simgwas/simulation_parameters.tsv", reps = 400, subset = 's200-m25')
    output:
        "results/simgwas/simulated_sum_stats/whole_genome_sum_stats/400_reps/randomised/s200-m25.done"
    shell: "touch {output}"

rule compile_test_files:
    input:
        ldsc_files = lambda wildcards: get_test_files("results/simgwas/simulation_parameters.tsv", reps = wildcards.no_reps, filetype = 'ldsc', subset = f"a_blocks == \'{wildcards.effect_blocks}\'"),
        sumher_files = lambda wildcards: get_test_files("results/simgwas/simulation_parameters.tsv", reps = wildcards.no_reps, filetype = 'sumher', subset = f"a_blocks == \'{wildcards.effect_blocks}\'"),
        hoeffdings_files = lambda wildcards: get_test_files("results/simgwas/simulation_parameters.tsv", reps = wildcards.no_reps, filetype = 'hoeffdings', subset = f"a_blocks == \'{wildcards.effect_blocks}\'"),
        gps_files = lambda wildcards: get_test_files("results/simgwas/simulation_parameters.tsv", reps = wildcards.no_reps, filetype = 'gps', subset = f"a_blocks == \'{wildcards.effect_blocks}\'"),
        theo_files = lambda wildcards: get_test_files("results/simgwas/simulation_parameters.tsv", reps = wildcards.no_reps, subset = f"a_blocks == \'{wildcards.effect_blocks}\'")
    output:
        ldsc_out = "results/ldsc/simgwas/{no_reps}_reps/randomised/compiled_{effect_blocks}_results.tsv",
        sumher_out = "results/ldak/ldak-thin/simgwas/{no_reps}_reps/randomised/rg/compiled_{effect_blocks}_results.tsv",
        hoeffdings_out = "results/hoeffdings/simgwas/{no_reps}_reps/randomised/compiled_{effect_blocks}_results.tsv",
        gps_out = "results/gps/simgwas/{no_reps}_reps/randomised/compiled_{effect_blocks}_results.tsv",
        theo_out = "results/ldsc/simgwas/{no_reps}_reps/randomised/compiled_{effect_blocks}_theo_rg.tsv",
        done_out = "results/{effect_blocks}_{no_reps}_reps.done"
    run:
        ldsc_daf = compile_ldsc_results_into_daf(input.ldsc_files)
        ldsc_daf.to_csv(output.ldsc_out, sep = '\t', index = False)

        sumher_daf = compile_sumher_results_into_daf(input.sumher_files)
        sumher_daf.to_csv(output.sumher_out, sep = '\t', index = False)

        hoeffdings_daf = compile_hoeffdings_results_into_daf(input.hoeffdings_files)
        hoeffdings_daf.to_csv(output.hoeffdings_out, sep = '\t', index = False)

        gps_daf = compile_gps_results_into_daf(input.gps_files)
        gps_daf.to_csv(output.gps_out, sep = '\t', index = False)

        theo_daf = compile_theo_rg_results_into_daf(input.theo_files)
        theo_daf.to_csv(output[0], sep = '\t', index = False)

        shell("touch {output.done_out}")
