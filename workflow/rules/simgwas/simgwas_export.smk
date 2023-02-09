import re
import os
from numpy import nan

include: 'simgwas_export_functions.py'

sample_sizes = [
                (500, 10000),
                (1000, 10000),
                (5000, 10000),
                (10000, 10000),
                (100000, 100000),
                ]

ncases = [500, 1000, 5000, 10000, 100000]
ncontrols = [10000, 10000, 10000, 10000, 100000]
sub_ncases = [1000, 5000, 10000]
sub_ncontrols = [10000, 10000, 10000]

s400_shared_blocks = ['s0', 's100', 's200', 's300', 's400']
m25_shared_blocks = ['m0', 'm5', 'm10', 'm15', 'm20', 'm25']
m50_shared_blocks = ['m0', 'm10', 'm20', 'm30', 'm40', 'm50']
s200_m25_shared_blocks = ['s0-m0', 's100-m0', 's100-m15', 's100-m25', 's200-m0', 's200-m15', 's200-m25']

localrules: simulation_result_quintet

rule run_500_null_block_simulations:
    input:
        block_files = get_all_block_files(sample_sizes = [(500, 10000)], effect = 'null')

rule run_t1000_block_simulations:
    input:
        block_files = get_all_block_files(sample_sizes = sample_sizes, effect = 'tiny')

rule simulation_result_quintet:
    input:
        "results/ldsc/simgwas/{no_reps}_reps/randomised/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/rg/fixed_h2_free_rg_intercept/seed_{seed}_tags_{tag_A}-{tag_B}.log",
        "results/gps/simgwas/{no_reps}_reps/randomised/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/window_1000kb_step_50_r2_0_2/3000_permutations/seed_{seed}_tags_{tag_A}-{tag_B}_gps_pvalue.tsv",
        "results/gps/simgwas/{no_reps}_reps/randomised/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/window_1000kb_step_50_r2_0_2/seed_{seed}_tags_{tag_A}-{tag_B}_li_gps_pvalue.tsv",
        "results/hoeffdings/simgwas/{no_reps}_reps/randomised/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/window_1000kb_step_50_r2_0_2/seed_{seed}_tags_{tag_A}-{tag_B}_hoeffdings.tsv",
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
        li_gps_files = lambda wildcards: get_test_files("results/simgwas/simulation_parameters.tsv", reps = wildcards.no_reps, filetype = 'li_gps', subset = f"a_blocks == \'{wildcards.effect_blocks}\'"),
        theo_files = lambda wildcards: get_test_files("results/simgwas/simulation_parameters.tsv", reps = wildcards.no_reps, filetype = 'theo', subset = f"a_blocks == \'{wildcards.effect_blocks}\'")
    output:
        ldsc_out = "results/ldsc/simgwas/{no_reps}_reps/randomised/compiled_{effect_blocks}_results.tsv",
        sumher_out = "results/ldak/ldak-thin/simgwas/{no_reps}_reps/randomised/rg/compiled_{effect_blocks}_results.tsv",
        hoeffdings_out = "results/hoeffdings/simgwas/{no_reps}_reps/randomised/compiled_{effect_blocks}_results.tsv",
        gps_out = "results/gps/simgwas/{no_reps}_reps/randomised/compiled_{effect_blocks}_results.tsv",
        li_gps_out = "results/gps/simgwas/{no_reps}_reps/randomised/compiled_{effect_blocks}_li_gps_results.tsv",
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

        li_gps_daf = compile_li_gps_results_into_daf(input.li_gps_files)
        li_gps_daf.to_csv(output.li_gps_out, sep = '\t', index = False)

        theo_daf = compile_theo_rg_results_into_daf(input.theo_files)
        theo_daf.to_csv(output[0], sep = '\t', index = False)

        shell("touch {output.done_out}")

rule run_s400_li_gps_simulations:
    input:
        li_gps_files = lambda wildcards: get_test_files("results/simgwas/simulation_parameters.tsv", reps = 400, filetype = 'li_gps', subset = f"a_blocks == \'s400\' & shared_blocks in {s400_shared_blocks} & ncases_A in {ncases} & ncontrols_A in {ncontrols}")

rule run_m25_li_gps_simulations:
    input:
        li_gps_files = lambda wildcards: get_test_files("results/simgwas/simulation_parameters.tsv", reps = 400, filetype = 'li_gps', subset = f"a_blocks == \'m25\' & shared_blocks in {m25_shared_blocks} & ncases_A in {ncases} & ncontrols_A in {ncontrols}")

rule run_m50_li_gps_simulations:
    input:
        li_gps_files = lambda wildcards: get_test_files("results/simgwas/simulation_parameters.tsv", reps = 400, filetype = 'li_gps', subset = f"a_blocks == \'m50\' & shared_blocks in {m50_shared_blocks} & ncases_A in {ncases} & ncontrols_A in {ncontrols}")

rule run_s200_m25_li_gps_simulations:
    input:
        li_gps_files = lambda wildcards: get_test_files("results/simgwas/simulation_parameters.tsv", reps = 400, filetype = 'li_gps', subset = f"a_blocks == \'s200-m25\' & shared_blocks in {s200_m25_shared_blocks} & ncases_A in {ncases} & ncontrols_A in {ncontrols}")

rule run_t1000_simulations_and_analyses:
    input:
        input_files = lambda wildcards: get_all_simulation_done_files("results/simgwas/simulation_parameters.tsv", reps = 400, subset = 't1000')
    output:
        "results/simgwas/simulated_sum_stats/whole_genome_sum_stats/400_reps/randomised/t1000.done"
    shell: "touch {output}"
