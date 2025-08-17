import pandas as pd
import sys

exec(open('workflow/rules/simgwas/simgwas_experiments_functions.py', 'r').read())

effect_blocks = ['m25', 'm50', 's400', 's200-m25', 'all']

s400_shared_blocks = ['s0', 's50', 's100', 's150', 's200', 's250', 's300', 's350', 's400']

m25_shared_blocks = ['m0', 'm5', 'm10', 'm15', 'm20', 'm25']

m50_shared_blocks = ['m0', 'm10', 'm20', 'm30', 'm40', 'm50']

s200_m25_shared_blocks = ['s0-m0', 's100-m0', 's100-m15', 's100-m25', 's200-m0', 's200-m15', 's200-m25']

if not sys.argv[1]:
    raise Exception("Must specify effect block as first argument")
elif sys.argv[1] not in effect_blocks:
    raise Exception(f"Must specify effect block in {', '.join(effect_blocks)}")
else:
    effect_block = sys.argv[1]

# s400
if effect_block in ['s400', 'all']:
    s400_ldsc_files = get_existing_test_files("results/simgwas/simulation_parameters.tsv", reps = 400, filetype = 'ldsc', subset = "a_blocks == \'s400\' & shared_blocks in @s400_shared_blocks")

    s400_ldsc_daf = compile_ldsc_results_into_daf(s400_ldsc_files)

    s400_ldsc_daf.to_csv("results/ldsc/simgwas/400_reps/randomised/manuscript_compiled_s400_results.tsv", sep = '\t', index = False)

    s400_sumher_files = get_existing_test_files("results/simgwas/simulation_parameters.tsv", reps = 400, filetype = 'sumher', subset = "a_blocks == \'s400\' & shared_blocks in @s400_shared_blocks")

    s400_sumher_daf = compile_sumher_results_into_daf(s400_sumher_files)

    s400_sumher_daf.to_csv("results/ldak/ldak-thin/simgwas/400_reps/randomised/rg/manuscript_compiled_s400_results.tsv", sep = '\t', index = False)

    s400_gps_files = get_existing_test_files("results/simgwas/simulation_parameters.tsv", reps = 400, filetype = 'gps', subset = "a_blocks == \'s400\' & shared_blocks in @s400_shared_blocks")

    s400_gps_daf = compile_gps_results_into_daf(s400_gps_files)

    s400_gps_daf.to_csv("results/gps/simgwas/400_reps/randomised/manuscript_compiled_s400_results.tsv", sep = '\t', index = False)

    s400_hoeffdings_files = get_existing_test_files("results/simgwas/simulation_parameters.tsv", reps = 400, filetype = 'hoeffdings', subset = "a_blocks == \'s400\' & shared_blocks in @s400_shared_blocks")

    s400_hoeffdings_daf = compile_hoeffdings_results_into_daf(s400_hoeffdings_files)

    s400_hoeffdings_daf.to_csv("results/hoeffdings/simgwas/400_reps/randomised/manuscript_compiled_s400_results.tsv", sep = '\t', index = False)

    s400_theo_files = get_existing_test_files("results/simgwas/simulation_parameters.tsv", reps = 400, filetype = 'theo', subset = "a_blocks == \'s400\' & shared_blocks in @s400_shared_blocks")

    s400_theo_daf = compile_theo_results_into_daf(s400_theo_files)

    s400_theo_daf.to_csv("results/ldsc/simgwas/400_reps/randomised/manuscript_compiled_s400_theo_results.tsv", sep = '\t', index = False)

# m25
if effect_block in ['m25', 'all']:
    m25_ldsc_files = get_existing_test_files("results/simgwas/simulation_parameters.tsv", reps = 400, filetype = 'ldsc', subset = "a_blocks == \'m25\' & shared_blocks in @m25_shared_blocks")

    m25_ldsc_daf = compile_ldsc_results_into_daf(m25_ldsc_files)

    m25_ldsc_daf.to_csv("results/ldsc/simgwas/400_reps/randomised/manuscript_compiled_m25_results.tsv", sep = '\t', index = False)

    m25_sumher_files = get_existing_test_files("results/simgwas/simulation_parameters.tsv", reps = 400, filetype = 'sumher', subset = "a_blocks == \'m25\' & shared_blocks in @m25_shared_blocks")

    m25_sumher_daf = compile_sumher_results_into_daf(m25_sumher_files)

    m25_sumher_daf.to_csv("results/ldak/ldak-thin/simgwas/400_reps/randomised/rg/manuscript_compiled_m25_results.tsv", sep = '\t', index = False)

    m25_gps_files = get_existing_test_files("results/simgwas/simulation_parameters.tsv", reps = 400, filetype = 'gps', subset = "a_blocks == \'m25\' & shared_blocks in @m25_shared_blocks")

    m25_gps_daf = compile_gps_results_into_daf(m25_gps_files)

    m25_gps_daf.to_csv("results/gps/simgwas/400_reps/randomised/manuscript_compiled_m25_results.tsv", sep = '\t', index = False)

    m25_hoeffdings_files = get_existing_test_files("results/simgwas/simulation_parameters.tsv", reps = 400, filetype = 'hoeffdings', subset = "a_blocks == \'m25\' & shared_blocks in @m25_shared_blocks")

    m25_hoeffdings_daf = compile_hoeffdings_results_into_daf(m25_hoeffdings_files)

    m25_hoeffdings_daf.to_csv("results/hoeffdings/simgwas/400_reps/randomised/manuscript_compiled_m25_results.tsv", sep = '\t', index = False)

    m25_theo_files = get_existing_test_files("results/simgwas/simulation_parameters.tsv", reps = 400, filetype = 'theo', subset = "a_blocks == \'m25\' & shared_blocks in @m25_shared_blocks")

    m25_theo_daf = compile_theo_results_into_daf(m25_theo_files)

    m25_theo_daf.to_csv("results/ldsc/simgwas/400_reps/randomised/manuscript_compiled_m25_theo_results.tsv", sep = '\t', index = False)

# m50
if effect_block in ['m50', 'all']:
    m50_ldsc_files = get_existing_test_files("results/simgwas/simulation_parameters.tsv", reps = 400, filetype = 'ldsc', subset = "a_blocks == \'m50\' & shared_blocks in @m50_shared_blocks")

    m50_ldsc_daf = compile_ldsc_results_into_daf(m50_ldsc_files)

    m50_ldsc_daf.to_csv("results/ldsc/simgwas/400_reps/randomised/manuscript_compiled_m50_results.tsv", sep = '\t', index = False)

    m50_sumher_files = get_existing_test_files("results/simgwas/simulation_parameters.tsv", reps = 400, filetype = 'sumher', subset = "a_blocks == \'m50\' & shared_blocks in @m50_shared_blocks")

    m50_sumher_daf = compile_sumher_results_into_daf(m50_sumher_files)

    m50_sumher_daf.to_csv("results/ldak/ldak-thin/simgwas/400_reps/randomised/rg/manuscript_compiled_m50_results.tsv", sep = '\t', index = False)

    m50_gps_files = get_existing_test_files("results/simgwas/simulation_parameters.tsv", reps = 400, filetype = 'gps', subset = "a_blocks == \'m50\' & shared_blocks in @m50_shared_blocks")

    m50_gps_daf = compile_gps_results_into_daf(m50_gps_files)

    m50_gps_daf.to_csv("results/gps/simgwas/400_reps/randomised/manuscript_compiled_m50_results.tsv", sep = '\t', index = False)

    m50_hoeffdings_files = get_existing_test_files("results/simgwas/simulation_parameters.tsv", reps = 400, filetype = 'hoeffdings', subset = "a_blocks == \'m50\' & shared_blocks in @m50_shared_blocks")

    m50_hoeffdings_daf = compile_hoeffdings_results_into_daf(m50_hoeffdings_files)

    m50_hoeffdings_daf.to_csv("results/hoeffdings/simgwas/400_reps/randomised/manuscript_compiled_m50_results.tsv", sep = '\t', index = False)

    m50_theo_files = get_existing_test_files("results/simgwas/simulation_parameters.tsv", reps = 400, filetype = 'theo', subset = "a_blocks == \'m50\' & shared_blocks in @m50_shared_blocks")

    m50_theo_daf = compile_theo_results_into_daf(m50_theo_files)

    m50_theo_daf.to_csv("results/ldsc/simgwas/400_reps/randomised/manuscript_compiled_m50_theo_results.tsv", sep = '\t', index = False)

# s200-m25
if effect_block in ['s200-m25', 'all']:
    s200_m25_ldsc_files = get_existing_test_files("results/simgwas/simulation_parameters.tsv", reps = 400, filetype = 'ldsc', subset = "a_blocks == \'s200-m25\' & shared_blocks in @s200_m25_shared_blocks")

    s200_m25_ldsc_daf = compile_ldsc_results_into_daf(s200_m25_ldsc_files)

    s200_m25_ldsc_daf.to_csv("results/ldsc/simgwas/400_reps/randomised/manuscript_compiled_s200_m25_results.tsv", sep = '\t', index = False)

    s200_m25_sumher_files = get_existing_test_files("results/simgwas/simulation_parameters.tsv", reps = 400, filetype = 'sumher', subset = "a_blocks == \'s200-m25\' & shared_blocks in @s200_m25_shared_blocks")

    s200_m25_sumher_daf = compile_sumher_results_into_daf(s200_m25_sumher_files)

    s200_m25_sumher_daf.to_csv("results/ldak/ldak-thin/simgwas/400_reps/randomised/rg/manuscript_compiled_s200_m25_results.tsv", sep = '\t', index = False)

    s200_m25_gps_files = get_existing_test_files("results/simgwas/simulation_parameters.tsv", reps = 400, filetype = 'gps', subset = "a_blocks == \'s200-m25\' & shared_blocks in @s200_m25_shared_blocks")

    s200_m25_gps_daf = compile_gps_results_into_daf(s200_m25_gps_files)

    s200_m25_gps_daf.to_csv("results/gps/simgwas/400_reps/randomised/manuscript_compiled_s200_m25_results.tsv", sep = '\t', index = False)

    s200_m25_hoeffdings_files = get_existing_test_files("results/simgwas/simulation_parameters.tsv", reps = 400, filetype = 'hoeffdings', subset = "a_blocks == \'s200-m25\' & shared_blocks in @s200_m25_shared_blocks")

    s200_m25_hoeffdings_daf = compile_hoeffdings_results_into_daf(s200_m25_hoeffdings_files)

    s200_m25_hoeffdings_daf.to_csv("results/hoeffdings/simgwas/400_reps/randomised/manuscript_compiled_s200_m25_results.tsv", sep = '\t', index = False)

    s200_m25_theo_files = get_existing_test_files("results/simgwas/simulation_parameters.tsv", reps = 400, filetype = 'theo', subset = "a_blocks == \'s200-m25\' & shared_blocks in @s200_m25_shared_blocks")

    s200_m25_theo_daf = compile_theo_results_into_daf(s200_m25_theo_files)

    s200_m25_theo_daf.to_csv("results/ldsc/simgwas/400_reps/randomised/manuscript_compiled_s200_m25_theo_results.tsv", sep = '\t', index = False)
