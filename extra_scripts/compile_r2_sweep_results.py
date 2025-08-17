import sys

exec(open('workflow/rules/simgwas/simgwas_export_functions.py', 'r').read())

r2_values = [0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.5, 0.8]

daf = pd.DataFrame(columns = ['ncases.A', 'ncontrols.A', 'ncases.B', 'ncontrols.B', 'blocks.A',
       'blocks.B', 'shared_blocks', 'tag_pair', 'seed', 'gps', 'pval'])

for x in r2_values:
    gps_daf = compile_gps_results_into_daf(get_test_files("results/simgwas/simulation_parameters.tsv", reps = 400, filetype = 'gps', subset = f"a_blocks == \'m50\' & shared_blocks == \'m25\' & ncases_A == 1000", r2 = x))
    gps_daf.drop(columns = ['n', 'loc', 'loc.sd', 'scale', 'scale.sd', 'shape', 'shape.sd'], inplace = True)
    gps_daf = gps_daf.assign(r2 = x, dist = 'gev')

    # Not getting test files here correctly
    li_gps_daf = compile_li_gps_results_into_daf(get_test_files("results/simgwas/simulation_parameters.tsv", reps = 400, filetype = 'li_gps', subset = f"a_blocks == \'m50\' & shared_blocks == \'m25\' & ncases_A == 1000", r2 = x))
    li_gps_daf = li_gps_daf.assign(r2 = x, dist = 'exp')

    daf = pd.concat([daf, gps_daf, li_gps_daf])

daf.to_csv('results/gps/simgwas/400_reps/randomised/1000_10000_1000_10000/m50_m50_m25/r2_sweep/r2_sweep.tsv', sep = '\t', index = False)

daf = pd.DataFrame(columns = ['ncases.A', 'ncontrols.A', 'ncases.B', 'ncontrols.B', 'blocks.A',
       'blocks.B', 'shared_blocks', 'tag_pair', 'seed', 'gps', 'pval'])

for x in r2_values:
    gps_daf = compile_gps_results_into_daf(get_test_files("results/simgwas/simulation_parameters.tsv", reps = 400, filetype = 'gps', subset = f"a_blocks == \'s400\' & shared_blocks == \'s200\' & ncases_A == 10000", r2 = x))
    gps_daf.drop(columns = ['n', 'loc', 'loc.sd', 'scale', 'scale.sd', 'shape', 'shape.sd'], inplace = True)
    gps_daf = gps_daf.assign(r2 = x, dist = 'gev')

    # Not getting test files here correctly
    li_gps_daf = compile_li_gps_results_into_daf(get_test_files("results/simgwas/simulation_parameters.tsv", reps = 400, filetype = 'li_gps', subset = f"a_blocks == \'s400\' & shared_blocks == \'s200\' & ncases_A == 10000", r2 = x))
    li_gps_daf = li_gps_daf.assign(r2 = x, dist = 'exp')

    daf = pd.concat([daf, gps_daf, li_gps_daf])

daf.to_csv('results/gps/simgwas/400_reps/randomised/10000_10000_10000_10000/s400_s400_s200/r2_sweep/r2_sweep.tsv', sep = '\t', index = False)

daf = pd.DataFrame(columns = ['ncases.A', 'ncontrols.A', 'ncases.B', 'ncontrols.B', 'blocks.A',
       'blocks.B', 'shared_blocks', 'tag_pair', 'seed', 'gps', 'pval'])

for x in r2_values:
    gps_daf = compile_gps_results_into_daf(get_test_files("results/simgwas/simulation_parameters.tsv", reps = 400, filetype = 'gps', subset = f"a_blocks == \'m50\' & shared_blocks == \'m0\' & ncases_A == 1000", r2 = x))
    gps_daf.drop(columns = ['n', 'loc', 'loc.sd', 'scale', 'scale.sd', 'shape', 'shape.sd'], inplace = True)
    gps_daf = gps_daf.assign(r2 = x, dist = 'gev')

    # Not getting test files here correctly
    li_gps_daf = compile_li_gps_results_into_daf(get_test_files("results/simgwas/simulation_parameters.tsv", reps = 400, filetype = 'li_gps', subset = f"a_blocks == \'m50\' & shared_blocks == \'m0\' & ncases_A == 1000", r2 = x))
    li_gps_daf = li_gps_daf.assign(r2 = x, dist = 'exp')

    daf = pd.concat([daf, gps_daf, li_gps_daf])

daf.to_csv('results/gps/simgwas/400_reps/randomised/1000_10000_1000_10000/m50_m50_m0/r2_sweep/r2_sweep.tsv', sep = '\t', index = False)

daf = pd.DataFrame(columns = ['ncases.A', 'ncontrols.A', 'ncases.B', 'ncontrols.B', 'blocks.A',
       'blocks.B', 'shared_blocks', 'tag_pair', 'seed', 'gps', 'pval'])

for x in r2_values:
    gps_daf = compile_gps_results_into_daf(get_test_files("results/simgwas/simulation_parameters.tsv", reps = 400, filetype = 'gps', subset = f"a_blocks == \'s400\' & shared_blocks == \'s0\' & ncases_A == 10000", r2 = x))
    gps_daf.drop(columns = ['n', 'loc', 'loc.sd', 'scale', 'scale.sd', 'shape', 'shape.sd'], inplace = True)
    gps_daf = gps_daf.assign(r2 = x, dist = 'gev')

    # Not getting test files here correctly
    li_gps_daf = compile_li_gps_results_into_daf(get_test_files("results/simgwas/simulation_parameters.tsv", reps = 400, filetype = 'li_gps', subset = f"a_blocks == \'s400\' & shared_blocks == \'s0\' & ncases_A == 10000", r2 = x))
    li_gps_daf = li_gps_daf.assign(r2 = x, dist = 'exp')

    daf = pd.concat([daf, gps_daf, li_gps_daf])

daf.to_csv('results/gps/simgwas/400_reps/randomised/10000_10000_10000_10000/s400_s400_s0/r2_sweep/r2_sweep.tsv', sep = '\t', index = False)
