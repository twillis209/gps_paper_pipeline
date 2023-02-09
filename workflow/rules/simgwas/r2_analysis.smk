rule run_r2_sweep:
    input:
        gps_files = lambda wildcards: get_test_files("results/simgwas/simulation_parameters.tsv", reps = 400, filetype = 'gps', subset = f"a_blocks == \'s400\' & shared_blocks == \'s200\' & ncases_A == 10000"),
        li_gps_files = lambda wildcards: get_test_files("results/simgwas/simulation_parameters.tsv", reps = 400, filetype = 'li_gps', subset = f"a_blocks == \'s400\' & shared_blocks == \'s200\' & ncases_A == 10000")
