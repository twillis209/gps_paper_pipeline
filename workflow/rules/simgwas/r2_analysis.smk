r2_values = [0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.5, 0.8]
#r2_values = [0.5, 0.8]

def get_test_files_for_r2_sweep(ncases, blocks, shared_blocks, r2_values):
    return list(chain(*[get_test_files("results/simgwas/simulation_parameters.tsv", reps = 400, filetype = 'gps', subset = f"a_blocks == \'{blocks}\' & shared_blocks == \'{shared_blocks}\' & ncases_A == {ncases}", r2 = x) for x in r2_values]))+list(chain(*[get_test_files("results/simgwas/simulation_parameters.tsv", reps = 400, filetype = 'li_gps', subset = f"a_blocks == \'{blocks}\' & shared_blocks == \'{shared_blocks}\' & ncases_A == {ncases}", r2 = x) for x in r2_values]))

rule run_r2_sweep:
    input:
        lambda wildcards: list(chain(*[get_test_files("results/simgwas/simulation_parameters.tsv", reps = 400, filetype = y, subset = f"a_blocks == \'{wildcards.blocks}\' & shared_blocks == \'{wildcards.shared_blocks}\' & ncases_A == {wildcards.ncases}", r2 = x) for x in r2_values for y in ['gps', 'li_gps']]))
    output:
        "results/gps/simgwas/400_reps/randomised/{ncases}_{ncontrols}_{ncases}_{ncontrols}/{blocks}_{blocks}_{shared_blocks}/r2_sweep/r2_sweep.tsv"
    run:
        daf = pd.DataFrame(columns = ['ncases.A', 'ncontrols.A', 'ncases.B', 'ncontrols.B', 'blocks.A',
                                      'blocks.B', 'shared_blocks', 'tag_pair', 'seed', 'gps', 'pval'])

        for x in r2_values:
            gps_daf = compile_gps_results_into_daf(get_test_files("results/simgwas/simulation_parameters.tsv", reps = 400, filetype = 'gps', subset = f"a_blocks == \'{wildcards.blocks}\' & shared_blocks == \'{wildcards.shared_blocks}\' & ncases_A == {wildcards.ncases}", r2 = x))
            gps_daf.drop(columns = ['n', 'loc', 'loc.sd', 'scale', 'scale.sd', 'shape', 'shape.sd'], inplace = True)
            gps_daf = gps_daf.assign(r2 = x, dist = 'gev')

            li_gps_daf = compile_li_gps_results_into_daf(get_test_files("results/simgwas/simulation_parameters.tsv", reps = 400, filetype = 'li_gps', subset = f"a_blocks == \'{wildcards.blocks}\' & shared_blocks == \'{wildcards.shared_blocks}\' & ncases_A == {wildcards.ncases}", r2 = x))
            li_gps_daf = li_gps_daf.assign(r2 = x, dist = 'exp')

            daf = pd.concat([daf, gps_daf, li_gps_daf])

        daf.to_csv(output[0], sep = '\t', index = False)

rule run_r2_sweeps:
    input:
        "results/gps/simgwas/400_reps/randomised/10000_10000_10000_10000/s400_s400_s200/r2_sweep/r2_sweep.tsv",
        "results/gps/simgwas/400_reps/randomised/1000_10000_1000_10000/m50_m50_m25/r2_sweep/r2_sweep.tsv",
        "results/gps/simgwas/400_reps/randomised/10000_10000_10000_10000/s400_s400_s0/r2_sweep/r2_sweep.tsv",
        "results/gps/simgwas/400_reps/randomised/1000_10000_1000_10000/m50_m50_m0/r2_sweep/r2_sweep.tsv"


