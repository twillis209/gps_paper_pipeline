import pandas as pd

def get_one_chrom_simulation_done_files(simulation_pars_file, reps, subset = None):
    daf = pd.read_csv(simulation_pars_file, sep = '\t')

    if subset:
        daf = daf.query('a_blocks == @subset & b_blocks == @subset')

    daf['r2'] = daf['r2'].apply(lambda x: str(x).replace('.', '_'))

    return [f"results/simgwas/done/{reps}_reps/randomised/{row.chr}/{row.ncases_A}_{row.ncontrols_A}_{row.ncases_B}_{row.ncontrols_B}/{row.a_blocks}_{row.b_blocks}_{row.shared_blocks}/window_{row.window}_step_{row.step}_r2_{row.r2}/seed_{row.seed}_tags_{row.tag_A}-{row.tag_B}.done" for row in daf.itertuples()]
