import pandas as pd
import os

daf = pd.read_csv("results/simgwas/simulation_parameters.tsv", sep = '\t')

subset = 's400'

reps = 400

daf = daf.query('a_blocks == @subset & b_blocks == @subset')

file_bools = [os.path.exists(f"results/simgwas/done/{reps}_reps/randomised/{row.ncases_A}_{row.ncontrols_A}_{row.ncases_B}_{row.ncontrols_B}/{row.a_blocks}_{row.b_blocks}_{row.shared_blocks}/seed_{row.seed}_tags_{row.tag_A}-{row.tag_B}.done") for row in daf.itertuples()]

daf['exists'] = file_bools

daf.to_csv('state_of_done_files.tsv', sep = '\t', index = False)
