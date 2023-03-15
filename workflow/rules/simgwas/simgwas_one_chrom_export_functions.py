import pandas as pd

def get_one_chrom_simulation_done_files(simulation_pars_file, reps, subset = None):
    daf = pd.read_csv(simulation_pars_file, sep = '\t')

    if subset:
        daf = daf.query('a_blocks == @subset & b_blocks == @subset')

    daf['r2'] = daf['r2'].apply(lambda x: str(x).replace('.', '_'))

    return [f"results/simgwas/done/{reps}_reps/randomised/{row.chr}/{row.ncases_A}_{row.ncontrols_A}_{row.ncases_B}_{row.ncontrols_B}/{row.a_blocks}_{row.b_blocks}_{row.shared_blocks}/window_{row.window}_step_{row.step}_r2_{row.r2}/seed_{row.seed}_tags_{row.tag_A}-{row.tag_B}.done" for row in daf.itertuples()]

def get_one_chrom_test_files(simulation_pars_file, reps, filetype, subset = None, draws = None):
    daf = pd.read_csv(simulation_pars_file, sep = '\t')

    if subset:
        daf = daf.query('a_blocks == @subset & b_blocks == @subset')

    daf['r2'] = daf['r2'].apply(lambda x: str(x).replace('.', '_'))

    if filetype == 'ldsc':
        files = [f"results/ldsc/simgwas/{reps}_reps/randomised/{row.chr}/{row.ncases_A}_{row.ncontrols_A}_{row.ncases_B}_{row.ncontrols_B}/{row.a_blocks}_{row.b_blocks}_{row.shared_blocks}/rg/fixed_h2_free_rg_intercept/seed_{row.seed}_tags_{row.tag_A}-{row.tag_B}.log" for row in daf.itertuples()]
    elif filetype == 'sumher':
        files = [f"results/ldak/ldak-thin/simgwas/{reps}_reps/randomised/rg/{row.chr}/{row.ncases_A}_{row.ncontrols_A}_{row.ncases_B}_{row.ncontrols_B}/{row.a_blocks}_{row.b_blocks}_{row.shared_blocks}/seed_{row.seed}_tags_{row.tag_A}-{row.tag_B}.cors" for row in daf.itertuples()]
    elif filetype == 'hoeffdings':
        files = [f"results/hoeffdings/simgwas/{reps}_reps/randomised/{row.chr}/{row.ncases_A}_{row.ncontrols_A}_{row.ncases_B}_{row.ncontrols_B}/{row.a_blocks}_{row.b_blocks}_{row.shared_blocks}/window_{row.window}_step_{row.step}_r2_{row.r2}/seed_{row.seed}_tags_{row.tag_A}-{row.tag_B}_hoeffdings.tsv" for row in daf.itertuples()]
    elif filetype == 'gps':
        files = [f"results/gps/simgwas/{reps}_reps/randomised/{row.chr}/{row.ncases_A}_{row.ncontrols_A}_{row.ncases_B}_{row.ncontrols_B}/{row.a_blocks}_{row.b_blocks}_{row.shared_blocks}/window_{row.window}_step_{row.step}_r2_{row.r2}/{draws}_permutations/seed_{row.seed}_tags_{row.tag_A}-{row.tag_B}_gps_pvalue.tsv" for row in daf.itertuples()]
    elif filetype == 'theo':
        files = [f"results/ldsc/simgwas/{reps}_reps/randomised/{row.chr}/{row.ncases_A}_{row.ncontrols_A}_{row.ncases_B}_{row.ncontrols_B}/{row.a_blocks}_{row.b_blocks}_{row.shared_blocks}/theoretical_rg/seed_{row.seed}_{row.tag_A}-{row.tag_B}_theo_rg.tsv" for row in daf.itertuples()]
    elif filetype == 'li_gps':
        files = [f"results/gps/simgwas/{reps}_reps/randomised/{row.chr}/{row.ncases_A}_{row.ncontrols_A}_{row.ncases_B}_{row.ncontrols_B}/{row.a_blocks}_{row.b_blocks}_{row.shared_blocks}/window_{row.window}_step_{row.step}_r2_{row.r2}/seed_{row.seed}_tags_{row.tag_A}-{row.tag_B}_li_gps_pvalue.tsv" for row in daf.itertuples()]
    elif filetype == 'done':
        files = [f"results/simgwas/done/{reps}_reps/randomised/{row.chr}/{row.ncases_A}_{row.ncontrols_A}_{row.ncases_B}_{row.ncontrols_B}/{row.a_blocks}_{row.b_blocks}_{row.shared_blocks}/seed_{row.seed}_tags_{row.tag_A}-{row.tag_B}.done" for row in daf.itertuples()]
    else:
        raise Exception(f"Invalid filetype specified: {filetype}")

    return files
