from scipy.stats import chi2
import os
import pandas as pd
import re
from numpy import nan

def get_all_block_files(sample_sizes, effect):

    block_files = []

    seed_label = f"{effect}_seed"

    for row in block_daf.itertuples():
        for sample_size in sample_sizes:
            if effect == 'tiny':
                block_files.append(f"results/simgwas/simulated_sum_stats/block_sum_stats/400_reps/{effect}/2_cv/{sample_size[0]}_{sample_size[1]}/chr{row.chr}/block_{row.block}_seed_{getattr(row, seed_label)}_sum_stats.tsv.gz")
            else:
                block_files.append(f"results/simgwas/simulated_sum_stats/block_sum_stats/400_reps/{effect}/1_cv/{sample_size[0]}_{sample_size[1]}/chr{row.chr}/block_{row.block}_seed_{getattr(row, seed_label)}_sum_stats.tsv.gz")


    return block_files

def get_all_simulation_done_files(simulation_pars_file, reps, subset = None):
    daf = pd.read_csv(simulation_pars_file, sep = '\t')

    if subset:
        daf = daf.query('a_blocks == @subset & b_blocks == @subset')

    return [f"results/simgwas/done/{reps}_reps/randomised/{row.ncases_A}_{row.ncontrols_A}_{row.ncases_B}_{row.ncontrols_B}/{row.a_blocks}_{row.b_blocks}_{row.shared_blocks}/seed_{row.seed}_tags_{row.tag_A}-{row.tag_B}.done" for row in daf.itertuples()]

def get_all_randomised_block_files(simulation_pars_file, reps):
    daf = pd.read_csv(simulation_pars_file, sep = '\t')

    a_block_files =  [f"results/simgwas/simulated_sum_stats/whole_genome_sum_stats/{reps}_reps/randomised/{row.ncases_A}_{row.ncontrols_A}_{row.ncases_B}_{row.ncontrols_B}/{row.a_blocks}_{row.b_blocks}_{row.shared_blocks}/seed_{row.seed}_sum_stats_A_tags_{row.tag_A}-{row.tag_B}_files.txt" for row in daf.itertuples()]
    b_block_files = [f"results/simgwas/simulated_sum_stats/whole_genome_sum_stats/{reps}_reps/randomised/{row.ncases_A}_{row.ncontrols_A}_{row.ncases_B}_{row.ncontrols_B}/{row.a_blocks}_{row.b_blocks}_{row.shared_blocks}/seed_{row.seed}_sum_stats_B_tags_{row.tag_A}-{row.tag_B}_files.txt" for row in daf.itertuples()]

    return a_block_files+b_block_files

def get_test_files(simulation_pars_file, reps, filetype, subset = None, draws = "3000", window = "1000kb", step = 50, r2 = 0.2):
    daf = pd.read_csv(simulation_pars_file, sep = '\t')

    str_r2 = str(r2).replace('.', '_')

    if subset:
        daf = daf.query(subset)

    if filetype == 'ldsc':
        files = [f"results/ldsc/simgwas/{reps}_reps/randomised/{row.ncases_A}_{row.ncontrols_A}_{row.ncases_B}_{row.ncontrols_B}/{row.a_blocks}_{row.b_blocks}_{row.shared_blocks}/rg/fixed_h2_free_rg_intercept/seed_{row.seed}_tags_{row.tag_A}-{row.tag_B}.log" for row in daf.itertuples()]
    elif filetype == 'sumher':
        files = [f"results/ldak/ldak-thin/simgwas/{reps}_reps/randomised/rg/{row.ncases_A}_{row.ncontrols_A}_{row.ncases_B}_{row.ncontrols_B}/{row.a_blocks}_{row.b_blocks}_{row.shared_blocks}/seed_{row.seed}_tags_{row.tag_A}-{row.tag_B}.cors" for row in daf.itertuples()]
    elif filetype == 'hoeffdings':
        files = [f"results/hoeffdings/simgwas/{reps}_reps/randomised/{row.ncases_A}_{row.ncontrols_A}_{row.ncases_B}_{row.ncontrols_B}/{row.a_blocks}_{row.b_blocks}_{row.shared_blocks}/window_{window}_step_{step}_r2_{str_r2}/seed_{row.seed}_tags_{row.tag_A}-{row.tag_B}_hoeffdings.tsv" for row in daf.itertuples()]
    elif filetype == 'gps':
        files = [f"results/gps/simgwas/{reps}_reps/randomised/{row.ncases_A}_{row.ncontrols_A}_{row.ncases_B}_{row.ncontrols_B}/{row.a_blocks}_{row.b_blocks}_{row.shared_blocks}/window_{window}_step_{step}_r2_{str_r2}/{draws}_permutations/seed_{row.seed}_tags_{row.tag_A}-{row.tag_B}_gps_pvalue.tsv" for row in daf.itertuples()]
    elif filetype == 'theo':
        files = [f"results/ldsc/simgwas/{reps}_reps/randomised/{row.ncases_A}_{row.ncontrols_A}_{row.ncases_B}_{row.ncontrols_B}/{row.a_blocks}_{row.b_blocks}_{row.shared_blocks}/theoretical_rg/seed_{row.seed}_{row.tag_A}-{row.tag_B}_theo_rg.tsv" for row in daf.itertuples()]
    elif filetype == 'li_gps':
        files = [f"results/gps/simgwas/{reps}_reps/randomised/{row.ncases_A}_{row.ncontrols_A}_{row.ncases_B}_{row.ncontrols_B}/{row.a_blocks}_{row.b_blocks}_{row.shared_blocks}/window_{window}_step_{step}_r2_{str_r2}/seed_{row.seed}_tags_{row.tag_A}-{row.tag_B}_li_gps_pvalue.tsv" for row in daf.itertuples()]
    elif filetype == 'mean_stat':
        files = [f"results/gps/simgwas/{reps}_reps/randomised/{row.ncases_A}_{row.ncontrols_A}_{row.ncases_B}_{row.ncontrols_B}/{row.a_blocks}_{row.b_blocks}_{row.shared_blocks}/window_{window}_step_{step}_r2_{str_r2}/seed_{row.seed}_tags_{row.tag_A}-{row.tag_B}/mean_stat/{draws}_permutations/pp_pert_1_pvalue.tsv" for row in daf.itertuples()]
    elif filetype == 'done':
        files = [f"results/simgwas/done/{reps}_reps/randomised/{row.ncases_A}_{row.ncontrols_A}_{row.ncases_B}_{row.ncontrols_B}/{row.a_blocks}_{row.b_blocks}_{row.shared_blocks}/seed_{row.seed}_tags_{row.tag_A}-{row.tag_B}.done" for row in daf.itertuples()]
    else:
        raise Exception(f"Invalid filetype specified: {filetype}")

    return files

def compile_sumher_results_into_daf(input_files):
    d = []

    for x in input_files:

        m = re.match(r"results/ldak/ldak-thin/simgwas/(?P<no_reps>\d+)_reps/randomised/rg(/chr\d+)?/(?P<ncases_A>\d+)_(?P<ncontrols_A>\d+)_(?P<ncases_B>\d+)_(?P<ncontrols_B>\d+)/(?P<a_blocks>[\w-]+)_(?P<b_blocks>[\w-]+)_(?P<shared_blocks>[\w-]+)/seed_(?P<seed>\w+)_tags_(?P<tag_a>\d+)-(?P<tag_b>\d+)\.cors", x)

        try:
            with open(x, 'r') as infile:

                line = infile.readline()

                while re.match("^Her1_All", line) is None:
                    line = infile.readline()

                h2_A, h2_A_se = re.match("Her1_All (-?\d+\.\d+|-?nan) (\d+\.\d+|nan)", line).groups()

                line = infile.readline()

                h2_B, h2_B_se = re.match("Her2_All (-?\d+\.\d+|-?nan) (\d+\.\d+|nan)", line).groups()

                line = infile.readline()

                cov, cov_se = re.match("Coher_All (-?\d+\.\d+|-?nan) (\d+\.\d+|nan)", line).groups()

                line = infile.readline()

                rg, rg_se = re.match("Cor_All (-?\d+\.\d+|-?nan) (\d+\.\d+|nan)", line).groups()

            rg_z = float(rg)/float(rg_se)

            rg_p = chi2.sf(rg_z**2, df = 1, loc = 0, scale = 1)

            d.append(
                {
                    'ncases.A' : m.group('ncases_A'),
                    'ncontrols.A' : m.group('ncontrols_A'),
                    'ncases.B' : m.group('ncases_B'),
                    'ncontrols.B' : m.group('ncontrols_B'),
                    'blocks.A' : m.group('a_blocks'),
                    'blocks.B' : m.group('b_blocks'),
                    'shared_blocks' : m.group('shared_blocks'),
                    'tag_pair' : f"{m.group('tag_a')}-{m.group('tag_b')}",
                    'seed' : f"{m.group('seed')}",
                    'h2.A' : float(h2_A),
                    'h2.A.se' : float(h2_A_se),
                    'h2.B' : float(h2_B),
                    'h2.B.se' : float(h2_B_se),
                    'gcov' : float(cov),
                    'gcov.se' : float(cov_se),
                    'rg' : float(rg),
                    'rg.se' : float(rg_se),
                    'rg.z' : float(rg_z),
                    'rg.p' : float(rg_p)
                }
            )
        except FileNotFoundError:
            continue

    return pd.DataFrame(d)

def compile_ldsc_results_into_daf(input_files):
    d = []

    h2_regex = r"Total Liability scale h2: (.+)\s+\((.+)\)"
    int_regex = r"Intercept: (.+)\s+\((.+)\)"

    gcov_regex = r"Total Liability scale gencov: (.+)\s+\((.+)\)"
    gcov_zprod_regex = r"Mean z1\*z2: (.+)"

    for x in input_files:

        m = re.match(r"results/ldsc/simgwas/(?P<no_reps>\d+)_reps/randomised/(chr\d+/)?(?P<ncases_A>\d+)_(?P<ncontrols_A>\d+)_(?P<ncases_B>\d+)_(?P<ncontrols_B>\d+)/(?P<a_blocks>[\w-]+)_(?P<b_blocks>[\w-]+)_(?P<shared_blocks>[\w-]+)/rg/fixed_h2_free_rg_intercept/seed_(?P<seed>\w+)_tags_(?P<tag_a>\d+)-(?P<tag_b>\d+)\.log", x)

        try:
            with open(x, 'r') as infile:
                line = infile.readline()

                # TODO fix these for the null case
                while re.match(h2_regex, line) is None and re.match('ERROR', line) is None:
                    line = infile.readline()

                if re.match('ERROR', line):
                    d.append(
                        {
                            'ncases.A' : m.group('ncases_A'),
                            'ncontrols.A' : m.group('ncontrols_A'),
                            'ncases.B' : m.group('ncases_B'),
                            'ncontrols.B' : m.group('ncontrols_B'),
                            #'odds_ratio.A': odds_ratios_A,
                            #'odds_ratio.B': odds_ratios_B,
                            'blocks.A' : m.group('a_blocks'),
                            'blocks.B' : m.group('b_blocks'),
                            'shared_blocks' : m.group('shared_blocks'),
                            'tag_pair' : f"{m.group('tag_a')}-{m.group('tag_b')}",
                            'seed' : f"{m.group('seed')}",
                            'h2.A' : nan,
                            'h2.A.se' : nan,
                            'h2.B' : nan,
                            'h2.B.se' : nan,
                            'gcov' : nan,
                            'gcov.se' : nan,
                            'rg' : nan,
                            'rg.se' : nan,
                            'rg.z' : nan,
                            'rg.p' : nan
                        }
                    )

                else:
                    h2_match_A = re.match(h2_regex, line)
                    h2_A = float(h2_match_A.group(1))
                    h2_A_se = float(h2_match_A.group(2))

                    line = infile.readline()
                    line = infile.readline()
                    line = infile.readline()

                    h2_int_A_match = re.match(int_regex, line)

                    if h2_int_A_match:
                        h2_int_A = float(h2_int_A_match.group(1))
                        h2_int_A_se = float(h2_int_A_match.group(2))
                    elif 'constrained to 1.' in line:
                        h2_int_A = 1.0
                        h2_int_A_se = nan
                    else:
                        raise Exception("No match for h2_B int_regex")

                    while re.match(h2_regex, line) is None:
                        line = infile.readline()

                    h2_match_B = re.match(h2_regex, line)
                    h2_B = float(h2_match_B.group(1))
                    h2_B_se = float(h2_match_B.group(2))

                    line = infile.readline()
                    line = infile.readline()
                    line = infile.readline()

                    h2_int_B_match = re.match(int_regex, line)

                    if h2_int_B_match:
                            h2_int_B = float(h2_int_B_match.group(1))
                            h2_int_B_se = float(h2_int_B_match.group(2))
                    elif 'constrained to 1.' in line:
                            h2_int_B = 1.0
                            h2_int_B_se = nan
                    else:
                            raise Exception("No match for h2_A int_regex")

                    while re.match(gcov_regex, line) is None:
                        line = infile.readline()

                    gcov_match = re.match(gcov_regex, line)
                    gcov = float(gcov_match.group(1))
                    gcov_se = float(gcov_match.group(2))

                    line = infile.readline()

                    gcov_zprod_match = re.match(gcov_zprod_regex, line)
                    gcov_zprod = float(gcov_zprod_match.group(1))

                    line = infile.readline()

                    gcov_int_match = re.match(int_regex, line)

                    if gcov_int_match:
                        gcov_int = float(gcov_int_match.group(1))
                        gcov_int_se = float(gcov_int_match.group(2))
                    elif 'constrained to 0.' in line:
                        gcov_int = 0.0
                        gcov_int_se = nan
                    else:
                        raise Exception("No match for gcov_int_regex")

                    line = infile.readline()

                    while re.match("^p1\s", line) is None:
                        line = infile.readline()

                    line = infile.readline()

                    rg, rg_se, rg_z, rg_p = [float(z) if z != 'NA' else nan for z in line.split()[2:6]]

                d.append(
                    {
                        'ncases.A' : m.group('ncases_A'),
                        'ncontrols.A' : m.group('ncontrols_A'),
                        'ncases.B' : m.group('ncases_B'),
                        'ncontrols.B' : m.group('ncontrols_B'),
                        #'odds_ratio.A': odds_ratios_A,
                        #'odds_ratio.B': odds_ratios_B,
                        'blocks.A' : m.group('a_blocks'),
                        'blocks.B' : m.group('b_blocks'),
                        'shared_blocks' : m.group('shared_blocks'),
                        'tag_pair' : f"{m.group('tag_a')}-{m.group('tag_b')}",
                        'seed' : f"{m.group('seed')}",
                        'h2.A' : h2_A,
                        'h2.A.se' : h2_A_se,
                        'h2.B' : h2_B,
                        'h2.B.se' : h2_B_se,
                        'gcov' : gcov,
                        'gcov.se' : gcov_se,
                        'rg' : rg,
                        'rg.se' : rg_se,
                        'rg.z' : rg_z,
                        'rg.p' : rg_p
                    }
                )
        except FileNotFoundError:
            continue
        except AttributeError:
            print(x)


    return pd.DataFrame(d)

def compile_hoeffdings_results_into_daf(input_files):
    d = []

    for x in input_files:
        m = re.match(r"results/hoeffdings/simgwas/(?P<no_reps>\d+)_reps/randomised(/chr\d+)?/(?P<ncases_A>\d+)_(?P<ncontrols_A>\d+)_(?P<ncases_B>\d+)_(?P<ncontrols_B>\d+)/(?P<a_blocks>[\w-]+)_(?P<b_blocks>[\w-]+)_(?P<shared_blocks>[\w-]+)/window_(?P<window>\d+kb)_step_(?P<step>\d+)_r2_(?P<r2>\d+_\d+)/seed_(?P<seed>\w+)_tags_(?P<tag_a>\d+)-(?P<tag_b>\d+)_hoeffdings\.tsv", x)

        try:
            with open(x, 'r') as infile:
                lines = [x.strip() for x in infile.readlines()]

            # NB: The 'n' here is no. of SNPs, not no. of permutations as with GPS
            _, _, n, Dn, scaled, pval = lines[1].split('\t')

            d.append(
                {
                    'ncases.A' : m.group('ncases_A'),
                    'ncontrols.A' : m.group('ncontrols_A'),
                    'ncases.B' : m.group('ncases_B'),
                    'ncontrols.B' : m.group('ncontrols_B'),
                    'blocks.A' : m.group('a_blocks'),
                    'blocks.B' : m.group('b_blocks'),
                    'shared_blocks' : m.group('shared_blocks'),
                    'tag_pair' : f"{m.group('tag_a')}-{m.group('tag_b')}",
                    'seed' : f"{m.group('seed')}",
                    'window': m.group('window'),
                    'step': m.group('step'),
                    'r2': float(m.group('r2').replace('_', '.')),
                    'hoeff.p' : pval
                }
            )

        except ValueError:
            _, _, n, Dn, scaled = lines[1].split('\t')
            pval = nan

            d.append(
                {
                    'ncases.A' : m.group('ncases_A'),
                    'ncontrols.A' : m.group('ncontrols_A'),
                    'ncases.B' : m.group('ncases_B'),
                    'ncontrols.B' : m.group('ncontrols_B'),
                    'blocks.A' : m.group('a_blocks'),
                    'blocks.B' : m.group('b_blocks'),
                    'shared_blocks' : m.group('shared_blocks'),
                    'tag_pair' : f"{m.group('tag_a')}-{m.group('tag_b')}",
                    'seed' : f"{m.group('seed')}",
                    'window': m.group('window'),
                    'step': m.group('step'),
                    'r2': float(m.group('r2').replace('_', '.')),
                    'hoeff.p' : pval
                }
            )
            continue
        except FileNotFoundError:
            continue

    return pd.DataFrame(d)

def compile_gps_results_into_daf(input_files):
    d = []

    for x in input_files:

        m = re.match(r"results/gps/simgwas/(?P<no_reps>\d+)_reps/randomised/(chr\d+/)?(?P<ncases_A>\d+)_(?P<ncontrols_A>\d+)_(?P<ncases_B>\d+)_(?P<ncontrols_B>\d+)/(?P<a_blocks>[\w-]+)_(?P<b_blocks>[\w-]+)_(?P<shared_blocks>[\w-]+)/window_(?P<window>\d+kb)_step_(?P<step>\d+)_r2_(?P<r2>\d+_\d+)/(?P<draws>\d+)_permutations/seed_(?P<seed>\w+)_tags_(?P<tag_a>\d+)-(?P<tag_b>\d+)_gps_pvalue\.tsv", x)

        try:
            with open(x, 'r') as infile:
                lines = [x.strip() for x in infile.readlines()]

            gps, n, loc, loc_sd, scale, scale_sd, shape, shape_sd, pval = lines[1].split('\t')

            d.append(
                {
                    'ncases.A' : m.group('ncases_A'),
                    'ncontrols.A' : m.group('ncontrols_A'),
                    'ncases.B' : m.group('ncases_B'),
                    'ncontrols.B' : m.group('ncontrols_B'),
                    'blocks.A' : m.group('a_blocks'),
                    'blocks.B' : m.group('b_blocks'),
                    'shared_blocks' : m.group('shared_blocks'),
                    'tag_pair' : f"{m.group('tag_a')}-{m.group('tag_b')}",
                    'seed' : f"{m.group('seed')}",
                    'window': m.group('window'),
                    'step': m.group('step'),
                    'r2': float(m.group('r2').replace('_', '.')),
                    'gps' : gps,
                    'n' : n,
                    'loc' : loc,
                    'loc.sd' : loc_sd,
                    'scale' : scale,
                    'scale.sd' : scale_sd,
                    'shape' : shape,
                    'shape.sd' : shape_sd,
                    'pval' : pval,
                }
            )
        except FileNotFoundError:
            continue

    return pd.DataFrame(d)

def compile_li_gps_results_into_daf(input_files):
    d = []

    for x in input_files:

        m = re.match(r"results/gps/simgwas/(?P<no_reps>\d+)_reps/randomised/(chr\d+/)?(?P<ncases_A>\d+)_(?P<ncontrols_A>\d+)_(?P<ncases_B>\d+)_(?P<ncontrols_B>\d+)/(?P<a_blocks>[\w-]+)_(?P<b_blocks>[\w-]+)_(?P<shared_blocks>[\w-]+)/window_(?P<window>\d+kb)_step_(?P<step>\d+)_r2_(?P<r2>\d+_\d+)/seed_(?P<seed>\w+)_tags_(?P<tag_a>\d+)-(?P<tag_b>\d+)_li_gps_pvalue\.tsv", x)

        try:
            with open(x, 'r') as infile:
                lines = [x.strip() for x in infile.readlines()]

            gps, pval = lines[1].split('\t')

            d.append(
                {
                    'ncases.A' : m.group('ncases_A'),
                    'ncontrols.A' : m.group('ncontrols_A'),
                    'ncases.B' : m.group('ncases_B'),
                    'ncontrols.B' : m.group('ncontrols_B'),
                    'blocks.A' : m.group('a_blocks'),
                    'blocks.B' : m.group('b_blocks'),
                    'shared_blocks' : m.group('shared_blocks'),
                    'tag_pair' : f"{m.group('tag_a')}-{m.group('tag_b')}",
                    'seed' : f"{m.group('seed')}",
                    'window': m.group('window'),
                    'step': m.group('step'),
                    'r2': float(m.group('r2').replace('_', '.')),
                    'gps' : gps,
                    'pval' : pval
                }
            )
        except FileNotFoundError:
            continue

    return pd.DataFrame(d)

def compile_theo_results_into_daf(input_files):
    d = []

    for x in input_files:

        m = re.match(r"results/ldsc/simgwas/(?P<no_reps>\d+)_reps/randomised/(chr\d+/)?(?P<ncases_A>\d+)_(?P<ncontrols_A>\d+)_(?P<ncases_B>\d+)_(?P<ncontrols_B>\d+)/(?P<a_blocks>[\w-]+)_(?P<b_blocks>[\w-]+)_(?P<shared_blocks>[\w-]+)/theoretical_rg/seed_(?P<seed>\w+)_(?P<tag_a>\d+)-(?P<tag_b>\d+)_theo_rg\.tsv", x)

        with open(x, 'r') as infile:
            lines = [x.strip() for x in infile.readlines()]

        _, _, _, _, _, h2_theo_obs_A, h2_theo_obs_B, h2_theo_liab_A, h2_theo_liab_B, V_A_A, V_A_B, C_A_AB, r_A_AB  = lines[1].split('\t')

        d.append(
            {
                'ncases.A' : m.group('ncases_A'),
                'ncontrols.A' : m.group('ncontrols_A'),
                'ncases.B' : m.group('ncases_B'),
                'ncontrols.B' : m.group('ncontrols_B'),
                'blocks.A' : m.group('a_blocks'),
                'blocks.B' : m.group('b_blocks'),
                'shared_blocks' : m.group('shared_blocks'),
                'tag_pair' : f"{m.group('tag_a')}-{m.group('tag_b')}",
                'seed' : f"{m.group('seed')}",
                "h2.theo.obs.A" : float(h2_theo_obs_A),
                "h2.theo.obs.B" : float(h2_theo_obs_B),
                "h2.theo.liab.A" : float(h2_theo_liab_A),
                "h2.theo.liab.B" : float(h2_theo_liab_B),
                "V_A.A" : float(V_A_A),
                "V_A.B" : float(V_A_B),
                "C_A.AB" : float(C_A_AB),
                "r_A.AB" : float(r_A_AB)
            }
        )

    return pd.DataFrame(d)
