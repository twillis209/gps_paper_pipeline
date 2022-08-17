from scipy.stats import chi2
import os
import pandas as pd
import re
from numpy import nan

def get_all_block_files(sample_sizes):

    block_files = []

    for row in block_daf.itertuples():
        for sample_size in sample_sizes:
            block_files.append(f"results/simgwas/simulated_sum_stats/block_sum_stats/400_reps/null/{sample_size[0]}_{sample_size[1]}/chr{row.chr}/block_{row.block}_seed_{row.null_seed}_sum_stats.tsv.gz")
            block_files.append(f"results/simgwas/simulated_sum_stats/block_sum_stats/400_reps/small/{sample_size[0]}_{sample_size[1]}/chr{row.chr}/block_{row.block}_seed_{row.small_seed}_sum_stats.tsv.gz")
            block_files.append(f"results/simgwas/simulated_sum_stats/block_sum_stats/400_reps/medium/{sample_size[0]}_{sample_size[1]}/chr{row.chr}/block_{row.block}_seed_{row.medium_seed}_sum_stats.tsv.gz")

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

def parse_effect_token_to_odds_ratios(token):
    return ','.join([str(odds_ratio_dict[re.match('[smlhv]', x).group()]) for x in token.split('-')])

def compile_theoretical_rg_results(input, output):
    with open(output[0], 'w') as outfile:
        outfile.write("ncases.A\tncontrols.A\tncases.B\tncontrols.B\tseed\ttag_pair\todds_ratio.A\todds_ratio.B\tblocks.A\tblocks.B\tshared_blocks\th2.theo.obs.A\th2.theo.obs.B\th2.theo.liab.A\th2.theo.liab.B\tV_A.A\tV_A.B\tC_A.AB\tr_A.AB\n")

        for x in input:
            with open(x, 'r') as infile:
                head, tail = os.path.split(x)

                ncases_A, ncontrols_A, ncases_B, ncontrols_B = re.match("results/ldsc/rg/whole_genome/randomised/theoretical_rg/(\d+)_(\d+)_(\d+)_(\d+)", head).groups()


                effect_blocks_A, effect_blocks_B, shared_effect_blocks, seed, tag_pair = re.match("([smlhv\d-]+)_([smlhv\d-]+)_([smlhv\d-]+)_seed_(\d+)_(\w{2})", tail).groups()

                odds_ratios_A = parse_effect_token_to_odds_ratios(effect_blocks_A)
                odds_ratios_B = parse_effect_token_to_odds_ratios(effect_blocks_B)

                lines = [y.strip() for y in infile.readlines()]

                h2_theo_obs_A, h2_theo_obs_B, h2_theo_liab_A, h2_theo_liab_B, V_A_A, V_A_B, C_A_AB, r_A_AB = lines[1].split('\t')[5:]
                outfile.write(f"{ncases_A}\t{ncontrols_A}\t{ncases_B}\t{ncontrols_B}\t{seed}\t{tag_pair}\t{odds_ratios_A}\t{odds_ratios_B}\t{effect_blocks_A}\t{effect_blocks_B}\t{shared_effect_blocks}\t{h2_theo_obs_A}\t{h2_theo_obs_B}\t{h2_theo_liab_A}\t{h2_theo_liab_B}\t{V_A_A}\t{V_A_B}\t{C_A_AB}\t{r_A_AB}\n")

    return

def get_existing_test_files(simulation_pars_file, reps, filetype, subset = None):
    daf = pd.read_csv(simulation_pars_file, sep = '\t')

    if subset:
        daf = daf.query('a_blocks == @subset & b_blocks == @subset')

    if filetype == 'ldsc':
        files = [f"results/ldsc/simgwas/{reps}_reps/randomised/{row.ncases_A}_{row.ncontrols_A}_{row.ncases_B}_{row.ncontrols_B}/{row.a_blocks}_{row.b_blocks}_{row.shared_blocks}/rg/fixed_h2_free_rg_intercept/seed_{row.seed}_tags_{row.tag_A}-{row.tag_B}.log" for row in daf.itertuples()]
    elif filetype == 'sumher':
        files = [f"results/ldak/ldak-thin/simgwas/{reps}_reps/randomised/rg/{row.ncases_A}_{row.ncontrols_A}_{row.ncases_B}_{row.ncontrols_B}/{row.a_blocks}_{row.b_blocks}_{row.shared_blocks}/seed_{row.seed}_tags_{row.tag_A}-{row.tag_B}.cors" for row in daf.itertuples()]
    elif filetype == 'hoeffdings':
        files = [f"results/hoeffdings/simgwas/{reps}_reps/randomised/{row.ncases_A}_{row.ncontrols_A}_{row.ncases_B}_{row.ncontrols_B}/{row.a_blocks}_{row.b_blocks}_{row.shared_blocks}/window_1000kb_step_50/seed_{row.seed}_tags_{row.tag_A}-{row.tag_B}_hoeffdings.tsv" for row in daf.itertuples()]
    elif filetype == 'gps':
        files = [f"results/gps/simgwas/{reps}_reps/randomised/{row.ncases_A}_{row.ncontrols_A}_{row.ncases_B}_{row.ncontrols_B}/{row.a_blocks}_{row.b_blocks}_{row.shared_blocks}/window_1000kb_step_50/3000_permutations/seed_{row.seed}_tags_{row.tag_A}-{row.tag_B}_gps_pvalue.tsv" for row in daf.itertuples()]
    else:
        raise Exception(f"Invalid filetype specified: {filetype}")

    return [x for x in files if os.path.exists(x)]

def compile_sumher_results_into_daf(input_files):
    daf = pd.DataFrame()

    d = []

    for x in input_files:

        m = re.match(r"results/ldak/ldak-thin/simgwas/(?P<no_reps>\d+)_reps/randomised/rg/(?P<ncases_A>\d+)_(?P<ncontrols_A>\d+)_(?P<ncases_B>\d+)_(?P<ncontrols_B>\d+)/(?P<a_blocks>[\w]+)_(?P<b_blocks>[\w]+)_(?P<shared_blocks>[\w]+)/seed_(?P<seed>\w+)_tags_(?P<tag_a>\d+)-(?P<tag_b>\d+)\.cors", x)

        #odds_ratios_A = parse_effect_token_to_odds_ratios(m.group('a_blocks'))
        #odds_ratios_B = parse_effect_token_to_odds_ratios(m.group('b_blocks'))

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
                #'odds_ratio.A': odds_ratios_A,
                #'odds_ratio.B': odds_ratios_B,
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

    return pd.DataFrame(d)

def compile_ldsc_results_into_daf(input_files):
    daf = pd.DataFrame()

    d = []

    h2_regex = r"Total Liability scale h2: (.+)\s+\((.+)\)"
    int_regex = r"Intercept: (.+)\s+\((.+)\)"

    gcov_regex = r"Total Liability scale gencov: (.+)\s+\((.+)\)"
    gcov_zprod_regex = r"Mean z1\*z2: (.+)"

    for x in input_files:

        m = re.match(r"results/ldsc/simgwas/(?P<no_reps>\d+)_reps/randomised/(?P<ncases_A>\d+)_(?P<ncontrols_A>\d+)_(?P<ncases_B>\d+)_(?P<ncontrols_B>\d+)/(?P<a_blocks>[\w]+)_(?P<b_blocks>[\w]+)_(?P<shared_blocks>[\w]+)/rg/fixed_h2_free_rg_intercept/seed_(?P<seed>\w+)_tags_(?P<tag_a>\d+)-(?P<tag_b>\d+)\.log", x)

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

    return pd.DataFrame(d)

def compile_hoeffdings_results_into_daf(input_files):
    daf = pd.DataFrame()

    d = []

    for x in input_files:

        m = re.match(r"results/hoeffdings/simgwas/(?P<no_reps>\d+)_reps/randomised/(?P<ncases_A>\d+)_(?P<ncontrols_A>\d+)_(?P<ncases_B>\d+)_(?P<ncontrols_B>\d+)/(?P<a_blocks>[\w]+)_(?P<b_blocks>[\w]+)_(?P<shared_blocks>[\w]+)/window_1000kb_step_50/seed_(?P<seed>\w+)_tags_(?P<tag_a>\d+)-(?P<tag_b>\d+)_hoeffdings\.tsv", x)

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
                'hoeff.p' : pval
            }
        )

    return pd.DataFrame(d)


def compile_gps_results_into_daf(input_files):
    daf = pd.DataFrame()

    d = []

    for x in input_files:

        m = re.match(r"results/gps/simgwas/(?P<no_reps>\d+)_reps/randomised/(?P<ncases_A>\d+)_(?P<ncontrols_A>\d+)_(?P<ncases_B>\d+)_(?P<ncontrols_B>\d+)/(?P<a_blocks>[\w]+)_(?P<b_blocks>[\w]+)_(?P<shared_blocks>[\w]+)/window_1000kb_step_50/3000_permutations/seed_(?P<seed>\w+)_tags_(?P<tag_a>\d+)-(?P<tag_b>\d+)_gps_pvalue\.tsv", x)

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

    return pd.DataFrame(d)
