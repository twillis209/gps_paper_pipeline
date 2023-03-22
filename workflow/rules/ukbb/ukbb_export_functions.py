import pandas as pd

def compile_ukbb_gps_results_into_daf(input_files):
    d = []

    for x in input_files:
        m = re.match(r"results/gps/(?P<snp_set>\w+)/(?P<variant_set>\w+)/window_(?P<window>\d+kb)_step_(?P<step>\d+)_r2_(?P<r2>0_\d+)/(?P<trait_A>\w+)-(?P<trait_B>\w+)_(?P<draws>\d+)_permutations_gps_pvalue.tsv", x)

        try:
            with open(x, 'r') as infile:
                lines = [x.strip() for x in infile.readlines()]

            gps, n, loc, loc_sd, scale, scale_sd, shape, shape_sd, pval = lines[1].split('\t')

            d.append(
                {
                    'snp_set' : m.group('snp_set'),
                    'variant_set' : m.group('variant_set'),
                    'window' : m.group('window'),
                    'step' : m.group('step'),
                    'r2' : m.group('r2').replace('_', '.'),
                    'trait_A' : m.group('trait_A'),
                    'trait_B' : m.group('trait_B'),
                    'draws' : m.group('draws'),
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

def compile_ukbb_li_gps_results_into_daf(input_files):
    d = []

    for x in input_files:
        m = re.match(r"results/gps/(?P<snp_set>\w+)/(?P<variant_set>\w+)/window_(?P<window>\d+kb)_step_(?P<step>\d+)_r2_(?P<r2>0_\d+)/(?P<trait_A>\w+)-(?P<trait_B>\w+)_li_gps_pvalue.tsv", x)

        try:
            with open(x, 'r') as infile:
                lines = [x.strip() for x in infile.readlines()]

            gps, pval = lines[1].split('\t')

            d.append(
                {
                    'snp_set' : m.group('snp_set'),
                    'variant_set' : m.group('variant_set'),
                    'window' : m.group('window'),
                    'step' : m.group('step'),
                    'r2' : m.group('r2').replace('_', '.'),
                    'trait_A' : m.group('trait_A'),
                    'trait_B' : m.group('trait_B'),
                    'gps' : gps,
                    'pval' : pval
                }
            )
        except FileNotFoundError:
            continue

    return pd.DataFrame(d)

def compile_ukbb_ldsc_results_into_daf(input_files):
    d = []

    h2_regex = r"Total \w+ scale h2: (.+)\s+\((.+)\)"
    int_regex = r"Intercept: (.+)\s+\((.+)\)"

    gcov_regex = r"Total \w+ scale gencov: (.+)\s+\((.+)\)"
    gcov_zprod_regex = r"Mean z1\*z2: (.+)"

    for x in input_files:

        m = re.match(r"results/ldsc/rg/ukbb/(?P<snp_set>\w+)/fixed_h2_free_rg_intercept/(?P<trait_A>\w+)-(?P<trait_B>\w+)\.log", x)

        try:
            with open(x, 'r') as infile:
                line = infile.readline()

                # TODO fix these for the null case
                while re.match(h2_regex, line) is None and re.match('ERROR', line) is None:
                    line = infile.readline()

                if re.match('ERROR', line):
                    d.append(
                        {
                            'trait_A' : m.group('trait_A'),
                            'trait_B' : m.group('trait_B'),
                            'snp_set' : m.group('snp_set'),
                            'h2.A.obs.ldsc' : nan,
                            'h2.A.obs.se.ldsc' : nan,
                            'h2.B.obs.ldsc' : nan,
                            'h2.B.obs.se.ldsc' : nan,
                            'gcov.obs.ldsc' : nan,
                            'gcov.obs.se.ldsc' : nan,
                            'rg.ldsc' : nan,
                            'rg.se.ldsc' : nan,
                            'rg.z.ldsc' : nan,
                            'rg.p.ldsc' : nan
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
                            'trait_A' : m.group('trait_A'),
                            'trait_B' : m.group('trait_B'),
                            'snp_set' : m.group('snp_set'),
                            'h2.A.obs.ldsc' : h2_A,
                            'h2.A.obs.se.ldsc' : h2_A_se,
                            'h2.B.obs.ldsc' : h2_B,
                            'h2.B.obs.se.ldsc' : h2_B_se,
                            'gcov.obs.ldsc' : gcov,
                            'gcov.obs.se.ldsc' : gcov_se,
                            'rg.ldsc' : rg,
                            'rg.se.ldsc' : rg_se,
                            'rg.z.ldsc' : rg_z,
                            'rg.p.ldsc' : rg_p
                        }
                    )
        except FileNotFoundError:
            continue
        except AttributeError:
            print(x)

    return pd.DataFrame(d)

def compile_ukbb_sumher_results_into_daf(input_files):
    d = []

    for x in input_files:
        m = re.match(r"results/ldak/ldak-thin/(?P<snp_set>\w+)/rg/(?P<trait_A>\w+)-(?P<trait_B>\w+).cors.full", x)

        with open(x, 'r') as infile:
            line = infile.readline()
            line = infile.readline()

        # Category Trait1_Her SD Trait2_Her SD Both_Coher SD Correlation SD
        _, h2_A, h2_A_se, h2_B, h2_B_se, cov, cov_se, rg, rg_se = line.split()

        rg_z = float(rg)/float(rg_se)

        rg_p = chi2.sf(rg_z**2, df = 1, loc = 0, scale = 1)

        d.append(
            {
                'trait_A' : m.group('trait_A'),
                'trait_B' : m.group('trait_B'),
                'snp_set' : m.group('snp_set'),
                'h2.A.obs.sr' : float(h2_A),
                'h2.A.obs.se.sr' : float(h2_A_se),
                'h2.B.obs.sr' : float(h2_B),
                'h2.B.obs.se.sr' : float(h2_B_se),
                'gcov.obs.sr' : float(cov),
                'gcov.obs.se.sr' : float(cov_se),
                'rg.sr' : float(rg),
                'rg.se.sr' : float(rg_se),
                'rg.z.sr' : rg_z,
                'rg.p.sr' : rg_p
            }
        )

    return pd.DataFrame(d)

def compile_ukbb_hoeffdings_results_into_daf(input_files):
    d = []

    for x in input_files:
        m = re.match(r"results/(?P<snp_set>\w+)/(?P<variant_set>\w+)/window_(?P<window>\d+kb)_step_(?P<step>\d+)_r2_(?P<r2>0_\d+)/(?P<trait_A>\w+)-(?P<trait_B>\w+)_hoeffdings.tsv", x)

        with open(x, 'r') as infile:
            line = infile.readline()
            line = infile.readline()

        trait_A, trait_B, n, Dn, scaled, pvalue = line.split()

        d.append(
            {
                'trait_A' : m.group('trait_A'),
                'trait_B' : m.group('trait_B'),
                'snp_set' : m.group('snp_set'),
                'variant_set' : m.group('variant_set'),
                'window' : m.group('window'),
                'step' : m.group('step'),
                'r2' : m.group('r2'),
                'pval': pvalue
            }
        )

    return pd.DataFrame(d)
