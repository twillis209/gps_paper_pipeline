import pandas as pd
import re
from numpy import nan

exec(open('workflow/rules/ukbb/traits.smk', 'r').read())

def compile_ukbb_ldsc_results_into_daf(input_files):
    d = []

    h2_regex = r"Total \w+ scale h2: (.+)\s+\((.+)\)"
    int_regex = r"Intercept: (.+)\s+\((.+)\)"

    gcov_regex = r"Total \w+ scale gencov: (.+)\s+\((.+)\)"
    gcov_zprod_regex = r"Mean z1\*z2: (.+)"

    for x in input_files:
        print(x)
        m = re.match(r"results/ldsc/rg/ukbb/(?P<snp_set>\w+)/fixed_h2_free_rg_intercept/(?P<trait_A>\w+)-(?P<trait_B>\w+)\.log", x)

        try:
            with open(x, 'r') as infile:
                line = infile.readline()
                print('File opened and first line read')

                # TODO fix these for the null case
                # Probably nature of content has changed so we don't get the regex?
                while re.match(h2_regex, line) is None and re.match('ERROR', line) is None:
                    line = infile.readline()

                print('All lines read')

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

    return pd.DataFrame(d)


ukbb_sans_mhc =  [f"results/ldsc/rg/ukbb/ukbb_sans_mhc/fixed_h2_free_rg_intercept/{trait_pair}.log" for trait_pair in ukbb_trait_pairs]
ukbb_with_mhc =  [f"results/ldsc/rg/ukbb/ukbb_with_mhc/fixed_h2_free_rg_intercept/{trait_pair}.log" for trait_pair in ukbb_trait_pairs]

sans_mhc_daf = compile_ukbb_ldsc_results_into_daf(ukbb_sans_mhc)
with_mhc_daf = compile_ukbb_ldsc_results_into_daf(ukbb_with_mhc)
