import re
import os
from numpy import nan

tag_pairs = [tags[i]+tags[i+1] for i in range(0, 19, 2)]

odds_ratio_dict = {"s": 1.05, "m": 1.2, 'l': 1.4, 'v': 2, 'r': 'random', 'n' : 1, 'i' : 1.1}

medium_effect_tuples = [f"m50_m50_m{x}" for x in [0, 10, 20, 30, 40, 50]]+[f"m25_m25_m{x}" for x in [0, 5, 10, 15, 20, 25]]

small_effect_tuples = [f"s400_s400_s{x}" for x in range(0, 101, 10)]

sample_sizes = [1000, 5000, 10000, 50000, 100000, 250000]

medium_effect_rg_estimate_files = [f"results/ldsc/rg/whole_genome/randomised/{size}_{size}_{size}_{size}/{effect_tuple}_seed_%d_{tag_pair}.log" for size in sample_sizes for tag_pair in tag_pairs for effect_tuple in medium_effect_tuples]

medium_effect_gps_files = [f"results/gps/simgwas/randomised/window_1000kb_step_50/{size}_{size}_{size}_{size}/3000_permutations/{effect_tuple}_seed_%d_tags_{tag_pair}_gps_pvalue.tsv" for size in sample_sizes for tag_pair in tag_pairs for effect_tuple in medium_effect_tuples]

medium_effect_hoeffdings_files = [f"results/hoeffdings/simgwas/randomised/window_1000kb_step_50/{size}_{size}_{size}_{size}/{effect_tuple}_seed_%d_tags_{tag_pair}_hoeffdings.tsv" for size in sample_sizes for tag_pair in tag_pairs for effect_tuple in medium_effect_tuples]

small_effect_rg_estimate_files = [f"results/ldsc/rg/whole_genome/randomised/{size}_{size}_{size}_{size}/{effect_tuple}_seed_%d_{tag_pair}.log" for size in sample_sizes for tag_pair in tag_pairs for effect_tuple in small_effect_tuples]

small_effect_gps_files = [f"results/gps/simgwas/randomised/window_1000kb_step_50/{size}_{size}_{size}_{size}/3000_permutations/{effect_tuple}_seed_%d_tags_{tag_pair}_gps_pvalue.tsv" for size in sample_sizes for tag_pair in tag_pairs for effect_tuple in small_effect_tuples]

small_effect_hoeffdings_files = [f"results/hoeffdings/simgwas/randomised/window_1000kb_step_50/{size}_{size}_{size}_{size}/{effect_tuple}_seed_%d_tags_{tag_pair}_hoeffdings.tsv" for size in sample_sizes for tag_pair in tag_pairs for effect_tuple in small_effect_tuples]

rule compile_medium_rg_estimates:
    input:
        [x[0] % x[1] for x in zip(medium_effect_rg_estimate_files, range(100, 100+len(medium_effect_rg_estimate_files)))]
    output:
        "results/ldsc/rg/whole_genome/randomised/compiled_m_rg_estimates.tsv"
    run:
        with open(output[0], 'w') as outfile:
            outfile.write("ncases.A\tncontrols.A\tncases.B\tncontrols.B\todds_ratio.A\todds_ratio.B\tblocks.A\tblocks.B\tno_shared_blocks\ttag_pair\th2.A\th2.A.se\th2.B\th2.B.se\tgcov\tgcov.se\tgcov.int\tgcov.int.se\tgcov.zprod\trg\trg.se\trg.p\n")
            for x in input:
                head, tail = os.path.split(x)
                head_res = re.match("results/ldsc/rg/whole_genome/randomised/(\d+)_(\d+)_(\d+)_(\d+)", head)
                ncases_A = int(head_res.group(1))
                ncontrols_A = int(head_res.group(2))
                ncases_B = int(head_res.group(3))
                ncontrols_B = int(head_res.group(4))

                odds_ratio_A = 1.2
                odds_ratio_B = 1.2

                tail_res = re.match("m(\d+)_m(\d+)_m(\d+)_seed_(\d+)_(\w{2})", tail)

                effect_blocks_A, effect_blocks_B, shared_effect_blocks, seed, tag_pair = tail_res.groups()

                h2_regex = r"Total Liability scale h2: (.+)\s+\((.+)\)"

                gcov_regex = r"Total Liability scale gencov: (.+)\s+\((.+)\)"
                gcov_zprod_regex = r"Mean z1\*z2: (.+)"
                gcov_int_regex = r"Intercept: (.+)\s+\((.+)\)"

                with open(x, 'r') as infile:
                    line = infile.readline()

                    # TODO fix these for the null case
                    while re.match(h2_regex, line) is None:
                        line = infile.readline()

                    h2_match_A = re.match(h2_regex, line)
                    h2_A = float(h2_match_A.group(1))
                    h2_A_se = float(h2_match_A.group(2))

                    line = infile.readline()

                    while re.match(h2_regex, line) is None:
                        line = infile.readline()

                    h2_match_B = re.match(h2_regex, line)
                    h2_B = float(h2_match_B.group(1))
                    h2_B_se = float(h2_match_B.group(2))

                    line = infile.readline()

                    gcov_regex = r"Total Liability scale gencov: (.+)\s+\((.+)\)"
                    gcov_zprod_regex = r"Mean z1\*z2: (.+)"
                    gcov_int_regex = r"Intercept: (.+)\s+\((.+)\)"

                    while re.match(gcov_regex, line) is None:
                        line = infile.readline()

                    gcov_match = re.match(gcov_regex, line)
                    gcov = float(gcov_match.group(1))
                    gcov_se = float(gcov_match.group(2))

                    line = infile.readline()

                    gcov_zprod_match = re.match(gcov_zprod_regex, line)
                    gcov_zprod = float(gcov_zprod_match.group(1))

                    line = infile.readline()

                    gcov_int_match = re.match(gcov_int_regex, line)
                    gcov_int = float(gcov_int_match.group(1))
                    gcov_int_se = float(gcov_int_match.group(2))

                    line = infile.readline()

                    while re.match("^p1\s", line) is None:
                        line = infile.readline()

                    line = infile.readline()

                    rg, rg_se, rg_z, rg_p = [float(z) if z != 'NA' else nan for z in line.split()[2:6]]

                    outfile.write(f"{ncases_A}\t{ncontrols_A}\t{ncases_B}\t{ncontrols_B}\t{odds_ratio_A}\t{odds_ratio_B}\t{effect_blocks_A}\t{effect_blocks_B}\t{shared_effect_blocks}\t{tag_pair}\t{h2_A:.4}\t{h2_A_se:.4}\t{h2_B:.4}\t{h2_B_se:.4}\t{gcov:.4}\t{gcov_se:.4}\t{gcov_int:.4}\t{gcov_int_se:.4}\t{gcov_zprod:.4}\t{rg:.4}\t{rg_se:.4}\t{rg_p:.4}\n")

rule run_medium_effect_gps_hoeffdings:
    input:
        [x[0] % x[1] for x in zip(medium_effect_gps_files, range(100, 100+len(medium_effect_gps_files)))]+
        [x[0] % x[1] for x in zip(medium_effect_hoeffdings_files, range(100, 100+len(medium_effect_hoeffdings_files)))]
#        [x[0] % x[1] for x in zip(medium_effect_rg_estimate_files, range(100, 100+len(medium_effect_rg_estimate_files)))]+

"""
# TODO need to remove permute output for incorrect files as well
# TODO GPS
# TODO Hoeffdings
# TODO Check for temporary sum stats files

rule remove_erroneous_outputs:
    input:
     input_files = [x.strip() for x in open('files_to_delete.txt', 'r').readlines()]
    output:
     'files_deleted.done'
    run:
        for x in input.input_files:
            shell("rm {x}")

        shell("touch {output}")
"""
