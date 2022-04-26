import re
import os
from numpy import nan

tag_pairs = [tags[i]+tags[i+1] for i in range(0, 19, 2)]

odds_ratio_dict = {"s": 1.05, "m": 1.2, 'l': 1.4, 'v': 2, 'r': 'random', 'n' : 1, 'i' : 1.1}

# 22
medium_effect_tuples = [f"m50_m50_m{x}" for x in [0, 10, 20, 30, 40, 50]]+[f"m25_m25_m{x}" for x in [0, 5, 10, 15, 20, 25]]+[f"m50_m50_m{x}" for x in [2, 5, 7, 12, 15, 17]]+[f"m25_m25_m{x}" for x in [2, 7, 12, 17]]

small_effect_tuples = [f"s400_s400_s{x}" for x in range(0, 101, 10)]

sample_sizes = [(1000, 1000, 1000, 1000),
                (5000, 5000, 5000, 5000),
                (10000, 10000, 10000, 10000),
                (50000, 50000, 50000, 50000),
                (100000, 100000, 100000, 100000),
                (250000, 250000, 250000, 250000),
                (500, 10000, 500, 10000),
                (1000, 10000, 1000, 10000),
                (2000, 10000, 2000, 10000),
                (5000, 10000, 5000, 10000)]

medium_effect_rg_estimate_files  = [f"results/ldsc/rg/whole_genome/randomised/{size[0]}_{size[1]}_{size[2]}_{size[3]}/{effect_tuple}_seed_%d_{tag_pair}.log" for size in sample_sizes for tag_pair in tag_pairs for effect_tuple in medium_effect_tuples]

medium_effect_rg_estimate_no_overlap_files  = [f"results/ldsc/rg/whole_genome/randomised/{size[0]}_{size[1]}_{size[2]}_{size[3]}/{effect_tuple}_seed_%d_{tag_pair}_intercept_0.log" for size in sample_sizes for tag_pair in tag_pairs for effect_tuple in medium_effect_tuples]

medium_effect_gps_files = [f"results/gps/simgwas/randomised/window_1000kb_step_50/{size[0]}_{size[1]}_{size[2]}_{size[3]}/3000_permutations/{effect_tuple}_seed_%d_tags_{tag_pair}_gps_pvalue.tsv" for size in sample_sizes for tag_pair in tag_pairs for effect_tuple in medium_effect_tuples]

medium_effect_hoeffdings_files = [f"results/hoeffdings/simgwas/randomised/window_1000kb_step_50/{size[0]}_{size[1]}_{size[2]}_{size[3]}/{effect_tuple}_seed_%d_tags_{tag_pair}_hoeffdings.tsv" for size in sample_sizes for tag_pair in tag_pairs for effect_tuple in medium_effect_tuples]

small_effect_rg_estimate_files = [f"results/ldsc/rg/whole_genome/randomised/{size[0]}_{size[1]}_{size[2]}_{size[3]}/{effect_tuple}_seed_%d_{tag_pair}.log" for size in sample_sizes for tag_pair in tag_pairs for effect_tuple in small_effect_tuples]

small_effect_gps_files = [f"results/gps/simgwas/randomised/window_1000kb_step_50/{size[0]}_{size[1]}_{size[2]}_{size[3]}/3000_permutations/{effect_tuple}_seed_%d_tags_{tag_pair}_gps_pvalue.tsv" for size in sample_sizes for tag_pair in tag_pairs for effect_tuple in small_effect_tuples]

small_effect_hoeffdings_files = [f"results/hoeffdings/simgwas/randomised/window_1000kb_step_50/{size[0]}_{size[1]}_{size[2]}_{size[3]}/{effect_tuple}_seed_%d_tags_{tag_pair}_hoeffdings.tsv" for size in sample_sizes for tag_pair in tag_pairs for effect_tuple in small_effect_tuples]

rule compile_medium_rg_estimates:
    input:
        [y for y in [x[0] % x[1] for x in zip(medium_effect_rg_estimate_files, range(100, 100+len(medium_effect_rg_estimate_files)))] if os.path.exists(y)]
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

rule run_medium_rg_gps_hoeffdings:
    input:
        [x[0] % x[1] for x in zip(medium_effect_rg_estimate_files, range(100, 100+len(medium_effect_rg_estimate_files)))]+
        [x[0] % x[1] for x in zip(medium_effect_gps_files, range(100, 100+len(medium_effect_gps_files)))]+
        [x[0] % x[1] for x in zip(medium_effect_hoeffdings_files, range(100, 100+len(medium_effect_hoeffdings_files)))]

rule compile_medium_gps_results:
    input:
        [x[0] % x[1] for x in zip(medium_effect_gps_files, range(100, 100+len(medium_effect_gps_files)))]
    output:
        "results/gps/simgwas/randomised/window_1000kb_step_50/compiled_m_gps_results.tsv"
    run:
        with open(output[0], 'w') as outfile:
            outfile.write("ncases.A\tncontrols.A\tncases.B\tncontrols.B\todds_ratio.A\todds_ratio.B\tblocks.A\tblocks.B\tno_shared_blocks\ttag_pair\tgps\tn\tgps.p\n")
            for x in input:
                head, tail = os.path.split(x)
                head_res = re.match("results/gps/simgwas/randomised/window_1000kb_step_50/(\d+)_(\d+)_(\d+)_(\d+)/3000_permutations", head)
                ncases_A = int(head_res.group(1))
                ncontrols_A = int(head_res.group(2))
                ncases_B = int(head_res.group(3))
                ncontrols_B = int(head_res.group(4))

                odds_ratio_A = 1.2
                odds_ratio_B = 1.2

                tail_res = re.match("m(\d+)_m(\d+)_m(\d+)_seed_(\d+)_tags_(\w{2})_gps_pvalue.tsv", tail)

                effect_blocks_A, effect_blocks_B, shared_effect_blocks, seed, tag_pair = tail_res.groups()

                with open(x, 'r') as infile:
                    lines = [x.strip() for x in infile.readlines()]
                    gps, n, loc, loc_sd, scale, scale_sd, shape, shape_sd, pval = lines[1].split('\t')
                    outfile.write(f"{ncases_A}\t{ncontrols_A}\t{ncases_B}\t{ncontrols_B}\t{odds_ratio_A}\t{odds_ratio_B}\t{effect_blocks_A}\t{effect_blocks_B}\t{shared_effect_blocks}\t{tag_pair}\t{gps}\t{n}\t{pval}\n")

rule compile_medium_effect_theoretical_rg:
    input:
        [x[0] % x[1] for x in zip([f"results/ldsc/rg/whole_genome/randomised/theoretical_rg/{size[0]}_{size[1]}_{size[2]}_{size[3]}/{effect_tuple}_seed_%d_{tag_pair}_theo_rg.tsv" for size in sample_sizes for tag_pair in tag_pairs for effect_tuple in medium_effect_tuples], range(100, 100+len(medium_effect_rg_estimate_files)))]
    output:
        compiled_rg_file = "results/ldsc/rg/whole_genome/randomised/theoretical_rg/compiled_m_theoretical_rg.tsv"
    run:
        with open(output.compiled_rg_file, 'w') as outfile:
            outfile.write("ncases.A\tncontrols.A\tncases.B\tncontrols.B\tseed\ttag_pair\todds_ratio.A\todds_ratio.B\tblocks.A\tblocks.B\tno_shared_blocks\th2.theo.obs.A\th2.theo.obs.B\th2.theo.liab.A\th2.theo.liab.B\tV_A.A\tV_A.B\tC_A.AB\tr_A.AB\n")

            for x in input:
                with open(x, 'r') as infile:
                    head, tail = os.path.split(x)

                    ncases_A, ncontrols_A, ncases_B, ncontrols_B = re.match("results/ldsc/rg/whole_genome/randomised/theoretical_rg/(\d+)_(\d+)_(\d+)_(\d+)", head).groups()

                    effect_blocks_A, effect_blocks_B, shared_effect_blocks, seed, tag_pair = re.match("m(\d+)_m(\d+)_m(\d+)_seed_(\d+)_(\w{2})_theo_rg.tsv", tail).groups()

                    lines = [y.strip() for y in infile.readlines()]

                    h2_theo_obs_A, h2_theo_obs_B, h2_theo_liab_A, h2_theo_liab_B, V_A_A, V_A_B, C_A_AB, r_A_AB = lines[1].split('\t')[5:]
                    outfile.write(f"{ncases_A}\t{ncontrols_A}\t{ncases_B}\t{ncontrols_B}\t{seed}\t{tag_pair}\t1.2\t1.2\t{effect_blocks_A}\t{effect_blocks_B}\t{shared_effect_blocks}\t{h2_theo_obs_A}\t{h2_theo_obs_B}\t{h2_theo_liab_A}\t{h2_theo_liab_B}\t{V_A_A}\t{V_A_B}\t{C_A_AB}\t{r_A_AB}\n")

rule compile_medium_hoeffdings_results:
    input:
        [y for y in [x[0] % x[1] for x in zip(medium_effect_hoeffdings_files, range(100, 100+len(medium_effect_hoeffdings_files)))] if os.path.exists(y)]
    output:
        "results/hoeffdings/simgwas/randomised/window_1000kb_step_50/compiled_m_hoeffdings_results.tsv"
    run:
        with open(output[0], 'w') as outfile:
            outfile.write("ncases.A\tncontrols.A\tncases.B\tncontrols.B\todds_ratio.A\todds_ratio.B\tblocks.A\tblocks.B\tno_shared_blocks\ttag_pair\thoeff.p\n")
            #        "results/gps/simgwas/randomised/window_1000kb_step_50/{size}_{size}_{size}_{size}/3000_permutations/{effect_tuple}_seed_%d_tags_{tag_pair}_gps_pvalue.tsv"
            for x in input:
                head, tail = os.path.split(x)
                head_res = re.match("results/hoeffdings/simgwas/randomised/window_1000kb_step_50/(\d+)_(\d+)_(\d+)_(\d+)", head)
                ncases_A = int(head_res.group(1))
                ncontrols_A = int(head_res.group(2))
                ncases_B = int(head_res.group(3))
                ncontrols_B = int(head_res.group(4))

                odds_ratio_A = 1.2
                odds_ratio_B = 1.2

                tail_res = re.match("m(\d+)_m(\d+)_m(\d+)_seed_(\d+)_tags_(\w{2})_hoeffdings.tsv", tail)

                effect_blocks_A, effect_blocks_B, shared_effect_blocks, seed, tag_pair = tail_res.groups()

                with open(x, 'r') as infile:
                    lines = [x.strip() for x in infile.readlines()]
                    # NB: The 'n' here is no. of SNPs, not no. of permutations as with GPS
                    _, _, n, Dn, scaled, pval = lines[1].split('\t')
                    outfile.write(f"{ncases_A}\t{ncontrols_A}\t{ncases_B}\t{ncontrols_B}\t{odds_ratio_A}\t{odds_ratio_B}\t{effect_blocks_A}\t{effect_blocks_B}\t{shared_effect_blocks}\t{tag_pair}\t{pval}\n")

rule compile_medium_stats:
    input:
        "results/ldsc/rg/whole_genome/randomised/compiled_m_rg_estimates.tsv",
        "results/gps/simgwas/randomised/window_1000kb_step_50/compiled_m_gps_results.tsv",
        "results/ldsc/rg/whole_genome/randomised/theoretical_rg/compiled_m_theoretical_rg.tsv",
        "results/hoeffdings/simgwas/randomised/window_1000kb_step_50/compiled_m_hoeffdings_results.tsv"

rule compile_mn_pair_block_numbers:
    input:
        [f"results/ldsc/rg/whole_genome/randomised/theoretical_rg/500_10000_500_10000/block_files/m{x[1]}_m{x[2]}_m{x[3]}_seed_{x[0]}_shared_mn.tsv" for x in [(1558,25,25,0),
            (1570,25,25,2),
            (1559,25,25,5),
            (1571,25,25,7),
            (1560,25,25,10),
            (1572,25,25,12),
            (1561,25,25,15),
            (1573,25,25,17),
            (1562,25,25,20),
            (1563,25,25,25),
            (1552,50,50,0),
            (1564,50,50,2),
            (1565,50,50,5),
            (1566,50,50,7),
            (1553,50,50,10),
            (1567,50,50,12),
            (1568,50,50,15),
            (1569,50,50,17),
            (1554,50,50,20),
            (1555,50,50,30),
            (1556,50,50,40),
            (1557,50,50,50)]]
