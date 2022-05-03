def compile_rg_results(input, output):
    with open(output[0], 'w') as outfile:
        outfile.write("ncases.A\tncontrols.A\tncases.B\tncontrols.B\todds_ratio.A\todds_ratio.B\th2.int.con\tgcov.int.con\tblocks.A\tblocks.B\tno_shared_blocks\ttag_pair\tseed\th2.A\th2.A.se\th2.int.A\th2.int.A.se\th2.B\th2.B.se\th2.int.B\th2.int.B.se\tgcov\tgcov.se\tgcov.int\tgcov.int.se\tgcov.zprod\trg\trg.se\trg.p\n")
        for x in input:
            head, tail = os.path.split(x)
            head_res = re.match("results/ldsc/rg/(\w+)_h2_(\w+)_rg_intercept/whole_genome/randomised/(\d+)_(\d+)_(\d+)_(\d+)", head)
            h2_intercept = 'fixed' if head_res.group(1) == 'fixed' else 'free'
            gcov_intercept = 'fixed' if head_res.group(2) == 'fixed' else 'free'
            ncases_A = int(head_res.group(3))
            ncontrols_A = int(head_res.group(4))
            ncases_B = int(head_res.group(5))
            ncontrols_B = int(head_res.group(6))
    
            tail_effect_res = re.match("(\w)\d+_(\w)\d+_\w\d+_seed_\d+_\w{2}", tail)
    
            effect_A, effect_B = tail_effect_res.groups()
    
            # NB: currently assuming effect size is same
            odds_ratio_A = odds_ratio_dict[effect_A]
            odds_ratio_B = odds_ratio_dict[effect_B]
    
            tail_res = re.match("\w(\d+)_\w(\d+)_\w(\d+)_seed_(\d+)_(\w{2})", tail)
    
            effect_blocks_A, effect_blocks_B, shared_effect_blocks, seed, tag_pair = tail_res.groups()
    
            h2_regex = r"Total Liability scale h2: (.+)\s+\((.+)\)"
            int_regex = r"Intercept: (.+)\s+\((.+)\)"
    
            gcov_regex = r"Total Liability scale gencov: (.+)\s+\((.+)\)"
            gcov_zprod_regex = r"Mean z1\*z2: (.+)"

            with open(x, 'r') as infile:
                line = infile.readline()

                # TODO fix these for the null case
                while re.match(h2_regex, line) is None:
                    line = infile.readline()

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

                outfile.write(f"{ncases_A}\t{ncontrols_A}\t{ncases_B}\t{ncontrols_B}\t{odds_ratio_A}\t{odds_ratio_B}\t{h2_intercept}\t{gcov_intercept}\t{effect_blocks_A}\t{effect_blocks_B}\t{shared_effect_blocks}\t{tag_pair}\t{seed}\t{h2_A:.4}\t{h2_A_se:.4}\t{h2_int_A:.4}\t{h2_int_A_se:.4}\t{h2_B:.4}\t{h2_B_se:.4}\t{h2_int_B:.4}\t{h2_int_B_se:.4}\t{gcov:.4}\t{gcov_se:.4}\t{gcov_int:.4}\t{gcov_int_se:.4}\t{gcov_zprod:.4}\t{rg:.4}\t{rg_se:.4}\t{rg_p:.4}\n")

    return

def compile_gps_results(input, output):
    with open(output[0], 'w') as outfile:
        outfile.write("ncases.A\tncontrols.A\tncases.B\tncontrols.B\todds_ratio.A\todds_ratio.B\tblocks.A\tblocks.B\tno_shared_blocks\ttag_pair\tseed\tgps\tn\tgps.p\n")
        for x in input:
            head, tail = os.path.split(x)
            head_res = re.match("results/gps/simgwas/randomised/window_1000kb_step_50/(\d+)_(\d+)_(\d+)_(\d+)/3000_permutations", head)
            ncases_A = int(head_res.group(1))
            ncontrols_A = int(head_res.group(2))
            ncases_B = int(head_res.group(3))
            ncontrols_B = int(head_res.group(4))

            tail_effect_res = re.match("(\w)\d+_(\w)\d+_\w\d+_seed_\d+_\w{2}", tail)

            effect_A, effect_B = tail_effect_res.groups()

            # NB: currently assuming effect size is same
            odds_ratio_A = odds_ratio_dict[effect_A]
            odds_ratio_B = odds_ratio_dict[effect_B]

            tail_res = re.match("\w(\d+)_\w(\d+)_\w(\d+)_seed_(\d+)_tags_(\w{2})_gps_pvalue.tsv", tail)

            effect_blocks_A, effect_blocks_B, shared_effect_blocks, seed, tag_pair = tail_res.groups()

            with open(x, 'r') as infile:
                lines = [x.strip() for x in infile.readlines()]
                gps, n, loc, loc_sd, scale, scale_sd, shape, shape_sd, pval = lines[1].split('\t')
                outfile.write(f"{ncases_A}\t{ncontrols_A}\t{ncases_B}\t{ncontrols_B}\t{odds_ratio_A}\t{odds_ratio_B}\t{effect_blocks_A}\t{effect_blocks_B}\t{shared_effect_blocks}\t{tag_pair}\t{seed}\t{gps}\t{n}\t{pval}\n")

    return

def compile_theoretical_rg_results(input, output):
    with open(output.compiled_rg_file, 'w') as outfile:
        outfile.write("ncases.A\tncontrols.A\tncases.B\tncontrols.B\tseed\ttag_pair\todds_ratio.A\todds_ratio.B\tblocks.A\tblocks.B\tno_shared_blocks\th2.theo.obs.A\th2.theo.obs.B\th2.theo.liab.A\th2.theo.liab.B\tV_A.A\tV_A.B\tC_A.AB\tr_A.AB\n")

        for x in input:
            with open(x, 'r') as infile:
                head, tail = os.path.split(x)

                ncases_A, ncontrols_A, ncases_B, ncontrols_B = re.match("results/ldsc/rg/whole_genome/randomised/theoretical_rg/(\d+)_(\d+)_(\d+)_(\d+)", head).groups()

                tail_effect_res = re.match("(\w)\d+_(\w)\d+_\w\d+_seed_\d+_\w{2}", tail)
                effect_A, effect_B = tail_effect_res.groups()

                # NB: currently assuming effect size is same
                odds_ratio_A = odds_ratio_dict[effect_A]
                odds_ratio_B = odds_ratio_dict[effect_B]

                effect_blocks_A, effect_blocks_B, shared_effect_blocks, seed, tag_pair = re.match("\w(\d+)_\w(\d+)_\w(\d+)_seed_(\d+)_(\w{2})_theo_rg.tsv", tail).groups()

                lines = [y.strip() for y in infile.readlines()]

                h2_theo_obs_A, h2_theo_obs_B, h2_theo_liab_A, h2_theo_liab_B, V_A_A, V_A_B, C_A_AB, r_A_AB = lines[1].split('\t')[5:]
                outfile.write(f"{ncases_A}\t{ncontrols_A}\t{ncases_B}\t{ncontrols_B}\t{seed}\t{tag_pair}\t{odds_ratio_A}\t{odds_ratio_B}\t{effect_blocks_A}\t{effect_blocks_B}\t{shared_effect_blocks}\t{h2_theo_obs_A}\t{h2_theo_obs_B}\t{h2_theo_liab_A}\t{h2_theo_liab_B}\t{V_A_A}\t{V_A_B}\t{C_A_AB}\t{r_A_AB}\n")

    return

def compile_hoeffdings_results(input, output):
    with open(output[0], 'w') as outfile:
        outfile.write("ncases.A\tncontrols.A\tncases.B\tncontrols.B\todds_ratio.A\todds_ratio.B\tblocks.A\tblocks.B\tno_shared_blocks\ttag_pair\tseed\thoeff.p\n")
        for x in input:
            head, tail = os.path.split(x)
            head_res = re.match("results/hoeffdings/simgwas/randomised/window_1000kb_step_50/(\d+)_(\d+)_(\d+)_(\d+)", head)
            ncases_A = int(head_res.group(1))
            ncontrols_A = int(head_res.group(2))
            ncases_B = int(head_res.group(3))
            ncontrols_B = int(head_res.group(4))

            tail_effect_res = re.match("(\w)\d+_(\w)\d+_\w\d+_seed_\d+_\w{2}", tail)

            effect_A, effect_B = tail_effect_res.groups()

            # NB: currently assuming effect size is same
            odds_ratio_A = odds_ratio_dict[effect_A]
            odds_ratio_B = odds_ratio_dict[effect_B]

            tail_res = re.match("\w(\d+)_\w(\d+)_\w(\d+)_seed_(\d+)_tags_(\w{2})_hoeffdings.tsv", tail)

            effect_blocks_A, effect_blocks_B, shared_effect_blocks, seed, tag_pair = tail_res.groups()

            with open(x, 'r') as infile:
                lines = [x.strip() for x in infile.readlines()]
                # NB: The 'n' here is no. of SNPs, not no. of permutations as with GPS
                _, _, n, Dn, scaled, pval = lines[1].split('\t')
                outfile.write(f"{ncases_A}\t{ncontrols_A}\t{ncases_B}\t{ncontrols_B}\t{odds_ratio_A}\t{odds_ratio_B}\t{effect_blocks_A}\t{effect_blocks_B}\t{shared_effect_blocks}\t{tag_pair}\t{seed}\t{pval}\n")

    return

def compile_sumher_results(input, output):
    with open(output[0], 'w') as outfile:

        outfile.write("ncases.A\tncontrols.A\tncases.B\tncontrols.B\todds_ratio.A\todds_ratio.B\tblocks.A\tblocks.B\tno_shared_blocks\ttag_pair\tseed\th2.A\th2.A.se\th2.B\th2.B.se\tgcov\tgcov.se\trg\trg.se\n")
        for x in input:
            head, tail = os.path.split(x)
            head_res = re.match("results/simgwas/ldak/ldak-thin/rg/(\d+)_(\d+)_(\d+)_(\d+)", head)
            ncases_A = int(head_res.group(1))
            ncontrols_A = int(head_res.group(2))
            ncases_B = int(head_res.group(3))
            ncontrols_B = int(head_res.group(4))

            tail_effect_res = re.match("(\w)\d+_(\w)\d+_\w\d+_seed_\d+_\w{2}", tail)

            effect_A, effect_B = tail_effect_res.groups()

            # NB: currently assuming effect size is same
            odds_ratio_A = odds_ratio_dict[effect_A]
            odds_ratio_B = odds_ratio_dict[effect_B]

            tail_res = re.match("\w(\d+)_\w(\d+)_\w(\d+)_seed_(\d+)_(\w{2}).cors.full", tail)

            effect_blocks_A, effect_blocks_B, shared_effect_blocks, seed, tag_pair = tail_res.groups()

            with open(x, 'r') as infile:
                line = infile.readline()
                line = infile.readline()

                # Category Trait1_Her SD Trait2_Her SD Both_Coher SD Correlation SD
                _, h2_A, h2_A_se, h2_B, h2_B_se, cov, cov_se, rg, rg_se = line.split()
                outfile.write(f"{ncases_A}\t{ncontrols_A}\t{ncases_B}\t{ncontrols_B}\t{odds_ratio_A}\t{odds_ratio_B}\t{effect_blocks_A}\t{effect_blocks_B}\t{shared_effect_blocks}\t{tag_pair}\t{seed}\t{float(h2_A)}\t{float(h2_A_se)}\t{float(h2_B)}\t{float(h2_B_se)}\t{float(cov)}\t{float(cov_se)}\t{float(rg)}\t{float(rg_se)}\n")

    return
