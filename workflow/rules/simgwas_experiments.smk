import re
import os
from numpy import nan

tag_pairs = [tags[i]+tags[i+1] for i in range(0, 19, 2)]

odds_ratio_dict = {"s": 1.05, "m": 1.2, 'l': 1.4, 'v': 2, 'r': 'random', 'n' : 1, 'i' : 1.1}

medium_effect_tuples = [f"m50_m50_m{x}" for x in [0, 10, 20, 30, 40, 50]]+[f"m25_m25_m{x}" for x in [0, 5, 10, 15, 20, 25]]

sample_sizes = [1000, 5000, 10000, 50000, 100000, 250000]

rule compile_medium_rg_estimates:
    input:
        [x[0] % x[1] for x in zip([f"results/ldsc/rg/whole_genome/randomised/{size}_{size}_{size}_{size}/{effect_tuple}_seed_%s_{tag_pair}.log" for size in sample_sizes for tag_pair in tag_pairs for effect_tuple in medium_effect_tuples], range(100, 1541))][0]

#rule compile_small_h2_theo_calculations:
#    input:
#        ["results/ldsc/h2/whole_genome/1-s0:%d_theo_h2.tsv" % j for j in [9, 19, 39, 79, 99, 119]],
#        ["results/ldsc/h2/whole_genome/1-s0:119+2-s0:%d_theo_h2.tsv" % j for j in [9, 19, 39, 79, 99, 139]],
#        ["results/ldsc/h2/whole_genome/1-s0:119+2-s0:139+3-s0:%d_theo_h2.tsv" % j for j in [9, 19, 39, 79, 99, 119]],
#        ["results/ldsc/h2/whole_genome/1-s0:119+2-s0:139+3-s0:119+4-s0:%d_theo_h2.tsv" % j for j in [9, 19, 39, 79, 99, 118]],
#        ["results/ldsc/h2/whole_genome/1-s0:119+2-s0:139+3-s0:119+4-s0:118+5-s0:%d_theo_h2.tsv" % j for j in [9, 19, 39, 79, 99, 107]],
#        ["results/ldsc/h2/whole_genome/1-s0:119+2-s0:139+3-s0:119+4-s0:118+5-s0:107+6-s0:%d_theo_h2.tsv" % j for j in [9, 19, 39, 79, 99, 105]],
#        ["results/ldsc/h2/whole_genome/1-s0:119+2-s0:139+3-s0:119+4-s0:118+5-s0:107+6-s0:105+7-s0:%d_theo_h2.tsv" % j for j in [9, 19, 39, 79, 91]],
#    output:
#        "results/ldsc/h2/whole_genome/compiled_s_theo_h2.tsv"
#    run:
#        for i,x in enumerate(input):
#            if i == 0:
#                shell("cat %s > %s" % (x, output[0]))
#            else:
#                shell("cat %s | tail -n +2  >> %s" % (x, output[0]))

#rule compile_medium_h2_theo_calculations:
#    input:
#        ["results/ldsc/h2/whole_genome/1-m0:%d_theo_h2.tsv" % j for j in [9, 19, 39, 79, 99, 119]],
#        ["results/ldsc/h2/whole_genome/1-m0:119+2-m0:%d_theo_h2.tsv" % j for j in [9, 19, 39, 79, 99, 139]],
#        ["results/ldsc/h2/whole_genome/1-m0:119+2-m0:139+3-m0:%d_theo_h2.tsv" % j for j in [9, 19, 39, 79, 99, 119]],
#        ["results/ldsc/h2/whole_genome/1-m0:119+2-m0:139+3-m0:119+4-m0:%d_theo_h2.tsv" % j for j in [9, 19, 39, 79, 99, 118]]
#    output:
#        "results/ldsc/h2/whole_genome/compiled_m_theo_h2.tsv"
#    run:
#        for i,x in enumerate(input):
#            if i == 0:
#                shell("cat %s > %s" % (x, output[0]))
#            else:
#                shell("cat %s | tail -n +2  >> %s" % (x, output[0]))
#
#rule compile_small_rg_theo_calculations:
#    input:
#        ["results/ldsc/rg/whole_genome/%s_theo_rg.tsv" % x for x in small_effect_block_pairs]
#    output:
#        "results/ldsc/rg/whole_genome/compiled_s_theo_rg.tsv"
#    run:
#        for i,x in enumerate(input):
#            if i == 0:
#                shell("cat %s > %s" % (x, output[0]))
#            else:
#                shell("cat %s | tail -n +2  >> %s" % (x, output[0]))
#
#rule compile_medium_theo_rg_calculations:
#    input:
#        ["results/ldsc/rg/whole_genome/%s_theo_rg.tsv" % x for x in medium_effect_tuples]
#    output:
#        "results/ldsc/rg/whole_genome/compiled_m_theo_rg.tsv"
#    run:
#        for i,x in enumerate(input):
#            if i == 0:
#                shell("cat %s > %s" % (x, output[0]))
#            else:
#                shell("cat %s | tail -n +2  >> %s" % (x, output[0]))

"""
rule compile_small_rg_estimates:
    input:
        [x for x in ["results/ldsc/rg/whole_genome/%d_%d_%d_%d/%s_%s.log" % (i,i,i,i,e,s) for i in sample_sizes for s in tag_pairs for e in
         small_effect_block_pairs] if x not in missing_rg_files]
    output:
        "results/ldsc/rg/whole_genome/compiled_s_rg_estimates.tsv"
    run:
        with open(output[0], 'w') as outfile:
            outfile.write("ncases.A\tncontrols.A\tncases.B\tncontrols.B\todds_ratio.A\todds_ratio.B\tblocks.A\tblocks.B\tno_shared_blocks\ttag_pair\th2.A\th2.A.se\th2.B\th2.B.se\trg\trg.se\trg.p\n")
            for x in input:
                head, tail = os.path.split(x)
                head_res = re.match("results/ldsc/rg/whole_genome/(\d+)_(\d+)_(\d+)_(\d+)", head)
                ncases_A = int(head_res.group(1))
                ncontrols_A = int(head_res.group(2))
                ncases_B = int(head_res.group(3))
                ncontrols_B = int(head_res.group(4))

                effect_blocks_wc_A, effect_blocks_wc_B, tag_pair = tail.split('_')

                tag_pair = tag_pair.split('.log')[0]

                odds_ratios_A = []
                odds_ratios_B = []

                effect_blocks_A = []
                effect_blocks_B = []

                for y in effect_blocks_wc_A.split('+'):
                    if effect_blocks_wc_A != 'null':
                        block_match = re.match('^(\d+)-(.+)', y)

                        chrom = int(block_match.group(1))

                        range_match = re.match('([smlvhi])(\d+):(\d+)', block_match.group(2))

                        odds_ratios_A.append(odds_ratio_dict[range_match.group(1)])

                        effect_blocks_A += ["%d-%s%d" % (chrom, range_match.group(1), z) for z in range(int(range_match.group(2)), int(range_match.group(3))+1) if z in block_dict[chrom]]
                    else:
                        odds_ratios_A.append(odds_ratio_dict['n'])

                for y in effect_blocks_wc_B.split('+'):
                    if effect_blocks_wc_B != 'null':
                        block_match = re.match('^(\d+)-(.+)', y)

                        chrom = int(block_match.group(1))

                        range_match = re.match('([smlvhi])(\d+):(\d+)', block_match.group(2))

                        odds_ratios_B.append(odds_ratio_dict[range_match.group(1)])

                        effect_blocks_B += ["%d-%s%d" % (chrom, range_match.group(1), z) for z in range(int(range_match.group(2)), int(range_match.group(3))+1) if z in block_dict[chrom]]
                    else:
                        odds_ratios_B.append(odds_ratio_dict['n'])

                odds_ratios_A = ','.join(set([str(x) for x in odds_ratios_A]))
                odds_ratios_B = ','.join(set([str(x) for x in odds_ratios_B]))

                no_shared_blocks = len([z for z in effect_blocks_A if z in effect_blocks_B])

                h2_regex = r"Total Liability scale h2: (.+)\s+\((.+)\)"

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

                    while re.match("^p1\s", line) is None:
                        line = infile.readline()

                    line = infile.readline()

                    rg, rg_se, rg_z, rg_p = [float(z) if z != 'NA' else nan for z in line.split()[2:6]]
                    outfile.write(f"{ncases_A}\t{ncontrols_A}\t{ncases_B}\t{ncontrols_B}\t{odds_ratios_A}\t{odds_ratios_B}\t{effect_blocks_wc_A}\t{effect_blocks_wc_B}\t{no_shared_blocks}\t{tag_pair}\t{h2_A:.4}\t{h2_A_se:.4}\t{h2_B:.4}\t{h2_B_se:.4}\t{rg:.4}\t{rg_se:.4}\t{rg_p:.4}\n")
"""

#    output:
#        "results/ldsc/rg/whole_genome/compiled_m_rg_estimates.tsv"
#    run:
#        with open(output[0], 'w') as outfile:
#            outfile.write("ncases.A\tncontrols.A\tncases.B\tncontrols.B\todds_ratio.A\todds_ratio.B\tblocks.A\tblocks.B\tno_shared_blocks\ttag_pair\th2.A\th2.A.se\th2.B\th2.B.se\tgcov\tgcov.se\tgcov.int\tgcov.int.se\trg\trg.se\trg.p\n")
#            for x in input:
#                head, tail = os.path.split(x)
#                head_res = re.match("results/ldsc/rg/whole_genome/(\d+)_(\d+)_(\d+)_(\d+)", head)
#                ncases_A = int(head_res.group(1))
#                ncontrols_A = int(head_res.group(2))
#                ncases_B = int(head_res.group(3))
#                ncontrols_B = int(head_res.group(4))
#
#                effect_blocks_wc_A, effect_blocks_wc_B, tag_pair = tail.split('_')
#
#                tag_pair = tag_pair.split('.log')[0]
#
#                odds_ratios_A = []
#                odds_ratios_B = []
#
#                effect_blocks_A = []
#                effect_blocks_B = []
#
#                for y in effect_blocks_wc_A.split('+'):
#                    if effect_blocks_wc_A != 'null':
#                        block_match = re.match('^(\d+)-(.+)', y)
#
#                        chrom = int(block_match.group(1))
#
#                        range_match = re.match('([smlvhi])(\d+):(\d+)', block_match.group(2))
#
#                        odds_ratios_A.append(odds_ratio_dict[range_match.group(1)])
#
#                        effect_blocks_A += ["%d-%s%d" % (chrom, range_match.group(1), z) for z in range(int(range_match.group(2)), int(range_match.group(3))+1) if z in block_dict[chrom]]
#                    else:
#                        odds_ratios_A.append(odds_ratio_dict['n'])
#
#                for y in effect_blocks_wc_B.split('+'):
#                    if effect_blocks_wc_B != 'null':
#                        block_match = re.match('^(\d+)-(.+)', y)
#
#                        chrom = int(block_match.group(1))
#
#                        range_match = re.match('([smlvhi])(\d+):(\d+)', block_match.group(2))
#
#                        odds_ratios_B.append(odds_ratio_dict[range_match.group(1)])
#
#                        effect_blocks_B += ["%d-%s%d" % (chrom, range_match.group(1), z) for z in range(int(range_match.group(2)), int(range_match.group(3))+1) if z in block_dict[chrom]]
#                    else:
#                        odds_ratios_B.append(odds_ratio_dict['n'])
#
#                odds_ratios_A = ','.join(set([str(x) for x in odds_ratios_A]))
#                odds_ratios_B = ','.join(set([str(x) for x in odds_ratios_B]))
#
#                no_shared_blocks = len([z for z in effect_blocks_A if z in effect_blocks_B])
#
#                h2_regex = r"Total Liability scale h2: (.+)\s+\((.+)\)"
#
#                with open(x, 'r') as infile:
#                    line = infile.readline()
#
#                    # TODO fix these for the null case
#                    while re.match(h2_regex, line) is None:
#                        line = infile.readline()
#
#                    h2_match_A = re.match(h2_regex, line)
#                    h2_A = float(h2_match_A.group(1))
#                    h2_A_se = float(h2_match_A.group(2))
#
#                    line = infile.readline()
#
#                    while re.match(h2_regex, line) is None:
#                        line = infile.readline()
#
#                    h2_match_B = re.match(h2_regex, line)
#                    h2_B = float(h2_match_B.group(1))
#                    h2_B_se = float(h2_match_B.group(2))
#
#                    line = infile.readline()
#
#                    while re.match("^p1\s", line) is None:
#                        line = infile.readline()
#
#                    line = infile.readline()
#
#                    rg, rg_se, rg_z, rg_p = [float(z) if z != 'NA' else nan for z in line.split()[2:6]]
#
#                    outfile.write(f"{ncases_A}\t{ncontrols_A}\t{ncases_B}\t{ncontrols_B}\t{odds_ratios_A}\t{odds_ratios_B}\t{effect_blocks_wc_A}\t{effect_blocks_wc_B}\t{no_shared_blocks}\t{tag_pair}\t{h2_A:.4}\t{h2_A_se:.4}\t{h2_B:.4}\t{h2_B_se:.4}\t{rg:.4}\t{rg_se:.4}\t{rg_p:.4}\n")
