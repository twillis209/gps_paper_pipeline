import re
import os
from numpy import nan

odds_ratio_dict = {"s": 1.05, "m": 1.2, 'l': 1.4, 'v': 2, 'r': 'random', 'n' : 1, 'i' : 1.1}

rule whole_genome_B0_1_h2_chr1_7_sweep:
    input:
        [[[["results/ldsc/h2/whole_genome/%d_%d/B0_1/1-{effect_size}0:%d_%s.log" % (i, i, j, k) for i in [10000, 50000, 100000, 250000]] for j in [9, 19, 39, 79, 99, 119]]] for k in ["a", "b", "c", "d", "e"]],
        [[[["results/ldsc/h2/whole_genome/%d_%d/B0_1/1-{effect_size}0:119+2-{effect_size}0:%d_%s.log" % (i, i, j, k) for i in [10000, 50000, 100000, 250000]] for j in [9, 19, 39, 79, 99, 139]]] for k in ["a", "b", "c", "d", "e"]],
        [[[["results/ldsc/h2/whole_genome/%d_%d/B0_1/1-{effect_size}0:119+2-{effect_size}0:139+3-{effect_size}0:%d_%s.log" % (i, i, j, k) for i in [10000, 50000, 100000, 250000]] for j in [9, 19, 39, 79, 99, 119]]] for k in ["a", "b", "c", "d", "e"]],
        [[[["results/ldsc/h2/whole_genome/%d_%d/B0_1/1-{effect_size}0:119+2-{effect_size}0:139+3-{effect_size}0:119+4-{effect_size}0:%d_%s.log" % (i, i, j, k) for i in [10000, 50000, 100000, 250000]] for j in [9, 19, 39, 79, 99, 118]]] for k in ["a", "b", "c", "d", "e"]],
        [[[["results/ldsc/h2/whole_genome/%d_%d/B0_1/1-{effect_size}0:119+2-{effect_size}0:139+3-{effect_size}0:119+4-{effect_size}0:118+5-{effect_size}0:%d_%s.log" % (i, i, j, k) for i in [10000, 50000, 100000, 250000]] for j in [9, 19, 39, 79, 99, 107]]] for k in ["a", "b", "c", "d", "e"]],
        [[[["results/ldsc/h2/whole_genome/%d_%d/B0_1/1-{effect_size}0:119+2-{effect_size}0:139+3-{effect_size}0:119+4-{effect_size}0:118+5-{effect_size}0:107+6-{effect_size}0:%d_%s.log" % (i, i, j, k) for i in [10000, 50000, 100000, 250000]] for j in [9, 19, 39, 79, 99, 105]]] for k in ["a", "b", "c", "d", "e"]],
        [[[["results/ldsc/h2/whole_genome/%d_%d/B0_1/1-{effect_size}0:119+2-{effect_size}0:139+3-{effect_size}0:119+4-{effect_size}0:118+5-{effect_size}0:107+6-{effect_size}0:105+7-{effect_size}0:%d_%s.log" % (i, i, j, k) for i in [10000, 50000, 100000, 250000]] for j in [9, 19, 39, 79, 91]]] for k in ["a", "b", "c", "d", "e"]]
    output:
        "results/ldsc/h2/compiled/whole_genome_B0_1_chr1-7_{effect_size}_sweep.tsv"
    run:
        with open(output[0], 'w') as outfile:
            outfile.write("ncases\tncontrols\todds_ratio\tno_blocks\tchroms\ttag\th2\tse\n")
            for x in input:
                head, tail = os.path.split(x)
                head_res = re.match("results/ldsc/h2/whole_genome/(\d+)_(\d+)/B0_1", head)
                ncases = int(head_res.group(1))
                ncontrols = int(head_res.group(2))
                # TODO handle null
                no_blocks = sum([int(y.split(':')[-1])+1 for y in re.findall('[smlvhr]0:\d+', tail)])
                tag = (tail.split('_')[-1]).split('.')[0]
                chroms = ','.join([y[:-1] for y in re.findall('\d+-', tail)])
                odds_ratio = odds_ratio_dict[wildcards.effect_size]

                with open(x, 'r') as infile:
                    for line in infile:
                        h2_regex = r"Total Liability scale h2: (-?\d+\.\d+)\s+\((\d+\.\d+)\)"
                        if re.search(h2_regex, line):
                            res = re.search(h2_regex, line)
                            h2 = float(res.group(1))
                            se = float(res.group(2))
                            break
                outfile.write("%d\t%d\t%.2f\t%d\t%s\t%s\t%.4f\t%.4f\n" % (ncases, ncontrols, odds_ratio, no_blocks, chroms, tag, h2, se))

rule compile_small_rg_theo_calculations:
    input:
        "results/ldsc/rg/whole_genome/1-s0:119+2-s0:139+3-s0:119_5-s0:108+6-s0:105+7-s0:91+8-s0:72_theo_rg.tsv",
        "results/ldsc/rg/whole_genome/1-s0:119+2-s0:139+3-s0:119+4-s0:9_4-s0:9+5-s0:108+6-s0:105+7-s0:91+8-s0:72_theo_rg.tsv",
        "results/ldsc/rg/whole_genome/1-s0:119+2-s0:139+3-s0:119+4-s0:19_4-s0:19+5-s0:108+6-s0:105+7-s0:91+8-s0:72_theo_rg.tsv",
        "results/ldsc/rg/whole_genome/1-s0:119+2-s0:139+3-s0:119+4-s0:29_4-s0:29+5-s0:108+6-s0:105+7-s0:91+8-s0:72_theo_rg.tsv",
        "results/ldsc/rg/whole_genome/1-s0:119+2-s0:139+3-s0:119+4-s0:39_4-s0:39+5-s0:108+6-s0:105+7-s0:91+8-s0:72_theo_rg.tsv",
        "results/ldsc/rg/whole_genome/1-s0:119+2-s0:139+3-s0:119+4-s0:49_4-s0:49+5-s0:108+6-s0:105+7-s0:91+8-s0:72_theo_rg.tsv"
    output:
        "results/ldsc/rg/whole_genome/compiled_s_theo_rg.tsv"
    run:
        for i,x in enumerate(input):
            if i == 0:
                shell("cat %s > %s" % (x, output[0]))
            else:
                shell("cat %s | tail -n +2  >> %s" % (x, output[0]))

rule compile_medium_theo_rg_calculations:
    input:
        "results/ldsc/rg/whole_genome/1-m0:49_1-m50:99_theo_rg.tsv",
        "results/ldsc/rg/whole_genome/1-m0:49_1-m40:89_theo_rg.tsv",
        "results/ldsc/rg/whole_genome/1-m0:49_1-m30:79_theo_rg.tsv",
        "results/ldsc/rg/whole_genome/1-m0:49_1-m20:69_theo_rg.tsv",
        "results/ldsc/rg/whole_genome/1-m0:49_1-m10:59_theo_rg.tsv",
        "results/ldsc/rg/whole_genome/1-m0:49_1-m0:49_theo_rg.tsv"
    output:
        "results/ldsc/rg/whole_genome/compiled_m_theo_rg.tsv"
    run:
        for i,x in enumerate(input):
            if i == 0:
                shell("cat %s > %s" % (x, output[0]))
            else:
                shell("cat %s | tail -n +2  >> %s" % (x, output[0]))

# TODO remove: 380 s variants are enough on h2 front
# 5-s0:108+6-s0:105+7-s0:91+8-s0:72
rule compile_small_rg_estimates:
    input:
        ["results/ldsc/rg/whole_genome/%d_%d_%d_%d/1-s0:119+2-s0:139+3-s0:119_5-s0:108+6-s0:105+7-s0:91+8-s0:72.log" % (i,i,i,i) for i in [10000, 50000, 100000, 250000]],
        ["results/ldsc/rg/whole_genome/%d_%d_%d_%d/1-s0:119+2-s0:139+3-s0:119+4-s0:9_4-s0:9+5-s0:108+6-s0:105+7-s0:91+8-s0:72.log" % (i,i,i,i) for i in [10000, 50000, 100000, 250000]],
        ["results/ldsc/rg/whole_genome/%d_%d_%d_%d/1-s0:119+2-s0:139+3-s0:119+4-s0:19_4-s0:19+5-s0:108+6-s0:105+7-s0:91+8-s0:72.log" % (i,i,i,i) for i in [10000, 50000, 100000, 250000]],
        # 10k of the following doesn't work
        ["results/ldsc/rg/whole_genome/%d_%d_%d_%d/1-s0:119+2-s0:139+3-s0:119+4-s0:29_4-s0:29+5-s0:108+6-s0:105+7-s0:91+8-s0:72.log" % (i,i,i,i) for i in [50000, 100000, 250000]],
        # 10k of the following doesn't work
        ["results/ldsc/rg/whole_genome/%d_%d_%d_%d/1-s0:119+2-s0:139+3-s0:119+4-s0:39_4-s0:39+5-s0:108+6-s0:105+7-s0:91+8-s0:72.log" % (i,i,i,i) for i in [50000, 100000, 250000]],
        ["results/ldsc/rg/whole_genome/%d_%d_%d_%d/1-s0:119+2-s0:139+3-s0:119+4-s0:49_4-s0:49+5-s0:108+6-s0:105+7-s0:91+8-s0:72.log" % (i,i,i,i) for i in [10000, 50000, 100000, 250000]]
    output:
        "results/ldsc/rg/whole_genome/compiled_s_rg_estimates.tsv"
    run:
        with open(output[0], 'w') as outfile:
            outfile.write("ncases.A\tncontrols.A\tncases.B\tncontrols.B\todds_ratio.A\todds_ratio.B\tblocks.A\tblocks.B\tno_shared_blocks\th2.A\th2.A.se\th2.B\th2.B.se\trg\trg.se\n")
            for x in input:
                head, tail = os.path.split(x)
                head_res = re.match("results/ldsc/rg/whole_genome/(\d+)_(\d+)_(\d+)_(\d+)", head)
                ncases_A = int(head_res.group(1))
                ncontrols_A = int(head_res.group(2))
                ncases_B = int(head_res.group(3))
                ncontrols_B = int(head_res.group(4))

                effect_blocks_wc_A, effect_blocks_wc_B = tail.split('.log')[0].split('_')

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

                h2_regex = r"Total Liability scale h2: (-?\d+\.\d+)\s+\((\d+\.\d+)\)"

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

                    outfile.write("%d\t%d\t%d\t%d\t%s\t%s\t%s\t%s\t%d\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n" % (ncases_A, ncontrols_A, ncases_B, ncontrols_B, odds_ratios_A, odds_ratios_B, effect_blocks_wc_A, effect_blocks_wc_B, no_shared_blocks, h2_A, h2_A_se, h2_B, h2_B_se, rg, rg_se))

rule compile_rg_estimates:
    input:
#        ["results/ldsc/rg/whole_genome/%d_%d_%d_%d/null_null.log" % (i,i,i,i) for i in [10000, 50000, 100000, 250000]],
        ["results/ldsc/rg/whole_genome/%d_%d_%d_%d/1-m0:49_1-m50:99.log" % (i,i,i,i) for i in [10000, 50000, 100000, 250000]],
        ["results/ldsc/rg/whole_genome/%d_%d_%d_%d/1-m0:49_1-m40:89.log" % (i,i,i,i) for i in [10000, 50000, 100000, 250000]],
        ["results/ldsc/rg/whole_genome/%d_%d_%d_%d/1-m0:49_1-m30:79.log" % (i,i,i,i) for i in [10000, 50000, 100000, 250000]],
        ["results/ldsc/rg/whole_genome/%d_%d_%d_%d/1-m0:49_1-m20:69.log" % (i,i,i,i) for i in [10000, 50000, 100000, 250000]],
        ["results/ldsc/rg/whole_genome/%d_%d_%d_%d/1-m0:49_1-m10:59.log" % (i,i,i,i) for i in [10000, 50000, 100000, 250000]],
        ["results/ldsc/rg/whole_genome/%d_%d_%d_%d/1-m0:49_1-m0:49.log" % (i,i,i,i) for i in [10000, 50000, 100000, 250000]]
    output:
        "results/ldsc/rg/whole_genome/compiled_m_rg_estimates.tsv"
    run:
        with open(output[0], 'w') as outfile:
            outfile.write("ncases.A\tncontrols.A\tncases.B\tncontrols.B\todds_ratio.A\todds_ratio.B\tblocks.A\tblocks.B\tno_shared_blocks\th2.A\th2.A.se\th2.B\th2.B.se\trg\trg.se\n")
            for x in input:
                head, tail = os.path.split(x)
                head_res = re.match("results/ldsc/rg/whole_genome/(\d+)_(\d+)_(\d+)_(\d+)", head)
                ncases_A = int(head_res.group(1))
                ncontrols_A = int(head_res.group(2))
                ncases_B = int(head_res.group(3))
                ncontrols_B = int(head_res.group(4))

                effect_blocks_wc_A, effect_blocks_wc_B = tail.split('.log')[0].split('_')

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

                h2_regex = r"Total Liability scale h2: (-?\d+\.\d+)\s+\((\d+\.\d+)\)"

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

                    outfile.write("%d\t%d\t%d\t%d\t%s\t%s\t%s\t%s\t%d\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n" % (ncases_A, ncontrols_A, ncases_B, ncontrols_B, odds_ratios_A, odds_ratios_B, effect_blocks_wc_A, effect_blocks_wc_B, no_shared_blocks, h2_A, h2_A_se, h2_B, h2_B_se, rg, rg_se))
