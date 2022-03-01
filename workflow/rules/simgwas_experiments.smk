import re as re

odds_ratio_dict = {"s": 1.05, "m": 1.2, 'l': 1.4, 'v': 2, 'r': 'random'}

rule compile_swept_chr1_h2_parameters:
    input:
        [[[["results/ldsc/h2/chr1/%d_%d/%s0:%d_%s.log" % (i, i, z, j, k) for i in [10000, 50000, 100000, 250000]] for j in [9, 19, 39, 79, 99, 119]] for z in ["s", "m", "l", "v"]] for k in ["a", "b", "c", "d", "e"]]
    output:
        "results/ldsc/h2/compiled/chr1_sweep.tsv"
    run:
        with open(output[0], 'w') as outfile:
            outfile.write("ncases\tncontrols\todds_ratio\tblocks\ttag\th2\tse\n")
            for x in input:
                res = re.match("results/ldsc/h2/chr1/([0-9]+)_([0-9]+)/([smlv])0:([0-9]+)_([abcde]).log", x)
                ncases = int(res.group(1))
                ncontrols = int(res.group(2))
                print(res.group(3))
                odds_ratio = odds_ratio_dict[res.group(3)]
                last_block = int(res.group(4))
                tag = res.group(5)
                with open(x, 'r') as infile:
                    for line in infile:
                        h2_regex = r"Total Liability scale h2: (-?\d+\.\d+)\s+\((\d+\.\d+)\)"
                        if re.search(h2_regex, line):
                            res = re.search(h2_regex, line)
                            h2 = float(res.group(1))
                            se = float(res.group(2))
                            break
                outfile.write("%d\t%d\t%.2f\t0-%d\t%s\t%.4f\t%.4f\n" % (ncases, ncontrols, odds_ratio, last_block, tag, h2, se))

rule compile_swept_chr1_B0_1_h2_parameters:
    input:
        [[[["results/ldsc/h2/chr1/%d_%d/B0_1/%s0:%d_%s.log" % (i, i, z, j, k) for i in [10000, 50000, 100000, 250000]] for j in [9, 19, 39, 79, 99, 119]] for z in ["s", "m", "l", "v"]] for k in ["a", "b", "c", "d", "e"]]
    output:
        "results/ldsc/h2/compiled/chr1_B0_1_sweep.tsv"
    run:
        with open(output[0], 'w') as outfile:
            outfile.write("ncases\tncontrols\todds_ratio\tblocks\ttag\th2\tse\n")
            for x in input:
                res = re.match("results/ldsc/h2/chr1/([0-9]+)_([0-9]+)/B0_1/([smlv])0:([0-9]+)_([abcde]).log", x)
                ncases = int(res.group(1))
                ncontrols = int(res.group(2))
                print(res.group(3))
                odds_ratio = odds_ratio_dict[res.group(3)]
                last_block = int(res.group(4))
                tag = res.group(5)
                with open(x, 'r') as infile:
                    for line in infile:
                        h2_regex = r"Total Liability scale h2: (-?\d+\.\d+)\s+\((\d+\.\d+)\)"
                        if re.search(h2_regex, line):
                            res = re.search(h2_regex, line)
                            h2 = float(res.group(1))
                            se = float(res.group(2))
                            break
                outfile.write("%d\t%d\t%.2f\t0-%d\t%s\t%.4f\t%.4f\n" % (ncases, ncontrols, odds_ratio, last_block, tag, h2, se))

rule compile_swept_chr1_random_effect_h2_parameters:
    input:
        [[[["results/ldsc/h2/chr1/%d_%d/%s0:%d_%s.log" % (i, i, z, j, k) for i in [10000, 50000, 100000, 250000]] for j in range(1,10)] for z in ["r"]] for k in ["a"]]
    output:
        "results/ldsc/h2/compiled/chr1_random_effect_sweep.tsv"
    run:
        with open(output[0], 'w') as outfile:
            outfile.write("ncases\tncontrols\todds_ratio\tblocks\ttag\th2\tse\n")
            for x in input:
                filename_res = re.match("results/ldsc/h2/chr1/([0-9]+)_([0-9]+)/([smlvr])0:([0-9]+)_([abcde]).log", x)
                ncases = int(filename_res.group(1))
                ncontrols = int(filename_res.group(2))
                odds_ratio = odds_ratio_dict[filename_res.group(3)]
                last_block = int(filename_res.group(4))
                tag = filename_res.group(5)
                with open(x, 'r') as infile:
                    for line in infile:
                        h2_regex = r"Total Liability scale h2:\s+(-?\d+\.\d+)\s+\((\d+\.\d+)\)"
                        if re.search(h2_regex, line):
                            h2_res = re.search(h2_regex, line)
                            h2 = float(h2_res.group(1))
                            se = float(h2_res.group(2))
                            break
                outfile.write("%d\t%d\t%s\t0-%d\t%s\t%.4f\t%.4f\n" % (ncases, ncontrols, odds_ratio, last_block, tag, h2, se))

rule compile_swept_chr1_random_effect_B0_1_h2_parameters:
    input:
        [[[["results/ldsc/h2/chr1/%d_%d/B0_1/%s0:%d_%s.log" % (i, i, z, j, k) for i in [10000, 50000, 100000, 250000]] for j in range(1,10)] for z in ["r"]] for k in ["a"]]
    output:
        "results/ldsc/h2/compiled/chr1_random_effect_B0_1_sweep.tsv"
    run:
        with open(output[0], 'w') as outfile:
            outfile.write("ncases\tncontrols\todds_ratio\tblocks\ttag\th2\tse\n")
            for x in input:
                filename_res = re.match("results/ldsc/h2/chr1/([0-9]+)_([0-9]+)/B0_1/([smlvr])0:([0-9]+)_([abcde]).log", x)
                ncases = int(filename_res.group(1))
                ncontrols = int(filename_res.group(2))
                odds_ratio = odds_ratio_dict[filename_res.group(3)]
                last_block = int(filename_res.group(4))
                tag = filename_res.group(5)
                with open(x, 'r') as infile:
                    for line in infile:
                        h2_regex = r"Total Liability scale h2:\s+(-?\d+\.\d+)\s+\((\d+\.\d+)\)"
                        if re.search(h2_regex, line):
                            h2_res = re.search(h2_regex, line)
                            h2 = float(h2_res.group(1))
                            se = float(h2_res.group(2))
                            break
                outfile.write("%d\t%d\t%s\t0-%d\t%s\t%.4f\t%.4f\n" % (ncases, ncontrols, odds_ratio, last_block, tag, h2, se))

rule compile_swept_whole_genome_h2_parameters:
    input:
        [[[["results/ldsc/h2/whole_genome/%d_%d/1-%s0:%d_%s.log" % (i, i, z, j, k) for i in [10000, 50000, 100000, 250000]] for j in [9, 19, 39, 79, 99, 119]] for z in ["s", "m", "l", "v"]] for k in ["a", "b", "c", "d", "e"]]
    output:
        "results/ldsc/h2/compiled/whole_genome_sweep.tsv"
    run:
        with open(output[0], 'w') as outfile:
            outfile.write("ncases\tncontrols\todds_ratio\tblocks\ttag\th2\tse\n")
            for x in input:
                # TODO need to generalise this regex
                res = re.match("results/ldsc/h2/whole_genome/(\d+)_(\d+)/(\d+)-([smlv])0:(\d+)_([abcde]).log", x)
                ncases = int(res.group(1))
                ncontrols = int(res.group(2))
                chr = int(res.group(3))
                odds_ratio = odds_ratio_dict[res.group(4)]
                last_block = int(res.group(5))
                tag = res.group(6)
                with open(x, 'r') as infile:
                    for line in infile:
                        h2_regex = r"Total Liability scale h2: (-?\d+\.\d+)\s+\((\d+\.\d+)\)"
                        if re.search(h2_regex, line):
                            res = re.search(h2_regex, line)
                            h2 = float(res.group(1))
                            se = float(res.group(2))
                            break
                outfile.write("%d\t%d\t%.2f\t0-%d\t%s\t%.4f\t%.4f\n" % (ncases, ncontrols, odds_ratio, last_block, tag, h2, se))

rule compile_swept_whole_genome_B0_1_h2_parameters:
    input:
        [[[["results/ldsc/h2/whole_genome/%d_%d/B0_1/1-%s0:%d_%s.log" % (i, i, z, j, k) for i in [10000, 50000, 100000, 250000]] for j in [9, 19, 39, 79, 99, 119]] for z in ["s", "m", "l", "v"]] for k in ["a", "b", "c", "d", "e"]]
    output:
        "results/ldsc/h2/compiled/whole_genome_B0_1_sweep.tsv"
    run:
        with open(output[0], 'w') as outfile:
            outfile.write("ncases\tncontrols\todds_ratio\tblocks\ttag\th2\tse\n")
            for x in input:
                # TODO need to generalise this regex
                res = re.match("results/ldsc/h2/whole_genome/(\d+)_(\d+)/B0_1/(\d+)-([smlv])0:(\d+)_([abcde]).log", x)
                ncases = int(res.group(1))
                ncontrols = int(res.group(2))
                chr = int(res.group(3))
                odds_ratio = odds_ratio_dict[res.group(4)]
                last_block = int(res.group(5))
                tag = res.group(6)
                with open(x, 'r') as infile:
                    for line in infile:
                        h2_regex = r"Total Liability scale h2: (-?\d+\.\d+)\s+\((\d+\.\d+)\)"
                        if re.search(h2_regex, line):
                            res = re.search(h2_regex, line)
                            h2 = float(res.group(1))
                            se = float(res.group(2))
                            break
                outfile.write("%d\t%d\t%.2f\t0-%d\t%s\t%.4f\t%.4f\n" % (ncases, ncontrols, odds_ratio, last_block, tag, h2, se))

rule compile_swept_whole_genome_random_effect_h2_parameters:
    input:
        [[[["results/ldsc/h2/whole_genome/%d_%d/1-%s0:%d_%s.log" % (i, i, z, j, k) for i in [10000, 50000, 100000, 250000]] for j in [9, 19, 39, 79, 99, 119]] for z in ["r"]] for k in ["a", "b", "c", "d", "e"]]
    output:
        "results/ldsc/h2/compiled/whole_genome_random_effect_sweep.tsv"
    run:
        with open(output[0], 'w') as outfile:
            outfile.write("ncases\tncontrols\todds_ratio\tblocks\ttag\th2\tse\n")
            for x in input:
                res = re.match("results/ldsc/h2/whole_genome/(\d+)_(\d+)/(\d+)-([smlvr])0:(\d+)_([abcde]).log", x)
                ncases = int(res.group(1))
                ncontrols = int(res.group(2))
                chr = int(res.group(3))
                odds_ratio = odds_ratio_dict[res.group(4)]
                last_block = int(res.group(5))
                tag = res.group(6)
                with open(x, 'r') as infile:
                    for line in infile:
                        h2_regex = r"Total Liability scale h2: (-?\d+\.\d+)\s+\((\d+\.\d+)\)"
                        if re.search(h2_regex, line):
                            res = re.search(h2_regex, line)
                            h2 = float(res.group(1))
                            se = float(res.group(2))
                            break
                outfile.write("%d\t%d\t%s\t0-%d\t%s\t%.4f\t%.4f\n" % (ncases, ncontrols, odds_ratio, last_block, tag, h2, se))

rule compile_swept_whole_genome_random_effect_B0_1_h2_parameters:
    input:
        [[[["results/ldsc/h2/whole_genome/%d_%d/B0_1/1-%s0:%d_%s.log" % (i, i, z, j, k) for i in [10000, 50000, 100000, 250000]] for j in [9, 19, 39, 79, 99, 119]] for z in ["r"]] for k in ["a", "b", "c", "d", "e"]]
    output:
        "results/ldsc/h2/compiled/whole_genome_random_effect_B0_1_sweep.tsv"
    run:
        with open(output[0], 'w') as outfile:
            outfile.write("ncases\tncontrols\todds_ratio\tblocks\ttag\th2\tse\n")
            for x in input:
                res = re.match("results/ldsc/h2/whole_genome/(\d+)_(\d+)/B0_1/(\d+)-([smlvr])0:(\d+)_([abcde]).log", x)
                ncases = int(res.group(1))
                ncontrols = int(res.group(2))
                chr = int(res.group(3))
                odds_ratio = odds_ratio_dict[res.group(4)]
                last_block = int(res.group(5))
                tag = res.group(6)
                with open(x, 'r') as infile:
                    for line in infile:
                        h2_regex = r"Total Liability scale h2: (-?\d+\.\d+)\s+\((\d+\.\d+)\)"
                        if re.search(h2_regex, line):
                            res = re.search(h2_regex, line)
                            h2 = float(res.group(1))
                            se = float(res.group(2))
                            break
                outfile.write("%d\t%d\t%s\t0-%d\t%s\t%.4f\t%.4f\n" % (ncases, ncontrols, odds_ratio, last_block, tag, h2, se))

# NB: Assuming single causal variant atm
rule compile_odds_ratios_for_random_effect_simgwas_simulation:
    input:
        effect_block_files = get_effect_block_files
    output:
        "results/simgwas/simulated_sum_stats/compiled_effects/chr{ch}/{ncases}_{ncontrols}/{effect_blocks}.tsv"
    run:
        out_daf = pd.DataFrame()
        for x in input:
            with open(x, 'r') as sim_file:
                daf = pd.read_csv(x, sep = '\t', compression = 'gzip')
                out_daf = pd.concat([out_daf, daf[daf["chosen_or"] != 0]])
        out_daf.to_csv(output[0], sep = '\t', index = False, header = True)
