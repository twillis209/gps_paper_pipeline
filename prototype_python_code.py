from collections import namedtuple
from string import ascii_lowercase
import random
import re
import pandas as pd

block_daf = pd.read_csv("resources/ldetect/available_blocks.tsv", sep = "\t")

block_dict = {}

for i in range(1,23):
    block_dict[i] = list(block_daf[block_daf["chr"] == i][block_daf["available"] == True]["block"])

tags = list(ascii_lowercase[:20])

tag_pairs = [tags[i]+tags[i+1] for i in range(0, 19, 2)]

odds_ratio_dict = {"s": 1.05, "m": 1.2, 'l': 1.4, 'v': 2, 'r': 'random', 'n' : 1, 'i' : 1.1}

effect_size_dict = {'s': 'small', 'm': 'medium', 'l': 'large', 'v': 'vlarge', 'h': 'huge', 'r': 'random', 'i': 'intermediate'}

medium_effect_tuples = [f"m50_m50_m{x}" for x in [0, 10, 20, 30, 40, 50]]+[f"m25_m25_m{x}" for x in [0, 5, 10, 15, 20, 25]]

small_effect_tuples = [f"s400_s400_s{x}" for x in range(0, 101, 10)]

sample_sizes = [1000, 5000, 10000, 50000, 100000, 250000]

Wildcards = namedtuple("Wildcards", "ncases_A ncontrols_A ncases_B ncontrols_B tag_pair effect_blocks_A effect_blocks_B shared_effect_blocks seed")

medium_effect_rg_estimate_files = [f"results/ldsc/rg/whole_genome/randomised/{size}_{size}_{size}_{size}/{effect_tuple}_seed_%d_{tag_pair}.log" for size in sample_sizes for tag_pair in tag_pairs for effect_tuple in medium_effect_tuples]

medium_effect_rg_estimate_tuples = [(size, tag_pair, effect_tuple) for size in sample_sizes for tag_pair in tag_pairs for effect_tuple in medium_effect_tuples]

medium_effect_rg_estimate_seeded_tuples = [(x[0][0], x[0][1], x[0][2].split('_')[0], x[0][2].split('_')[1], x[0][2].split('_')[2], x[1]) for x in zip(medium_effect_rg_estimate_tuples, range(100, 100+len(medium_effect_rg_estimate_files)))]

#medium_effect_rg_estimates_seeded_files = [x[0] % x[1] for x in zip(medium_effect_rg_estimate_files, range(100, 100+len(medium_effect_rg_estimate_files)))]

wildcards_list = []

for i,x in enumerate(medium_effect_rg_estimate_seeded_tuples):
    wildcards_list.append(Wildcards(x[0], x[0], x[0], x[0], x[1], x[2], x[3], x[4], x[5]))


def get_randomised_chrom_block_tuples_for_pair(wildcards):
    random.seed(wildcards.seed)

    shared_chrom_block_nos = []

    if wildcards.shared_effect_blocks != 'null':
        block_match = re.match('([smlvh])(\d+)', wildcards.shared_effect_blocks)

        if not block_match:
            raise ValueError("Invalid block format: %s" % wildcards.shared_effect_blocks)

        effect = effect_size_dict[block_match.group(1)]
        no_of_shared_blocks = int(block_match.group(2))

        i = 0

        while i < max(no_of_shared_blocks, 0):
            chrom = random.choice(list(block_dict.keys()))
            block_no = random.choice(block_dict[chrom])

            if (chrom, block_no) not in shared_chrom_block_nos:
                shared_chrom_block_nos.append((chrom, block_no))
                i += 1

    a_chrom_block_nos = []

    if wildcards.effect_blocks_A != 'null':
        block_match_a = re.match('([smlvh])(\d+)', wildcards.effect_blocks_A)

        if not block_match_a:
            raise ValueError("Invalid block format: %s" % wildcards.effect_blocks_A)

        effect_a = effect_size_dict[block_match_a.group(1)]
        no_of_blocks_a = int(block_match_a.group(2))

        i = 0

        while i < max(no_of_blocks_a-no_of_shared_blocks, 0):
            chrom = random.choice(list(block_dict.keys()))
            block_no = random.choice(block_dict[chrom])

            if (chrom, block_no) not in shared_chrom_block_nos and (chrom, block_no) not in a_chrom_block_nos:
                a_chrom_block_nos.append((chrom, block_no))
                i += 1

    b_chrom_block_nos = []

    if wildcards.effect_blocks_B != 'null':
        block_match_b = re.match('([smlvh])(\d+)', wildcards.effect_blocks_B)

        if not block_match_b:
            raise ValueError("Invalid block format: %s" % wildcards.effect_blocks_B)

        effect_b = effect_size_dict[block_match_b.group(1)]
        no_of_blocks_b = int(block_match_b.group(2))

        i = 0

        while i < max(no_of_blocks_b-no_of_shared_blocks, 0):
            chrom = random.choice(list(block_dict.keys()))
            block_no = random.choice(block_dict[chrom])

            if (chrom, block_no) not in shared_chrom_block_nos and (chrom, block_no) not in a_chrom_block_nos and (chrom, block_no) not in b_chrom_block_nos:
                b_chrom_block_nos.append((chrom, block_no))
                i += 1

    return (shared_chrom_block_nos, a_chrom_block_nos, b_chrom_block_nos)

def get_randomised_block_files_for_pair(wildcards):
    shared_chrom_block_nos, a_chrom_block_nos, b_chrom_block_nos = get_randomised_chrom_block_tuples_for_pair(wildcards)

    block_files = []

    a_block_files = []
    b_block_files = []

    for x in shared_chrom_block_nos:
        block_files.append(f"results/simgwas/simulated_sum_stats/chr{x[0]}/block_sum_stats/{effect}/{wildcards.ncases_A}_{wildcards.ncontrols_A}/block_{x[1]}_sum_stats.tsv.gz")
        a_block_files.append(f"results/simgwas/simulated_sum_stats/chr{x[0]}/block_sum_stats/{effect}/{wildcards.ncases_A}_{wildcards.ncontrols_A}/block_{x[1]}_sum_stats.tsv.gz")
        b_block_files.append(f"results/simgwas/simulated_sum_stats/chr{x[0]}/block_sum_stats/{effect}/{wildcards.ncases_B}_{wildcards.ncontrols_B}/block_{x[1]}_sum_stats.tsv.gz")
        if wildcards.ncases_A != wildcards.ncases_B or wildcards.ncontrols_A != wildcards.ncontrols_B:
            block_files.append(f"results/simgwas/simulated_sum_stats/chr{x[0]}/block_sum_stats/{effect}/{wildcards.ncases_B}_{wildcards.ncontrols_B}/block_{x[1]}_sum_stats.tsv.gz")

    for x in a_chrom_block_nos:
        block_files.append(f"results/simgwas/simulated_sum_stats/chr{x[0]}/block_sum_stats/{effect_a}/{wildcards.ncases_A}_{wildcards.ncontrols_A}/block_{x[1]}_sum_stats.tsv.gz")
        a_block_files.append(f"results/simgwas/simulated_sum_stats/chr{x[0]}/block_sum_stats/{effect_a}/{wildcards.ncases_A}_{wildcards.ncontrols_A}/block_{x[1]}_sum_stats.tsv.gz")

    for x in b_chrom_block_nos:
        b_block_files.append(f"results/simgwas/simulated_sum_stats/chr{x[0]}/block_sum_stats/{effect_b}/{wildcards.ncases_B}_{wildcards.ncontrols_B}/block_{x[1]}_sum_stats.tsv.gz")
        if f"results/simgwas/simulated_sum_stats/chr{x[0]}/block_sum_stats/{effect_b}/{wildcards.ncases_B}_{wildcards.ncontrols_B}/block_{x[1]}_sum_stats.tsv.gz" not in block_files:
            block_files.append(f"results/simgwas/simulated_sum_stats/chr{x[0]}/block_sum_stats/{effect_b}/{wildcards.ncases_B}_{wildcards.ncontrols_B}/block_{x[1]}_sum_stats.tsv.gz")

    for chrom in block_dict.keys():
        for block_no in block_dict[chrom]:
            if (chrom, block_no) not in a_chrom_block_nos and (chrom, block_no) not in shared_chrom_block_nos:
                block_files.append(f"results/simgwas/simulated_sum_stats/chr{chrom}/block_sum_stats/null/{wildcards.ncases_A}_{wildcards.ncontrols_A}/block_{block_no}_sum_stats.tsv.gz")
                a_block_files.append(f"results/simgwas/simulated_sum_stats/chr{chrom}/block_sum_stats/null/{wildcards.ncases_A}_{wildcards.ncontrols_A}/block_{block_no}_sum_stats.tsv.gz")
            if (chrom, block_no) not in b_chrom_block_nos and (chrom, block_no) not in shared_chrom_block_nos:
                b_block_files.append(f"results/simgwas/simulated_sum_stats/chr{chrom}/block_sum_stats/null/{wildcards.ncases_B}_{wildcards.ncontrols_B}/block_{block_no}_sum_stats.tsv.gz")
                if f"results/simgwas/simulated_sum_stats/chr{chrom}/block_sum_stats/null/{wildcards.ncases_B}_{wildcards.ncontrols_B}/block_{block_no}_sum_stats.tsv.gz" not in block_files:
                    block_files.append(f"results/simgwas/simulated_sum_stats/chr{chrom}/block_sum_stats/null/{wildcards.ncases_B}_{wildcards.ncontrols_B}/block_{block_no}_sum_stats.tsv.gz")

    return (block_files, a_block_files, b_block_files)

def broken_get_randomised_chrom_block_tuples_for_pair(wildcards):
    random.seed(wildcards.seed)

    shared_chrom_block_nos = []

    if wildcards.shared_effect_blocks != 'null':
        block_match = re.match('([smlvh])(\d+)', wildcards.shared_effect_blocks)

        if not block_match:
            raise ValueError("Invalid block format: %s" % wildcards.shared_effect_blocks)

        effect = effect_size_dict[block_match.group(1)]
        no_of_shared_blocks = int(block_match.group(2))

        for i in range(no_of_shared_blocks):
            chrom = random.choice(list(block_dict.keys()))
            block_no = random.choice(block_dict[chrom])

            if (chrom, block_no) not in shared_chrom_block_nos:
                shared_chrom_block_nos.append((chrom, block_no))

    a_chrom_block_nos = []

    if wildcards.effect_blocks_A != 'null':
        block_match_a = re.match('([smlvh])(\d+)', wildcards.effect_blocks_A)

        if not block_match_a:
            raise ValueError("Invalid block format: %s" % wildcards.effect_blocks_A)

        effect_a = effect_size_dict[block_match_a.group(1)]
        no_of_blocks_a = int(block_match_a.group(2))

        i = 0

        while i < max(no_of_blocks_a-no_of_shared_blocks, 0):
            chrom = random.choice(list(block_dict.keys()))
            block_no = random.choice(block_dict[chrom])

            if (chrom, block_no) not in shared_chrom_block_nos and (chrom, block_no) not in a_chrom_block_nos:
                a_chrom_block_nos.append((chrom, block_no))
                i += 1

    b_chrom_block_nos = []

    if wildcards.effect_blocks_B != 'null':
        block_match_b = re.match('([smlvh])(\d+)', wildcards.effect_blocks_B)

        if not block_match_b:
            raise ValueError("Invalid block format: %s" % wildcards.effect_blocks_B)

        effect_b = effect_size_dict[block_match_b.group(1)]
        no_of_blocks_b = int(block_match_b.group(2))

        i = 0

        while i < max(no_of_blocks_b-no_of_shared_blocks, 0):
            chrom = random.choice(list(block_dict.keys()))
            block_no = random.choice(block_dict[chrom])

            if (chrom, block_no) not in shared_chrom_block_nos and (chrom, block_no) not in a_chrom_block_nos and (chrom, block_no) not in b_chrom_block_nos:
                b_chrom_block_nos.append((chrom, block_no))
                i += 1

    return (shared_chrom_block_nos, a_chrom_block_nos, b_chrom_block_nos)

def broken_get_randomised_block_files_for_pair(wildcards):
    shared_chrom_block_nos, a_chrom_block_nos, b_chrom_block_nos = get_randomised_chrom_block_tuples_for_pair(wildcards)

    block_files = []

    a_block_files = []
    b_block_files = []

    for x in shared_chrom_block_nos:
        block_files.append(f"results/simgwas/simulated_sum_stats/chr{x[0]}/block_sum_stats/{effect}/{wildcards.ncases_A}_{wildcards.ncontrols_A}/block_{x[1]}_sum_stats.tsv.gz")
        a_block_files.append(f"results/simgwas/simulated_sum_stats/chr{x[0]}/block_sum_stats/{effect}/{wildcards.ncases_A}_{wildcards.ncontrols_A}/block_{x[1]}_sum_stats.tsv.gz")
        b_block_files.append(f"results/simgwas/simulated_sum_stats/chr{x[0]}/block_sum_stats/{effect}/{wildcards.ncases_B}_{wildcards.ncontrols_B}/block_{x[1]}_sum_stats.tsv.gz")
        if wildcards.ncases_A != wildcards.ncases_B or wildcards.ncontrols_A != wildcards.ncontrols_B:
            block_files.append(f"results/simgwas/simulated_sum_stats/chr{x[0]}/block_sum_stats/{effect}/{wildcards.ncases_B}_{wildcards.ncontrols_B}/block_{x[1]}_sum_stats.tsv.gz")

    for x in a_chrom_block_nos:
        block_files.append(f"results/simgwas/simulated_sum_stats/chr{x[0]}/block_sum_stats/{effect_a}/{wildcards.ncases_A}_{wildcards.ncontrols_A}/block_{x[1]}_sum_stats.tsv.gz")
        a_block_files.append(f"results/simgwas/simulated_sum_stats/chr{x[0]}/block_sum_stats/{effect_a}/{wildcards.ncases_A}_{wildcards.ncontrols_A}/block_{x[1]}_sum_stats.tsv.gz")

    for x in b_chrom_block_nos:
        b_block_files.append(f"results/simgwas/simulated_sum_stats/chr{x[0]}/block_sum_stats/{effect_b}/{wildcards.ncases_B}_{wildcards.ncontrols_B}/block_{x[1]}_sum_stats.tsv.gz")
        if f"results/simgwas/simulated_sum_stats/chr{x[0]}/block_sum_stats/{effect_b}/{wildcards.ncases_B}_{wildcards.ncontrols_B}/block_{x[1]}_sum_stats.tsv.gz" not in block_files:
            block_files.append(f"results/simgwas/simulated_sum_stats/chr{x[0]}/block_sum_stats/{effect_b}/{wildcards.ncases_B}_{wildcards.ncontrols_B}/block_{x[1]}_sum_stats.tsv.gz")

    for chrom in block_dict.keys():
        for block_no in block_dict[chrom]:
            if (chrom, block_no) not in a_chrom_block_nos and (chrom, block_no) not in shared_chrom_block_nos:
                block_files.append(f"results/simgwas/simulated_sum_stats/chr{chrom}/block_sum_stats/null/{wildcards.ncases_A}_{wildcards.ncontrols_A}/block_{block_no}_sum_stats.tsv.gz")
                a_block_files.append(f"results/simgwas/simulated_sum_stats/chr{chrom}/block_sum_stats/null/{wildcards.ncases_A}_{wildcards.ncontrols_A}/block_{block_no}_sum_stats.tsv.gz")
            if (chrom, block_no) not in b_chrom_block_nos and (chrom, block_no) not in shared_chrom_block_nos:
                b_block_files.append(f"results/simgwas/simulated_sum_stats/chr{chrom}/block_sum_stats/null/{wildcards.ncases_B}_{wildcards.ncontrols_B}/block_{block_no}_sum_stats.tsv.gz")
                if f"results/simgwas/simulated_sum_stats/chr{chrom}/block_sum_stats/null/{wildcards.ncases_B}_{wildcards.ncontrols_B}/block_{block_no}_sum_stats.tsv.gz" not in block_files:
                    block_files.append(f"results/simgwas/simulated_sum_stats/chr{chrom}/block_sum_stats/null/{wildcards.ncases_B}_{wildcards.ncontrols_B}/block_{block_no}_sum_stats.tsv.gz")

    return (block_files, a_block_files, b_block_files)

#with open('intended_shared_blocks_vs_actual_blocks.tsv', 'w') as outfile, open('files_to_delete.txt', 'w') as files_outfile:
#    for i,x in enumerate(wildcards_list):
#        if int(x.shared_effect_blocks[1:]) != len(broken_get_randomised_chrom_block_tuples_for_pair(x)[0]):
#            rg_file_path = f"results/ldsc/rg/whole_genome/randomised/{x.ncases_A}_{x.ncontrols_A}_{x.ncases_B}_{x.ncontrols_B}/{x.effect_blocks_A}_{x.effect_blocks_B}_{x.shared_effect_blocks}_seed_{x.seed}_{x.tag_pair}.log"
#            gps_file_path = f"results/gps/simgwas/randomised/window_1000kb_step_50/{x.ncases_A}_{x.ncontrols_A}_{x.ncases_B}_{x.ncontrols_B}/3000_permutations/{x.effect_blocks_A}_{x.effect_blocks_B}_{x.shared_effect_blocks}_seed_{x.seed}_tags_{x.tag_pair}_gps_pvalue.tsv"
#            hoeffdings_file_path = f"results/hoeffdings/simgwas/randomised/window_1000kb_step_50/{x.ncases_A}_{x.ncontrols_A}_{x.ncases_B}_{x.ncontrols_B}/{x.effect_blocks_A}_{x.effect_blocks_B}_{x.shared_effect_blocks}_seed_{x.seed}_tags_{x.tag_pair}_hoeffdings.tsv"
#            permute_file_path = f"results/gps/simgwas/randomised/window_1000kb_step_50/{x.ncases_A}_{x.ncontrols_A}_{x.ncases_B}_{x.ncontrols_B}/3000_permutations/{x.effect_blocks_A}_{x.effect_blocks_B}_{x.shared_effect_blocks}_seed_{x.seed}_tags_{x.tag_pair}.tsv"
#            #outfile.write("%s\t%s\t%s\t%s\t%s\t%d\n" % (x.tag_pair, x.seed, x.effect_blocks_A[1:], x.effect_blocks_B[1:], x.shared_effect_blocks[1:]))
#            files_outfile.write("%s\n%s\n%s\n%s\n" % (rg_file_path, gps_file_path, hoeffdings_file_path, permute_file_path))

with open('fixed_shared_blocks_vs_actual_blocks.tsv', 'w') as outfile:
    for i,x in enumerate(wildcards_list):
        outfile.write("%s\t%s\t%s\t%s\t%s\t%d\n" % (x.tag_pair, x.seed, x.effect_blocks_A[1:], x.effect_blocks_B[1:], x.shared_effect_blocks[1:], len(get_randomised_chrom_block_tuples_for_pair(x)[0])))

