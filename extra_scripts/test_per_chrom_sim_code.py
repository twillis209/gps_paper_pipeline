from collections import namedtuple
import random
import re
import pandas as pd

Wildcards = namedtuple("Wildcards", "seed chr ncases_A ncontrols_A ncases_B ncontrols_B effect_blocks_A effect_blocks_B shared_effect_blocks")

wildcards = Wildcards(1, "chr1", 10000, 10000, 10000, 10000, "s5", "s5", "s1")

exec(open('workflow/rules/simgwas/simgwas_randomised_one_chrom_functions.py', 'r').read())

effect_size_dict = {"s": "small", "m": "medium", "l": "large", "v": "vlarge", "h": "huge", "r": "random", "i": "infinitesimal", "t": "tiny"}

cv_per_block_dict = {"small": 1, "medium": 1, "large": 1, "vlarge": 1, "huge": 1, "tiny": 2, "null": 0}

cv_index_dict = {1: [500], 2: [250, 750]}

odds_ratio_dict = {"small": 1.05, "medium": 1.2, "null": 1, "tiny": 1.02}

block_daf = pd.read_csv("resources/ldetect/available_blocks.tsv", sep = "\t")

block_daf = block_daf.query('available == True')

get_randomised_chrom_block_tuples_for_chrom_for_pair(wildcards)
get_randomised_block_files_for_chrom_for_pair(wildcards)
