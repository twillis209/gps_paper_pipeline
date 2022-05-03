from collections import namedtuple
from string import ascii_lowercase
import random
import re
import pandas as pd
import os
from numpy import nan

block_daf = pd.read_csv("resources/ldetect/available_blocks.tsv", sep = "\t")

block_dict = {}

for i in range(1,23):
    block_dict[i] = list(block_daf[block_daf["chr"] == i][block_daf["available"] == True]["block"])

tags = list(ascii_lowercase[:20])

tag_pairs = [tags[i]+tags[i+1] for i in range(0, 19, 2)]

odds_ratio_dict = {"s": 1.05, "m": 1.2, 'l': 1.4, 'v': 2, 'r': 'random', 'n' : 1, 'i' : 1.1}

effect_size_dict = {'s': 'small', 'm': 'medium', 'l': 'large', 'v': 'vlarge', 'h': 'huge', 'r': 'random', 'i': 'intermediate'}

# 22
medium_effect_tuples = [f"m50_m50_m{x}" for x in [0, 10, 20, 30, 40, 50]]+[f"m25_m25_m{x}" for x in [0, 5, 10, 15, 20, 25]]+[f"m50_m50_m{x}" for x in [2, 5, 7, 12, 15, 17]]+[f"m25_m25_m{x}" for x in [2, 7, 12, 17]]

small_effect_tuples = [f"s400_s400_s{x}" for x in range(0, 401, 100)]+[f"s400_s400_s{x}" for x in range(50, 351, 100)]

# TODO need to fix this to remove the ncontrols < 10k jobs
medium_effect_sample_sizes = [(10000, 10000, 10000, 10000),
                (50000, 50000, 50000, 50000),
                (100000, 100000, 100000, 100000),
                (250000, 250000, 250000, 250000),
                (500, 10000, 500, 10000),
                (1000, 10000, 1000, 10000),
                (2000, 10000, 2000, 10000),
                (5000, 10000, 5000, 10000)]

small_effect_sample_sizes = [(500, 10000, 500, 10000),
                (1000, 10000, 1000, 10000),
                (2000, 10000, 2000, 10000),
                (5000, 10000, 5000, 10000),
                (10000, 10000, 10000, 10000),
                (50000, 50000, 50000, 50000),
                (100000, 100000, 100000, 100000),
                (250000, 250000, 250000, 250000)]

medium_effect_rg_estimate_files  = [f"results/ldsc/rg/whole_genome/randomised/{size[0]}_{size[1]}_{size[2]}_{size[3]}/{effect_tuple}_seed_%d_{tag_pair}.log" for size in medium_effect_sample_sizes for tag_pair in tag_pairs for effect_tuple in medium_effect_tuples]

medium_effect_rg_estimate_free_h2_free_rg_files = [f"results/ldsc/rg/free_h2_free_rg_intercept/whole_genome/randomised/{size[0]}_{size[1]}_{size[2]}_{size[3]}/{effect_tuple}_seed_%d_{tag_pair}.log" for size in medium_effect_sample_sizes for tag_pair in tag_pairs for effect_tuple in medium_effect_tuples]

medium_effect_rg_estimate_fixed_h2_free_rg_files = [f"results/ldsc/rg/fixed_h2_free_rg_intercept/whole_genome/randomised/{size[0]}_{size[1]}_{size[2]}_{size[3]}/{effect_tuple}_seed_%d_{tag_pair}.log" for size in medium_effect_sample_sizes for tag_pair in tag_pairs for effect_tuple in medium_effect_tuples]

medium_effect_rg_estimate_fixed_h2_fixed_rg_files = [f"results/ldsc/rg/fixed_h2_fixed_rg_intercept/whole_genome/randomised/{size[0]}_{size[1]}_{size[2]}_{size[3]}/{effect_tuple}_seed_%d_{tag_pair}.log" for size in medium_effect_sample_sizes for tag_pair in tag_pairs for effect_tuple in medium_effect_tuples]

medium_effect_theoretical_rg_files = [f"results/ldsc/rg/whole_genome/randomised/theoretical_rg/{size[0]}_{size[1]}_{size[2]}_{size[3]}/{effect_tuple}_seed_%d_{tag_pair}_theo_rg.tsv" for size in medium_effect_sample_sizes for tag_pair in tag_pairs for effect_tuple in medium_effect_tuples]

medium_effect_rg_estimate_free_h2_fixed_rg_files = [f"results/ldsc/rg/free_h2_fixed_rg_intercept/whole_genome/randomised/{size[0]}_{size[1]}_{size[2]}_{size[3]}/{effect_tuple}_seed_%d_{tag_pair}.log" for size in medium_effect_sample_sizes for tag_pair in tag_pairs for effect_tuple in medium_effect_tuples]

medium_effect_sumher_rg_files = [f"results/simgwas/ldak/ldak-thin/rg/{size[0]}_{size[1]}_{size[2]}_{size[3]}/{effect_tuple}_seed_%d_{tag_pair}.cors" for size in medium_effect_sample_sizes for tag_pair in tag_pairs for effect_tuple in medium_effect_tuples]

medium_effect_gps_files = [f"results/gps/simgwas/randomised/window_1000kb_step_50/{size[0]}_{size[1]}_{size[2]}_{size[3]}/3000_permutations/{effect_tuple}_seed_%d_tags_{tag_pair}_gps_pvalue.tsv" for size in medium_effect_sample_sizes for tag_pair in tag_pairs for effect_tuple in medium_effect_tuples]

medium_effect_hoeffdings_files = [f"results/hoeffdings/simgwas/randomised/window_1000kb_step_50/{size[0]}_{size[1]}_{size[2]}_{size[3]}/{effect_tuple}_seed_%d_tags_{tag_pair}_hoeffdings.tsv" for size in medium_effect_sample_sizes for tag_pair in tag_pairs for effect_tuple in medium_effect_tuples]

seeded_medium_effect_rg_estimate_free_h2_free_rg_files = [x[0] % x[1] for x in zip(medium_effect_rg_estimate_free_h2_free_rg_files, range(540, len(medium_effect_rg_estimate_free_h2_free_rg_files)))]

seeded_medium_effect_rg_estimate_fixed_h2_free_rg_files = [x[0] % x[1] for x in zip(medium_effect_rg_estimate_fixed_h2_free_rg_files, range(540, len(medium_effect_rg_estimate_fixed_h2_free_rg_files)))]

seeded_medium_effect_rg_estimate_fixed_h2_fixed_rg_files = [x[0] % x[1] for x in zip(medium_effect_rg_estimate_fixed_h2_fixed_rg_files, range(540, len(medium_effect_rg_estimate_fixed_h2_fixed_rg_files)))]

seeded_medium_effect_rg_estimate_free_h2_fixed_rg_files = [x[0] % x[1] for x in zip(medium_effect_rg_estimate_free_h2_fixed_rg_files, range(540, len(medium_effect_rg_estimate_free_h2_fixed_rg_files)))]

seeded_medium_effect_theoretical_rg_files = [x[0] % x[1] for x in zip(medium_effect_theoretical_rg_files, range(540, len(medium_effect_theoretical_rg_files)))]

seeded_medium_effect_sumher_rg_files = [x[0] % x[1] for x in zip(medium_effect_sumher_rg_files, range(540, len(medium_effect_sumher_rg_files)))]

seeded_medium_effect_gps_files = [x[0] % x[1] for x in zip(medium_effect_gps_files, range(540, len(medium_effect_gps_files)))]

seeded_medium_effect_hoeffdings_files = [x[0] % x[1] for x in zip(medium_effect_hoeffdings_files, range(540, len(medium_effect_hoeffdings_files)))]

small_effect_rg_estimate_free_h2_free_rg_files = [f"results/ldsc/rg/free_h2_free_rg_intercept/whole_genome/randomised/{size[0]}_{size[1]}_{size[2]}_{size[3]}/{effect_tuple}_seed_%d_{tag_pair}.log" for size in small_effect_sample_sizes for tag_pair in tag_pairs for effect_tuple in small_effect_tuples]

small_effect_rg_estimate_fixed_h2_free_rg_files = [f"results/ldsc/rg/fixed_h2_free_rg_intercept/whole_genome/randomised/{size[0]}_{size[1]}_{size[2]}_{size[3]}/{effect_tuple}_seed_%d_{tag_pair}.log" for size in small_effect_sample_sizes for tag_pair in tag_pairs for effect_tuple in small_effect_tuples]

small_effect_rg_estimate_fixed_h2_fixed_rg_files = [f"results/ldsc/rg/fixed_h2_fixed_rg_intercept/whole_genome/randomised/{size[0]}_{size[1]}_{size[2]}_{size[3]}/{effect_tuple}_seed_%d_{tag_pair}.log" for size in small_effect_sample_sizes for tag_pair in tag_pairs for effect_tuple in small_effect_tuples]

small_effect_rg_estimate_free_h2_fixed_rg_files = [f"results/ldsc/rg/free_h2_fixed_rg_intercept/whole_genome/randomised/{size[0]}_{size[1]}_{size[2]}_{size[3]}/{effect_tuple}_seed_%d_{tag_pair}.log" for size in small_effect_sample_sizes for tag_pair in tag_pairs for effect_tuple in small_effect_tuples]

small_effect_theoretical_rg_files = [f"results/ldsc/rg/whole_genome/randomised/theoretical_rg/{size[0]}_{size[1]}_{size[2]}_{size[3]}/{effect_tuple}_seed_%d_{tag_pair}_theo_rg.tsv" for size in small_effect_sample_sizes for tag_pair in tag_pairs for effect_tuple in small_effect_tuples]

small_effect_sumher_rg_files = [f"results/simgwas/ldak/ldak-thin/rg/{size[0]}_{size[1]}_{size[2]}_{size[3]}/{effect_tuple}_seed_%d_{tag_pair}.cors" for size in small_effect_sample_sizes for tag_pair in tag_pairs for effect_tuple in small_effect_tuples]

small_effect_gps_files = [f"results/gps/simgwas/randomised/window_1000kb_step_50/{size[0]}_{size[1]}_{size[2]}_{size[3]}/3000_permutations/{effect_tuple}_seed_%d_tags_{tag_pair}_gps_pvalue.tsv" for size in small_effect_sample_sizes for tag_pair in tag_pairs for effect_tuple in small_effect_tuples]

small_effect_hoeffdings_files = [f"results/hoeffdings/simgwas/randomised/window_1000kb_step_50/{size[0]}_{size[1]}_{size[2]}_{size[3]}/{effect_tuple}_seed_%d_tags_{tag_pair}_hoeffdings.tsv" for size in small_effect_sample_sizes for tag_pair in tag_pairs for effect_tuple in small_effect_tuples]

seeded_small_effect_rg_estimate_free_h2_free_rg_files = [x[0] % x[1] for x in zip(small_effect_rg_estimate_free_h2_free_rg_files, range(6000, 6000+len(small_effect_rg_estimate_free_h2_free_rg_files)))]

seeded_small_effect_rg_estimate_fixed_h2_free_rg_files = [x[0] % x[1] for x in zip(small_effect_rg_estimate_fixed_h2_free_rg_files, range(6000, 6000+len(small_effect_rg_estimate_fixed_h2_free_rg_files)))]

seeded_small_effect_rg_estimate_fixed_h2_fixed_rg_files = [x[0] % x[1] for x in zip(small_effect_rg_estimate_fixed_h2_fixed_rg_files, range(6000, 6000+len(small_effect_rg_estimate_fixed_h2_fixed_rg_files)))]

seeded_small_effect_rg_estimate_free_h2_fixed_rg_files = [x[0] % x[1] for x in zip(small_effect_rg_estimate_free_h2_fixed_rg_files, range(6000, 6000+len(small_effect_rg_estimate_free_h2_fixed_rg_files)))]

seeded_small_effect_theoretical_rg_files = [x[0] % x[1] for x in zip(small_effect_theoretical_rg_files, range(6000, 6000+len(small_effect_theoretical_rg_files)))]

seeded_small_effect_sumher_rg_files = [x[0] % x[1] for x in zip(small_effect_sumher_rg_files, range(6000, 6000+len(small_effect_sumher_rg_files)))]

seeded_small_effect_gps_files = [x[0] % x[1] for x in zip(small_effect_gps_files, range(6000, 6000+len(small_effect_gps_files)))]

seeded_small_effect_hoeffdings_files = [x[0] % x[1] for x in zip(small_effect_hoeffdings_files, range(6000, 6000+len(small_effect_hoeffdings_files)))]

exec(open('workflow/rules/simgwas_experiments_functions.smk').read())

#rule compile_small_rg_results:
#infiles = [y for y in seeded_small_effect_rg_estimate_free_h2_free_rg_files+
#         seeded_small_effect_rg_estimate_fixed_h2_free_rg_files+
#         seeded_small_effect_rg_estimate_fixed_h2_fixed_rg_files+
#         seeded_small_effect_rg_estimate_free_h2_fixed_rg_files if os.path.exists(y)]
#outfiles = "results/ldsc/rg/whole_genome/randomised/compiled_s_rg_results.tsv"
#
#compile_rg_results(infiles, outfiles)
#
##rule compile_medium_rg_results:
#infiles =  [y for y in seeded_medium_effect_rg_estimate_free_h2_free_rg_files+
#         seeded_medium_effect_rg_estimate_fixed_h2_free_rg_files+
#        seeded_medium_effect_rg_estimate_fixed_h2_fixed_rg_files+
#        seeded_medium_effect_rg_estimate_free_h2_fixed_rg_files if os.path.exists(y)]
#outfiles =  "results/ldsc/rg/whole_genome/randomised/compiled_m_rg_results.tsv"
#
#compile_rg_results(infiles, outfiles)

#Wildcards = namedtuple("Wildcards", "ncases_A ncontrols_A ncases_B ncontrols_B tag_pair effect_blocks_A effect_blocks_B shared_effect_blocks seed")
