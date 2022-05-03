import re
import os
from numpy import nan

include: 'simgwas_experiments_functions.smk'

tag_pairs = [tags[i]+tags[i+1] for i in range(0, 19, 2)]

odds_ratio_dict = {"s": 1.05, "m": 1.2, 'l': 1.4, 'v': 2, 'r': 'random', 'n' : 1, 'i' : 1.1}

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

rule run_small_effect_jobs:
    input:
        seeded_small_effect_rg_estimate_free_h2_free_rg_files,
        seeded_small_effect_rg_estimate_fixed_h2_free_rg_files,
        seeded_small_effect_rg_estimate_fixed_h2_fixed_rg_files,
        seeded_small_effect_rg_estimate_free_h2_fixed_rg_files,
        seeded_small_effect_gps_files,
        seeded_small_effect_hoeffdings_files,
        seeded_small_effect_sumher_rg_files

rule run_medium_effect_jobs:
    input:
        seeded_medium_effect_rg_estimate_free_h2_free_rg_files,
        seeded_medium_effect_rg_estimate_fixed_h2_free_rg_files,
        seeded_medium_effect_rg_estimate_fixed_h2_fixed_rg_files,
        seeded_medium_effect_rg_estimate_free_h2_fixed_rg_files,
        seeded_medium_effect_sumher_rg_files,
        seeded_medium_effect_gps_files,
        seeded_medium_effect_hoeffdings_files

rule compile_small_rg_results:
    input:
        [y for y in seeded_small_effect_rg_estimate_free_h2_free_rg_files+
         seeded_small_effect_rg_estimate_fixed_h2_free_rg_files+
         seeded_small_effect_rg_estimate_fixed_h2_fixed_rg_files+
         seeded_small_effect_rg_estimate_free_h2_fixed_rg_files if os.path.exists(y)]
    output:
        "results/ldsc/rg/whole_genome/randomised/compiled_s_rg_results.tsv"
    run:
        compile_rg_results(input, output)

rule compile_medium_rg_results:
    input:
        [y for y in seeded_medium_effect_rg_estimate_free_h2_free_rg_files+
         seeded_medium_effect_rg_estimate_fixed_h2_free_rg_files+
        seeded_medium_effect_rg_estimate_fixed_h2_fixed_rg_files+
        seeded_medium_effect_rg_estimate_free_h2_fixed_rg_files if os.path.exists(y)]
    output:
        "results/ldsc/rg/whole_genome/randomised/compiled_m_rg_results.tsv"
    run:
        compile_rg_results(input, output)

rule compile_small_effect_theoretical_rg:
    input:
        seeded_small_effect_theoretical_rg_files
    output:
        compiled_rg_file = "results/ldsc/rg/whole_genome/randomised/theoretical_rg/compiled_s_theoretical_rg.tsv"
    run:
        compile_theoretical_rg_results(input, output)

rule compile_medium_effect_theoretical_rg:
    input:
        seeded_medium_effect_theoretical_rg_files
    output:
        compiled_rg_file = "results/ldsc/rg/whole_genome/randomised/theoretical_rg/compiled_m_theoretical_rg.tsv"
    run:
        compile_theoretical_rg_results(input, output)

rule compile_small_effect_sumher_results:
    input:
        seeded_small_effect_sumher_rg_files
    output:
        compiled_rg_file = "results/simgwas/ldak/ldak-thin/rg/compiled_s_sumher_rg.tsv"
    run:
        compile_sumher_results(input, output)

rule compile_medium_effect_sumher_results:
    input:
        seeded_medium_effect_sumher_rg_files
    output:
        compiled_rg_file = "results/simgwas/ldak/ldak-thin/rg/compiled_m_sumher_rg.tsv"
    run:
        compile_sumher_results(input, output)

rule compile_small_gps_results:
    input:
        [y for y in seeded_small_effect_gps_files if os.path.exists(y)]
    output:
        "results/gps/simgwas/randomised/window_1000kb_step_50/compiled_s_gps_results.tsv"
    run:
        compile_gps_results(input, output)

rule compile_medium_gps_results:
    input:
        ancient([y for y in seeded_medium_effect_gps_files if os.path.exists(y)])
    output:
        "results/gps/simgwas/randomised/window_1000kb_step_50/compiled_m_gps_results.tsv"
    run:
        compile_gps_results(input, output)

rule compile_small_hoeffdings_results:
    input:
        [y for y in seeded_small_effect_hoeffdings_files if os.path.exists(y)]
    output:
        "results/hoeffdings/simgwas/randomised/window_1000kb_step_50/compiled_s_hoeffdings_results.tsv"
    run:
        compile_hoeffdings_results(input, output)

rule compile_medium_hoeffdings_results:
    input:
        ancient([y for y in seeded_medium_effect_hoeffdings_files if os.path.exists(y)])
    output:
        "results/hoeffdings/simgwas/randomised/window_1000kb_step_50/compiled_m_hoeffdings_results.tsv"
    run:
        compile_hoeffdings_results(input, output)
