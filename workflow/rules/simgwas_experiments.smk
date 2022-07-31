import re
import os
from numpy import nan

include: 'simgwas_experiments_functions.py'

odds_ratio_dict = {"s": 1.05, "m": 1.2, 'l': 1.4, 'v': 2, 'r': 'random', 'n' : 1, 'i' : 1.1}

sample_sizes = [(500, 10000),
                (1000, 10000),
                (2000, 10000),
                (5000, 10000),
                (10000, 10000),
                (50000, 50000),
                (100000, 100000),
                (250000, 250000)]

rule run_all_block_simulations:
    input:
        block_files = get_all_block_files(sample_sizes)

rule write_out_simulation_parameters_file:
    output:
        "results/simgwas/simulation_parameters.tsv"
    params:
        sample_sizes = sample_sizes
    script: "../scripts/simgwas/write_out_simulation_parameters.py"

rule test_new_randomised_simgwas_code:
    input:
        "results/ldsc/simgwas/400_reps/randomised/10000_10000_10000_10000/s200_s200_s10/rg/fixed_h2_free_rg_intercept/seed_5_tags_3-4.log",
        "results/gps/simgwas/400_reps/randomised/10000_10000_10000_10000/s200_s200_s10/window_1000kb_step_50/3000_permutations/seed_5_tags_3-4_gps_pvalue.tsv",
        "results/hoeffdings/simgwas/400_reps/randomised/10000_10000_10000_10000/s200_s200_s10/window_1000kb_step_50/seed_5_tags_3-4_hoeffdings.tsv",
        "results/ldak/ldak-thin/simgwas/400_reps/randomised/rg/10000_10000_10000_10000/s200_s200_s10/seed_5_tags_3-4.cors"




"""
        pruned_sum_stats = "results/simgwas/simulated_sum_stats/whole_genome_sum_stats/400_reps/randomised/10000_10000_10000_10000/s200_s200_s10/window_1000kb_step_50/seed_5_pruned_sum_stats_tags_3-4.tsv.gz",
sum_stats_A = "results/simgwas/simulated_sum_stats/whole_genome_sum_stats/400_reps/randomised/10000_10000_10000_10000/s200_s200_s10/seed_5_sum_stats_A_tag_3_of_3-4.tsv.gz",
sum_stats_B = "results/simgwas/simulated_sum_stats/whole_genome_sum_stats/400_reps/randomised/10000_10000_10000_10000/s200_s200_s10/seed_5_sum_stats_B_tag_4_of_3-4.tsv.gz"


"""

"""
medium_effect_rg_estimate_files  = [f"results/ldsc/rg/whole_genome/randomised/{size[0]}_{size[1]}_{size[2]}_{size[3]}/{effect_tuple}_seed_%d_{tag_pair}.log" for size in medium_effect_sample_sizes for tag_pair in tag_pairs for effect_tuple in medium_effect_tuples]

medium_effect_rg_estimate_free_h2_free_rg_files = [f"results/ldsc/rg/free_h2_free_rg_intercept/whole_genome/randomised/{size[0]}_{size[1]}_{size[2]}_{size[3]}/{effect_tuple}_seed_%d_{tag_pair}.log" for size in medium_effect_sample_sizes for tag_pair in tag_pairs for effect_tuple in medium_effect_tuples]

medium_effect_rg_estimate_fixed_h2_free_rg_files = [f"results/ldsc/rg/fixed_h2_free_rg_intercept/whole_genome/randomised/{size[0]}_{size[1]}_{size[2]}_{size[3]}/{effect_tuple}_seed_%d_{tag_pair}.log" for size in medium_effect_sample_sizes for tag_pair in tag_pairs for effect_tuple in medium_effect_tuples]

medium_effect_rg_estimate_fixed_h2_fixed_rg_files = [f"results/ldsc/rg/fixed_h2_fixed_rg_intercept/whole_genome/randomised/{size[0]}_{size[1]}_{size[2]}_{size[3]}/{effect_tuple}_seed_%d_{tag_pair}.log" for size in medium_effect_sample_sizes for tag_pair in tag_pairs for effect_tuple in medium_effect_tuples]

medium_effect_theoretical_rg_files = [f"results/ldsc/rg/whole_genome/randomised/theoretical_rg/{size[0]}_{size[1]}_{size[2]}_{size[3]}/{effect_tuple}_seed_%d_{tag_pair}_theo_rg.tsv" for size in medium_effect_sample_sizes for tag_pair in tag_pairs for effect_tuple in medium_effect_tuples]

medium_effect_rg_estimate_free_h2_fixed_rg_files = [f"results/ldsc/rg/free_h2_fixed_rg_intercept/whole_genome/randomised/{size[0]}_{size[1]}_{size[2]}_{size[3]}/{effect_tuple}_seed_%d_{tag_pair}.log" for size in medium_effect_sample_sizes for tag_pair in tag_pairs for effect_tuple in medium_effect_tuples]

medium_effect_sumher_rg_files = [f"results/ldak/ldak-thin/1000g/rg/{size[0]}_{size[1]}_{size[2]}_{size[3]}/{effect_tuple}_seed_%d_{tag_pair}.progress" for size in medium_effect_sample_sizes for tag_pair in tag_pairs for effect_tuple in medium_effect_tuples]

medium_effect_gps_files = [f"results/gps/simgwas/randomised/window_1000kb_step_50/{size[0]}_{size[1]}_{size[2]}_{size[3]}/3000_permutations/{effect_tuple}_seed_%d_tags_{tag_pair}_gps_pvalue.tsv" for size in medium_effect_sample_sizes for tag_pair in tag_pairs for effect_tuple in medium_effect_tuples]

medium_effect_hoeffdings_files = [f"results/hoeffdings/simgwas/randomised/window_1000kb_step_50/{size[0]}_{size[1]}_{size[2]}_{size[3]}/{effect_tuple}_seed_%d_tags_{tag_pair}_hoeffdings.tsv" for size in medium_effect_sample_sizes for tag_pair in tag_pairs for effect_tuple in medium_effect_tuples]

seeded_medium_effect_rg_estimate_free_h2_free_rg_files = [x[0] % x[1] for x in zip(medium_effect_rg_estimate_free_h2_free_rg_files, range(540, 540+len(medium_effect_rg_estimate_free_h2_free_rg_files)))]

seeded_medium_effect_rg_estimate_fixed_h2_free_rg_files = [x[0] % x[1] for x in zip(medium_effect_rg_estimate_fixed_h2_free_rg_files, range(540, 540+len(medium_effect_rg_estimate_fixed_h2_free_rg_files)))]

seeded_medium_effect_rg_estimate_fixed_h2_fixed_rg_files = [x[0] % x[1] for x in zip(medium_effect_rg_estimate_fixed_h2_fixed_rg_files, range(540, 540+len(medium_effect_rg_estimate_fixed_h2_fixed_rg_files)))]

seeded_medium_effect_rg_estimate_free_h2_fixed_rg_files = [x[0] % x[1] for x in zip(medium_effect_rg_estimate_free_h2_fixed_rg_files, range(540, 540+len(medium_effect_rg_estimate_free_h2_fixed_rg_files)))]

seeded_medium_effect_theoretical_rg_files = [x[0] % x[1] for x in zip(medium_effect_theoretical_rg_files, range(540, 540+len(medium_effect_theoretical_rg_files)))]

seeded_medium_effect_sumher_rg_files = [x[0] % x[1] for x in zip(medium_effect_sumher_rg_files, range(540, 540+len(medium_effect_sumher_rg_files)))]

seeded_medium_effect_gps_files = [x[0] % x[1] for x in zip(medium_effect_gps_files, range(540, 540+len(medium_effect_gps_files)))]

seeded_medium_effect_hoeffdings_files = [x[0] % x[1] for x in zip(medium_effect_hoeffdings_files, range(540, 540+len(medium_effect_hoeffdings_files)))]

small_effect_rg_estimate_free_h2_free_rg_files = [f"results/ldsc/rg/free_h2_free_rg_intercept/whole_genome/randomised/{size[0]}_{size[1]}_{size[2]}_{size[3]}/{effect_tuple}_seed_%d_{tag_pair}.log" for size in small_effect_sample_sizes for tag_pair in tag_pairs for effect_tuple in small_effect_tuples]

small_effect_rg_estimate_fixed_h2_free_rg_files = [f"results/ldsc/rg/fixed_h2_free_rg_intercept/whole_genome/randomised/{size[0]}_{size[1]}_{size[2]}_{size[3]}/{effect_tuple}_seed_%d_{tag_pair}.log" for size in small_effect_sample_sizes for tag_pair in tag_pairs for effect_tuple in small_effect_tuples]

small_effect_rg_estimate_fixed_h2_fixed_rg_files = [f"results/ldsc/rg/fixed_h2_fixed_rg_intercept/whole_genome/randomised/{size[0]}_{size[1]}_{size[2]}_{size[3]}/{effect_tuple}_seed_%d_{tag_pair}.log" for size in small_effect_sample_sizes for tag_pair in tag_pairs for effect_tuple in small_effect_tuples]

small_effect_rg_estimate_free_h2_fixed_rg_files = [f"results/ldsc/rg/free_h2_fixed_rg_intercept/whole_genome/randomised/{size[0]}_{size[1]}_{size[2]}_{size[3]}/{effect_tuple}_seed_%d_{tag_pair}.log" for size in small_effect_sample_sizes for tag_pair in tag_pairs for effect_tuple in small_effect_tuples]

small_effect_theoretical_rg_files = [f"results/ldsc/rg/whole_genome/randomised/theoretical_rg/{size[0]}_{size[1]}_{size[2]}_{size[3]}/{effect_tuple}_seed_%d_{tag_pair}_theo_rg.tsv" for size in small_effect_sample_sizes for tag_pair in tag_pairs for effect_tuple in small_effect_tuples]

small_effect_sumher_rg_files = [f"results/ldak/ldak-thin/1000g/rg/{size[0]}_{size[1]}_{size[2]}_{size[3]}/{effect_tuple}_seed_%d_{tag_pair}.progress" for size in small_effect_sample_sizes for tag_pair in tag_pairs for effect_tuple in small_effect_tuples]

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

mixed_effect_rg_estimate_fixed_h2_free_rg_files = [f"results/ldsc/rg/fixed_h2_free_rg_intercept/whole_genome/randomised/{size[0]}_{size[1]}_{size[2]}_{size[3]}/{effect_tuple}_seed_%d_{tag_pair}.log" for size in mixed_effect_sample_sizes for tag_pair in tag_pairs for effect_tuple in mixed_effect_tuples]

mixed_effect_theoretical_rg_files = [f"results/ldsc/rg/whole_genome/randomised/theoretical_rg/{size[0]}_{size[1]}_{size[2]}_{size[3]}/{effect_tuple}_seed_%d_{tag_pair}_theo_rg.tsv" for size in mixed_effect_sample_sizes for tag_pair in tag_pairs for effect_tuple in mixed_effect_tuples]

mixed_effect_sumher_rg_files = [f"results/ldak/ldak-thin/1000g/rg/{size[0]}_{size[1]}_{size[2]}_{size[3]}/{effect_tuple}_seed_%d_{tag_pair}.progress" for size in mixed_effect_sample_sizes for tag_pair in tag_pairs for effect_tuple in mixed_effect_tuples]

mixed_effect_gps_files = [f"results/gps/simgwas/randomised/window_1000kb_step_50/{size[0]}_{size[1]}_{size[2]}_{size[3]}/3000_permutations/{effect_tuple}_seed_%d_tags_{tag_pair}_gps_pvalue.tsv" for size in mixed_effect_sample_sizes for tag_pair in tag_pairs for effect_tuple in mixed_effect_tuples]

mixed_effect_hoeffdings_files = [f"results/hoeffdings/simgwas/randomised/window_1000kb_step_50/{size[0]}_{size[1]}_{size[2]}_{size[3]}/{effect_tuple}_seed_%d_tags_{tag_pair}_hoeffdings.tsv" for size in mixed_effect_sample_sizes for tag_pair in tag_pairs for effect_tuple in mixed_effect_tuples]

seeded_mixed_effect_rg_estimate_fixed_h2_free_rg_files = [x[0] % x[1] for x in zip(mixed_effect_rg_estimate_fixed_h2_free_rg_files, range(9000, 9000+len(mixed_effect_rg_estimate_fixed_h2_free_rg_files)))]

seeded_mixed_effect_theoretical_rg_files = [x[0] % x[1] for x in zip(mixed_effect_theoretical_rg_files, range(9000, 9000+len(mixed_effect_theoretical_rg_files)))]

seeded_mixed_effect_sumher_rg_files = [x[0] % x[1] for x in zip(mixed_effect_sumher_rg_files, range(9000, 9000+len(mixed_effect_sumher_rg_files)))]

seeded_mixed_effect_gps_files = [x[0] % x[1] for x in zip(mixed_effect_gps_files, range(9000, 9000+len(mixed_effect_gps_files)))]

seeded_mixed_effect_hoeffdings_files = [x[0] % x[1] for x in zip(mixed_effect_hoeffdings_files, range(9000, 9000+len(mixed_effect_hoeffdings_files)))]

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
        seeded_medium_effect_gps_files,
        seeded_medium_effect_hoeffdings_files,
        seeded_medium_effect_sumher_rg_files

rule run_mixed_effect_jobs:
    input:
        seeded_mixed_effect_rg_estimate_fixed_h2_free_rg_files,
        seeded_mixed_effect_gps_files,
        seeded_mixed_effect_hoeffdings_files,
        seeded_mixed_effect_sumher_rg_files

rule compile_small_rg_results:
    input:
        seeded_small_effect_rg_estimate_free_h2_free_rg_files+
        seeded_small_effect_rg_estimate_fixed_h2_free_rg_files+
        seeded_small_effect_rg_estimate_fixed_h2_fixed_rg_files+
        seeded_small_effect_rg_estimate_free_h2_fixed_rg_files
    output:
        "results/ldsc/rg/whole_genome/randomised/compiled_s_rg_results.tsv"
    run:
        compile_rg_results(input, output)

rule compile_medium_rg_results:
    input:
        seeded_medium_effect_rg_estimate_free_h2_free_rg_files+
        seeded_medium_effect_rg_estimate_fixed_h2_free_rg_files+
        seeded_medium_effect_rg_estimate_fixed_h2_fixed_rg_files+
        seeded_medium_effect_rg_estimate_free_h2_fixed_rg_files
    output:
        "results/ldsc/rg/whole_genome/randomised/compiled_m_rg_results.tsv"
    run:
        compile_rg_results(input, output)

rule compile_mixed_rg_results:
    input:
        seeded_mixed_effect_rg_estimate_fixed_h2_free_rg_files
    output:
        "results/ldsc/rg/whole_genome/randomised/compiled_mixed_rg_results.tsv"
    run:
        compile_rg_results(input, output)

rule compile_small_effect_theoretical_rg:
    input:
        seeded_small_effect_theoretical_rg_files
    output:
         "results/ldsc/rg/whole_genome/randomised/theoretical_rg/compiled_s_theoretical_rg.tsv"
    run:
        compile_theoretical_rg_results(input, output)

rule compile_medium_effect_theoretical_rg:
    input:
        seeded_medium_effect_theoretical_rg_files
    output:
        "results/ldsc/rg/whole_genome/randomised/theoretical_rg/compiled_m_theoretical_rg.tsv"
    run:
        compile_theoretical_rg_results(input, output)

rule compile_mixed_effect_theoretical_rg:
    input:
        seeded_mixed_effect_theoretical_rg_files
    output:
        "results/ldsc/rg/whole_genome/randomised/theoretical_rg/compiled_mixed_theoretical_rg.tsv"
    run:
        compile_theoretical_rg_results(input, output)

rule compile_small_effect_sumher_results:
    input:
        [x.replace('progress', 'log') for x in seeded_small_effect_sumher_rg_files]
    output:
        compiled_rg_file = "results/ldak/ldak-thin/1000g/rg/compiled_s_sumher_rg.tsv"
    run:
        compile_sumher_results_from_log_files(input, output)

rule compile_medium_effect_sumher_results:
    input:
        [x.replace('progress', 'log') for x in seeded_medium_effect_sumher_rg_files]
    output:
        compiled_rg_file = "results/ldak/ldak-thin/1000g/rg/compiled_m_sumher_rg.tsv"
    run:
        compile_sumher_results_from_log_files(input, output)

rule compile_mixed_effect_sumher_results:
    input:
        [x.replace('progress', 'log') for x in seeded_mixed_effect_sumher_rg_files]
    output:
        compiled_rg_file = "results/ldak/ldak-thin/1000g/rg/compiled_mixed_sumher_rg.tsv"
    run:
        compile_sumher_results_from_log_files(input, output)

rule compile_small_gps_results:
    input:
        seeded_small_effect_gps_files
    output:
        "results/gps/simgwas/randomised/window_1000kb_step_50/compiled_s_gps_results.tsv"
    run:
        compile_gps_results(input, output)

rule compile_medium_gps_results:
    input:
        seeded_medium_effect_gps_files
    output:
        "results/gps/simgwas/randomised/window_1000kb_step_50/compiled_m_gps_results.tsv"
    run:
        compile_gps_results(input, output)

rule compile_mixed_gps_results:
    input:
        seeded_mixed_effect_gps_files
    output:
        "results/gps/simgwas/randomised/window_1000kb_step_50/compiled_mixed_gps_results.tsv"
    run:
        compile_gps_results(input, output)

rule compile_small_hoeffdings_results:
    input:
        seeded_small_effect_hoeffdings_files
    output:
        "results/hoeffdings/simgwas/randomised/window_1000kb_step_50/compiled_s_hoeffdings_results.tsv"
    run:
        compile_hoeffdings_results(input, output)

rule compile_medium_hoeffdings_results:
    input:
        seeded_medium_effect_hoeffdings_files
    output:
        "results/hoeffdings/simgwas/randomised/window_1000kb_step_50/compiled_m_hoeffdings_results.tsv"
    run:
        compile_hoeffdings_results(input, output)

rule compile_mixed_hoeffdings_results:
    input:
        seeded_mixed_effect_hoeffdings_files
    output:
        "results/hoeffdings/simgwas/randomised/window_1000kb_step_50/compiled_mixed_hoeffdings_results.tsv"
    run:
        compile_hoeffdings_results(input, output)

rule compile_results:
    input:
        "results/ldsc/rg/whole_genome/randomised/compiled_s_rg_results.tsv",
        "results/ldsc/rg/whole_genome/randomised/compiled_m_rg_results.tsv",
        "results/ldsc/rg/whole_genome/randomised/compiled_mixed_rg_results.tsv",
        "results/ldsc/rg/whole_genome/randomised/theoretical_rg/compiled_s_theoretical_rg.tsv",
        "results/ldsc/rg/whole_genome/randomised/theoretical_rg/compiled_m_theoretical_rg.tsv",
        "results/ldsc/rg/whole_genome/randomised/theoretical_rg/compiled_mixed_theoretical_rg.tsv",
        "results/ldak/ldak-thin/1000g/rg/compiled_s_sumher_rg.tsv",
        "results/ldak/ldak-thin/1000g/rg/compiled_m_sumher_rg.tsv",
        "results/ldak/ldak-thin/1000g/rg/compiled_mixed_sumher_rg.tsv",
        "results/gps/simgwas/randomised/window_1000kb_step_50/compiled_s_gps_results.tsv",
        "results/gps/simgwas/randomised/window_1000kb_step_50/compiled_m_gps_results.tsv",
        "results/gps/simgwas/randomised/window_1000kb_step_50/compiled_mixed_gps_results.tsv",
        "results/hoeffdings/simgwas/randomised/window_1000kb_step_50/compiled_s_hoeffdings_results.tsv",
        "results/hoeffdings/simgwas/randomised/window_1000kb_step_50/compiled_m_hoeffdings_results.tsv",
        "results/hoeffdings/simgwas/randomised/window_1000kb_step_50/compiled_mixed_hoeffdings_results.tsv"
"""
