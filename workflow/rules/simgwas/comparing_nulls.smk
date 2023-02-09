import pandas as pd

rule comparing_gps_null_distributions:
    input:
        perm_file = "results/gps/simgwas/400_reps/randomised/10000_10000_10000_10000/s400_s400_s0/window_1000kb_step_50_r2_0_2/3000_permutations/seed_48801_tags_1-2.tsv",
        gps_files = [x[1] % (int(x[0])+48801) for x in enumerate([f"results/gps/simgwas/400_reps/randomised/10000_10000_10000_10000/s400_s400_s0/window_1000kb_step_50_r2_0_2/3000_permutations/seed_%d_tags_{x}-{x+1}.tsv" for x in range(1, 400, 2)])]
    output:
        "results/gps/simgwas/400_reps/randomised/10000_10000_10000_10000/s400_s400_s0/window_1000kb_step_50_r2_0_2/3000_permutations/comparing_perm_and_real_null.tsv"
    run:
        perm_gps = []

        sim_gps = []

        with open(input.perm_file, 'r') as perm_file:
            line = perm_file.readline()
            for i in range(len(input.gps_files)):
                line = perm_file.readline()
                perm_gps.append(float(line.strip()))

        for x in input.gps_files:
            with open(x, 'r') as gps_file:
                line = gps_file.readline()
                line = gps_file.readline()
                sim_gps.append(float(line.split('\t')[0]))

        pd.DataFrame({'sim_gps': sim_gps, 'perm_gps': perm_gps}).to_csv(output[0], sep = '\t', index = False)
