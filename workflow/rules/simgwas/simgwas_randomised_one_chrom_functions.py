import pandas as pd

def get_randomised_chrom_block_tuples_for_chrom_for_pair(wildcards):
    random.seed(wildcards.seed)

    chrom = int(wildcards.chr.replace('chr', ''))
    shared_chrom_block_nos = []
    shared_chrom_block_dict = {}

    if wildcards.shared_effect_blocks != 'null':
        for token in wildcards.shared_effect_blocks.split('-'):
            block_match = re.match('([tsmlvhi])(\d+)', token)

            if not block_match:
                raise ValueError("Invalid block format: %s" % token)

            shared_effect = effect_size_dict[block_match.group(1)]

            no_of_shared_blocks = int(block_match.group(2))

            no_of_cvs_per_block = cv_per_block_dict[shared_effect]

            no_of_shared_blocks /= no_of_cvs_per_block

            i = 0

            shared_chrom_block_dict[shared_effect] = []

            while i < max(no_of_shared_blocks, 0):
                block_no = random.choice(block_daf.query('chr == @chrom')[['block']].values.transpose().tolist()[0])

                if (chrom, block_no) not in shared_chrom_block_nos:
                    shared_chrom_block_nos.append((chrom, block_no))
                    shared_chrom_block_dict[shared_effect].append((chrom, block_no, no_of_cvs_per_block))
                    i += 1

    a_chrom_block_nos = []
    a_chrom_block_dict = {}

    if wildcards.effect_blocks_A != 'null':
        for token in wildcards.effect_blocks_A.split('-'):
            block_match_a = re.match('([tsmlvhi])(\d+)', token)

            if not block_match_a:
                raise ValueError("Invalid block format: %s" % token)

            effect_a = effect_size_dict[block_match_a.group(1)]

            no_of_blocks_a = int(block_match_a.group(2))

            no_of_cvs_per_block_a = cv_per_block_dict[effect_a]

            no_of_blocks_a /= no_of_cvs_per_block_a

            i = 0

            a_chrom_block_dict[effect_a] = []

            if effect_a not in shared_chrom_block_dict.keys():
                shared_chrom_block_dict[effect_a] = []

            while i < max(no_of_blocks_a-len(shared_chrom_block_dict[effect_a]), 0):
                block_no = random.choice(block_daf.query('chr == @chrom')[['block']].values.transpose().tolist()[0])

                if (chrom, block_no) not in shared_chrom_block_nos and (chrom, block_no) not in a_chrom_block_nos:
                    a_chrom_block_nos.append((chrom, block_no))
                    a_chrom_block_dict[effect_a].append((chrom, block_no, no_of_cvs_per_block_a))
                    i += 1

    b_chrom_block_nos = []
    b_chrom_block_dict = {}

    if wildcards.effect_blocks_B != 'null':
        for token in wildcards.effect_blocks_B.split('-'):
            block_match_b = re.match('([tismlvh])(\d+)', token)

            if not block_match_b:
                raise ValueError("Invalid block format: %s" % token)

            effect_b = effect_size_dict[block_match_b.group(1)]

            no_of_blocks_b = int(block_match_b.group(2))

            no_of_cvs_per_block_b = cv_per_block_dict[effect_b]

            no_of_blocks_b /= no_of_cvs_per_block_b

            i = 0

            b_chrom_block_dict[effect_b] = []

            if effect_b not in shared_chrom_block_dict.keys():
                shared_chrom_block_dict[effect_b] = []

            while i < max(no_of_blocks_b-len(shared_chrom_block_dict[effect_b]), 0):
                block_no = random.choice(block_daf.query('chr == @chrom')[['block']].values.transpose().tolist()[0])


                if (chrom, block_no) not in shared_chrom_block_nos and (chrom, block_no) not in a_chrom_block_nos and (chrom, block_no) not in b_chrom_block_nos:
                    b_chrom_block_nos.append((chrom, block_no))
                    b_chrom_block_dict[effect_b].append((chrom, block_no, no_of_cvs_per_block_b))
                    i += 1

    return (shared_chrom_block_nos, a_chrom_block_nos, b_chrom_block_nos, shared_chrom_block_dict, a_chrom_block_dict, b_chrom_block_dict)

def get_randomised_block_files_for_chrom_for_pair(wildcards):
    shared_chrom_block_nos, a_chrom_block_nos, b_chrom_block_nos, shared_chrom_block_dict, a_chrom_block_dict, b_chrom_block_dict = get_randomised_chrom_block_tuples_for_chrom_for_pair(wildcards)

    block_files = []

    a_block_files = []
    b_block_files = []

    chrom = int(wildcards.chr.replace('chr', ''))

    for k in shared_chrom_block_dict.keys():
        seed_label = k + '_seed'

        for v in shared_chrom_block_dict[k]:
            seed = block_daf.query('chr == @v[0] & block == @v[1]')[seed_label].values[0]

            a_file = f"results/simgwas/simulated_sum_stats/block_sum_stats/400_reps/{k}/{v[2]}_cv/{wildcards.ncases_A}_{wildcards.ncontrols_A}/chr{v[0]}/block_{v[1]}_seed_{seed}_sum_stats.tsv.gz"

            b_file = f"results/simgwas/simulated_sum_stats/block_sum_stats/400_reps/{k}/{v[2]}_cv/{wildcards.ncases_B}_{wildcards.ncontrols_B}/chr{v[0]}/block_{v[1]}_seed_{seed}_sum_stats.tsv.gz"

            block_files.append(a_file)

            a_block_files.append(a_file)
            b_block_files.append(b_file)

            if wildcards.ncases_A != wildcards.ncases_B or wildcards.ncontrols_A != wildcards.ncontrols_B:
                block_files.append(b_file)

    for k in a_chrom_block_dict.keys():
        seed_label = k + '_seed'

        for v in a_chrom_block_dict[k]:
            seed = block_daf.query('chr == @v[0] & block == @v[1]')[seed_label].values[0]

            a_file = f"results/simgwas/simulated_sum_stats/block_sum_stats/400_reps/{k}/{v[2]}_cv/{wildcards.ncases_A}_{wildcards.ncontrols_A}/chr{v[0]}/block_{v[1]}_seed_{seed}_sum_stats.tsv.gz"

            block_files.append(a_file)
            a_block_files.append(a_file)

    for k in b_chrom_block_dict.keys():
        seed_label = k + '_seed'

        for v in b_chrom_block_dict[k]:
            seed = block_daf.query('chr == @v[0] & block == @v[1]')[seed_label].values[0]

            b_file = f"results/simgwas/simulated_sum_stats/block_sum_stats/400_reps/{k}/{v[2]}_cv/{wildcards.ncases_B}_{wildcards.ncontrols_B}/chr{v[0]}/block_{v[1]}_seed_{seed}_sum_stats.tsv.gz"

            b_block_files.append(b_file)

            if b_file not in block_files:
                block_files.append(b_file)

    # Fill in null files
    for block in block_daf.query('chr == @chrom')['block']:
        seed = block_daf.query('chr == @chrom & block == @block')['null_seed'].values[0]

        if (chrom, block) not in a_chrom_block_nos and (chrom, block) not in shared_chrom_block_nos:
            a_file = f"results/simgwas/simulated_sum_stats/block_sum_stats/400_reps/null/0_cv/{wildcards.ncases_A}_{wildcards.ncontrols_A}/chr{chrom}/block_{block}_seed_{seed}_sum_stats.tsv.gz"
            block_files.append(a_file)
            a_block_files.append(a_file)

        if (chrom, block) not in b_chrom_block_nos and (chrom, block) not in shared_chrom_block_nos:
            b_file = f"results/simgwas/simulated_sum_stats/block_sum_stats/400_reps/null/0_cv/{wildcards.ncases_B}_{wildcards.ncontrols_B}/chr{chrom}/block_{block}_seed_{seed}_sum_stats.tsv.gz"
            b_block_files.append(b_file)

            if b_file not in block_files:
                block_files.append(b_file)

    all_files = block_files+a_block_files+b_block_files

    return all_files
