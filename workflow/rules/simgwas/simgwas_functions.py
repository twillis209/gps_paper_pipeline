def get_simulation_runtime(wildcards, attempt):
    return 15 + attempt*5

def get_null_block_files(wildcards):
    effect_block_files = get_effect_block_files(wildcards)

    null_block_files_to_omit = [re.sub('small|medium|large|vlarge|huge|intermediate|random_\d+-\d+', 'null', x) for x in effect_block_files]

    null_block_files = [f"results/simgwas/simulated_sum_stats/block_sum_stats/{{no_reps}}_reps/null/{{ncases}}_{{ncontrols}}/chr{row.chr}/block_{row.block}_seed_{row.null_seed}_sum_stats.tsv.gz" for row in block_daf.itertuples()]

    for x in null_block_files_to_omit:
        if x in null_block_files:
            null_block_files.remove(x)

    return null_block_files

def get_effect_block_files(wildcards):
    if wildcards.effect_blocks == 'null':
        return []

    block_file_format = "results/simgwas/simulated_sum_stats/block_sum_stats/{no_reps}_reps/%s/{ncases}_{ncontrols}/chr%d/block_%d_seed_%d_sum_stats.tsv.gz"

    effect_block_files = []

    for x in wildcards.effect_blocks.split('/')[-1].split('+'):
        block_match = re.match('^(\d+)-(.+)', x)

        chrom = int(block_match.group(1))

        if ':' in block_match.group(2):
            range_match = re.match('([smlvhi]|r\d+-\d+-)(\d+):(\d+)', block_match.group(2))

            # random effect
            # TODO no handling of seed in this atm
            if 'r' in range_match.group(1):
                effect = 'random_' + range_match.group(1)[1:-1]

                effect_block_files += [block_file_format % (chrom, effect, y) for y in range(int(range_match.group(2)), int(range_match.group(3))+1) if y in block_dict[chrom]]
            # non-random effect
            else:
                effect = effect_size_dict[range_match.group(1)]

                seed_col = effect + "_seed"

                for block in range(int(range_match.group(2)), int(range_match.group(3))+1):
                    if block in block_daf.query('chr == @chrom')['block']:

                        seed = block_daf.query('chr == @chrom & block == @block')[seed_col].values[0]

                        effect_block_files.append(block_file_format % (effect, chrom, block, seed))

        else:
            singleton_match = re.match('([smlvhi]|r\d+-\d+-)(\d+)', x)

            block = int(singleton_match.group(2))

            # random effect
            # TODO no handling of seed in this atm
            if 'r' in x:
                effect = 'random_' + singleton_match.group(1)[1:-1]

                if int(singleton_match.group(2)) in block_dict[chrom]:
                    effect_block_files += [block_file_format % (chrom, effect, int(singleton_match.group(2)))]
            # non-random effect
            else:
                effect = effect_size_dict[singleton_match.group(1)]

                seed_col = effect + "_seed"

                seed = block_daf.query('chr == @chrom & block == @block')[seed_col].values[0]

                if block in block_daf.query('chr == @chrom')['block']:
                    effect_block_files.append(block_file_format % (effect, chrom, block, seed))

    return effect_block_files
