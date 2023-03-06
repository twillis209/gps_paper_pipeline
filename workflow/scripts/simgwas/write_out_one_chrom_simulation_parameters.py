import pandas as pd
from collections import namedtuple

EffectBlocks = namedtuple("EffectBlocks", "a_blocks b_blocks shared_blocks")

medium_effect_tuples = [EffectBlocks(a_blocks = "m20", b_blocks = "m20", shared_blocks = f"m{x}") for x in range(0, 21, 5)]

effect_tuples = medium_effect_tuples

tags = [str(x) for x in range(1, 401)]

sim_dicts = []

seed = 1

sample_sizes = [(500, 10000),
                (1000, 10000),
                (5000, 10000),
                (10000, 10000),
                (100000, 100000)
                ]

r2_values = [0.2, 0.5, 0.8]

for sample_tuple in sample_sizes:
    for effect_tuple in effect_tuples:
        for r2 in r2_values:
            for i in range(len(tags))[::2]:
                sim_dicts.append(
                                {
                                    "chr" : "chr1",
                                    "a_blocks" : effect_tuple.a_blocks,
                                    "b_blocks" : effect_tuple.b_blocks,
                                    "shared_blocks" : effect_tuple.shared_blocks,
                                    "ncases_A" : sample_tuple[0],
                                    "ncontrols_A" : sample_tuple[1],
                                    "ncases_B" : sample_tuple[0],
                                    "ncontrols_B" : sample_tuple[1],
                                    "tag_A" : tags[i],
                                    "tag_B" : tags[i+1],
                                    "r2" : r2,
                                    "window" : "1000kb",
                                    "step" : 50,
                                    "seed" : seed
                                }
                            )
                seed += 1

daf = pd.DataFrame(sim_dicts)

daf = daf[["chr", "a_blocks", "b_blocks", "shared_blocks", "ncases_A", "ncontrols_A", "ncases_B", "ncontrols_B", "tag_A", "tag_B", "r2", "window", "step", "seed"]]

daf.to_csv('results/simgwas/one_chrom_simulation_parameters.tsv', sep = '\t', index = False)
