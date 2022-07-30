import pandas as pd
from collections import namedtuple

EffectBlocks = namedtuple("EffectBlocks", "a_blocks b_blocks shared_blocks")

small_effect_tuples = [EffectBlocks(a_blocks = "s400", b_blocks = "s400", shared_blocks = f"s{x}") for x in range(0, 401, 50)]

medium_effect_tuples = [EffectBlocks(a_blocks = "m25", b_blocks = "m25", shared_blocks = f"m{x}") for x in [0, 2, 5, 7, 10, 12, 15, 17, 20, 25]]+[EffectBlocks(a_blocks = "m50", b_blocks = "m50", shared_blocks = f"m{x}") for x in [0, 2, 5, 7, 10, 12, 15, 17, 20, 30, 40, 50]]

mixed_effect_tuples = [EffectBlocks(a_blocks = "s200-m25", b_blocks = "s200-m25", shared_blocks = f"s{x}-m{y}") for x in range(0, 201, 50) for y in range(0, 26, 5)]

effect_tuples = small_effect_tuples + medium_effect_tuples + mixed_effect_tuples

tags = [str(x) for x in range(1, 401)]

sim_dicts = []

seed = 1

#ncases_A ncontrols_A ncases_B ncontrols_B seed
for sample_tuple in snakemake.params.sample_sizes:
    for effect_tuple in effect_tuples:
        for i in range(len(tags))[::2]:
            sim_dicts.append(
                            {
                                "a_blocks" : effect_tuple.a_blocks,
                                "b_blocks" : effect_tuple.b_blocks,
                                "shared_blocks" : effect_tuple.shared_blocks,
                                "ncases_A" : sample_tuple[0],
                                "ncontrols_A" : sample_tuple[1],
                                "ncases_B" : sample_tuple[0],
                                "ncontrols_B" : sample_tuple[1],
                                "tag_A" : tags[i],
                                "tag_B" : tags[i+1],
                                "seed" : seed
                            }
                        )
            seed += 1

daf = pd.DataFrame(sim_dicts)

daf = daf[["a_blocks", "b_blocks", "shared_blocks", "ncases_A", "ncontrols_A", "ncases_B", "ncontrols_B", "tag_A", "tag_B", "seed"]]

daf.to_csv(snakemake.output[0], sep = '\t', index = False)
