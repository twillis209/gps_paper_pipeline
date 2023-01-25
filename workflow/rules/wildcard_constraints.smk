effect_blocks_regex = "(t|s|m|i)\d+(-(t|s|m|i)\d+)?"

wildcard_constraints:
    chr = "chr[0-9XY]{1,2}",
    no_reps = "\d+",
    ncases_A = "\d+",
    ncases_B = "\d+",
    ncontrols_A = "\d+",
    ncontrols_B = "\d+",
    tag_A = "\d+",
    tag_B = "\d+",
    pert = "\d+",
    ecdf = "naive|pp|pert",
    window = "\d+kb",
    step = "\d+",
    seed = "\d+",
    block = "\d+",
    effect_blocks_A = effect_blocks_regex,
    effect_blocks_B  = effect_blocks_regex,
    shared_effect_blocks  = effect_blocks_regex,
    effect_blocks = effect_blocks_regex,
    join = 'ukbb',
    snp_set = 'all_pruned_snps|sans_mhc|with_mhc'

