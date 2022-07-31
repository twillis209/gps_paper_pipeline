from snakemake.shell import shell

with open(snakemake.log.log, 'w') as logfile:
    z_column_name_A = f"zsim.{snakemake.wildcards.tag_A}"
    beta_column_name_A = f"betasim.{snakemake.wildcards.tag_A}"
    p_column_name_A = f"p.{snakemake.wildcards.tag_A}"

    header_string_A = "\t".join(["position", "a0", "a1", "id", "block", "TYPE", "EUR", z_column_name_A, beta_column_name_A, p_column_name_A, "ncases", "ncontrols", "chr", "block_effect_size"])

    z_column_index_A = 8 + int(snakemake.wildcards.tag_A)
    beta_column_index_A = 8 + (int(snakemake.wildcards.no_reps) * 2) + int(snakemake.wildcards.tag_A)
    p_column_index_A = 8 + (int(snakemake.wildcards.no_reps) * 3) + int(snakemake.wildcards.tag_A)

    ncases_column_index = 8  +  (int(snakemake.wildcards.no_reps)  *  4) + 2
    ncontrols_column_index = 8 + (int(snakemake.wildcards.no_reps) * 4) + 3
    chr_column_index = 8 + (int(snakemake.wildcards.no_reps) * 4) + 4
    block_effect_column_index = 8 + (int(snakemake.wildcards.no_reps) * 4) + 5

    cut_string_A = f"1-7,{z_column_index_A},{beta_column_index_A},{p_column_index_A},{ncases_column_index},{ncontrols_column_index},{chr_column_index},{block_effect_column_index}"

    shell("echo -e \"{header_string_A}\" > {snakemake.params.uncomp_sum_stats_A}")

    for x in snakemake.input.a_block_files:
        shell("zcat {x} | cut -f{cut_string_A} >> {snakemake.params.uncomp_sum_stats_A}")
        logfile.write(f"{x}\n")

    shell(f"gzip {snakemake.params.uncomp_sum_stats_A}")

    logfile.write("\n")

    z_column_name_B = f"zsim.{snakemake.wildcards.tag_B}"
    beta_column_name_B = f"betasim.{snakemake.wildcards.tag_B}"
    p_column_name_B = f"p.{snakemake.wildcards.tag_B}"

    header_string_B = "\t".join(["position", "a0", "a1", "id", "block", "TYPE", "EUR", z_column_name_B, beta_column_name_B, p_column_name_B, "ncases", "ncontrols", "chr", "block_effect_size"])

    z_column_index_B = 8 + int(snakemake.wildcards.tag_B)
    beta_column_index_B = 8 + (int(snakemake.wildcards.no_reps) * 2) + int(snakemake.wildcards.tag_B)
    p_column_index_B = 8 + (int(snakemake.wildcards.no_reps) * 3) + int(snakemake.wildcards.tag_B)

    cut_string_B = f"1-7,{z_column_index_B},{beta_column_index_B},{p_column_index_B},{ncases_column_index},{ncontrols_column_index},{chr_column_index},{block_effect_column_index}"

    shell("echo -e \"{header_string_B}\" > {snakemake.params.uncomp_sum_stats_B}")

    for x in snakemake.input.b_block_files:
        shell("zcat {x} | cut -f{cut_string_B} >> {snakemake.params.uncomp_sum_stats_B}")
        logfile.write(f"{x}\n")

    shell(f"gzip {snakemake.params.uncomp_sum_stats_B}")
