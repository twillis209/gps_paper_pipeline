from snakemake.shell import shell

def read_block_files(input_file):
    with open(input_file, "r") as fh:
        lines = fh.readlines()

    return [x.strip('\n') for x in lines]

z_column_index_A = 8 + int(snakemake.wildcards.tag_A)
beta_column_index_A = 8 + (int(snakemake.wildcards.no_reps) * 2) + int(snakemake.wildcards.tag_A)
p_column_index_A = 8 + (int(snakemake.wildcards.no_reps) * 3) + int(snakemake.wildcards.tag_A)

ncases_column_index = 8  +  (int(snakemake.wildcards.no_reps)  *  4) + 2
ncontrols_column_index = 8 + (int(snakemake.wildcards.no_reps) * 4) + 3
chr_column_index = 8 + (int(snakemake.wildcards.no_reps) * 4) + 4
block_effect_column_index = 8 + (int(snakemake.wildcards.no_reps) * 4) + 5

cut_string_A = f"1-7,{z_column_index_A},{beta_column_index_A},{p_column_index_A},{ncases_column_index},{ncontrols_column_index},{chr_column_index},{block_effect_column_index}"

a_block_files = read_block_files(snakemake.input.a_block_file)

for x in a_block_files:
    shell("cut -f{cut_string_A} {x} >> {snakemake.output.combined_sum_stats_A}")

line_count_A = shell("wc -l {snakemake.output.combined_sum_stats_A}", read = True).split()[0]

if int(line_count_A) < snakemake.params.line_count_bound:
    raise Exception(f"Too many lines in A file")

z_column_index_B = 8 + int(snakemake.wildcards.tag_B)
beta_column_index_B = 8 + (int(snakemake.wildcards.no_reps) * 2) + int(snakemake.wildcards.tag_B)
p_column_index_B = 8 + (int(snakemake.wildcards.no_reps) * 3) + int(snakemake.wildcards.tag_B)

cut_string_B = f"1-7,{z_column_index_B},{beta_column_index_B},{p_column_index_B},{ncases_column_index},{ncontrols_column_index},{chr_column_index},{block_effect_column_index}"

b_block_files = read_block_files(snakemake.input.b_block_file)

for x in b_block_files:
    shell("cut -f{cut_string_B} {x} >> {snakemake.output.combined_sum_stats_B}")

line_count_B = shell("wc -l {snakemake.output.combined_sum_stats_B}", read = True).split()[0]

if int(line_count_B) < snakemake.params.line_count_bound:
    raise Exception(f"Too many lines in B file")

