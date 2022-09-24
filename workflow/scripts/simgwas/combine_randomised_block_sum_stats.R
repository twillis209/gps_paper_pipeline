library(data.table)
setDTthreads(snakemake@threads)
library(magrittr)

a_block_file <- snakemake@input[['a_block_file']]
b_block_file <- snakemake@input[['b_block_file']]

a_block_files <- scan(a_block_file, what = character())
b_block_files <- scan(b_block_file, what = character())

ncases_column_index <- 8  +  (as.integer(snakemake@wildcards[['no_reps']])  *  4) + 2
ncontrols_column_index <- 8 + (as.integer(snakemake@wildcards[['no_reps']]) * 4) + 3
chr_column_index <- 8 + (as.integer(snakemake@wildcards[['no_reps']]) * 4) + 4
block_effect_column_index <- 8 + (as.integer(snakemake@wildcards[['no_reps']]) * 4) + 5

z_column_index_A <- 8 + as.integer(snakemake@wildcards[['tag_A']])
beta_column_index_A <- 8 + (as.integer(snakemake@wildcards[['no_reps']]) * 2) + as.integer(snakemake@wildcards[['tag_A']])
p_column_index_A <- 8 + (as.integer(snakemake@wildcards[['no_reps']]) * 3) + as.integer(snakemake@wildcards[['tag_A']])

a_cols <- c(1:7, z_column_index_A, beta_column_index_A, p_column_index_A, ncases_column_index, ncontrols_column_index, chr_column_index, block_effect_column_index)

a_block_files %>%
  lapply(., fread, sep = '\t', header = F, select = a_cols) %>%
  rbindlist -> a_block_dat

fwrite(a_block_dat, file = snakemake@output[['combined_sum_stats_A']], sep = '\t', col.names = F)

rm(a_block_dat)

z_column_index_B <- 8 + as.integer(snakemake@wildcards[['tag_B']])
beta_column_index_B <- 8 + (as.integer(snakemake@wildcards[['no_reps']]) * 2) + as.integer(snakemake@wildcards[['tag_B']])
p_column_index_B <- 8 + (as.integer(snakemake@wildcards[['no_reps']]) * 3) + as.integer(snakemake@wildcards[['tag_B']])

b_cols <- c(1:7, z_column_index_B, beta_column_index_B, p_column_index_B, ncases_column_index, ncontrols_column_index, chr_column_index, block_effect_column_index)

b_block_files %>%
  lapply(., fread, sep = '\t', header = F, select = b_cols) %>%
  rbindlist -> b_block_dat

fwrite(b_block_dat, file = snakemake@output[['combined_sum_stats_B']], sep = '\t', col.names = F)
