library(data.table)
library(parallel)
library(magrittr)

setDTthreads(1)

ncases_column_index <- 8 + (as.integer(snakemake@wildcards[['no_reps']])  *  4) + 2
ncontrols_column_index <- 8 + (as.integer(snakemake@wildcards[['no_reps']]) * 4) + 3
chr_column_index <- 8 + (as.integer(snakemake@wildcards[['no_reps']]) * 4) + 4
block_effect_column_index <- 8 + (as.integer(snakemake@wildcards[['no_reps']]) * 4) + 5

z_column_name_A = sprintf("zsim.%s", snakemake@wildcards[['tag_A']])
beta_column_name_A = sprintf("betasim.%s", snakemake@wildcards[['tag_A']])
p_column_name_A = sprintf("p.%s", snakemake@wildcards[['tag_A']])

z_column_index_A <- 8 + as.integer(snakemake@wildcards[['tag_A']])
beta_column_index_A <- 8 + (as.integer(snakemake@wildcards[['no_reps']]) * 2) + as.integer(snakemake@wildcards[['tag_A']])
p_column_index_A <- 8 + (as.integer(snakemake@wildcards[['no_reps']]) * 3) + as.integer(snakemake@wildcards[['tag_A']])

a_col_indices <- c(1:7, z_column_index_A, beta_column_index_A, p_column_index_A, ncases_column_index, ncontrols_column_index, chr_column_index, block_effect_column_index)

a_cols <- c("position", "a0", "a1", "id", "block", "TYPE", "EUR", z_column_name_A, beta_column_name_A, p_column_name_A, "ncases", "ncontrols", "chr", "block_effect_size")

mclapply(snakemake@input[['a_block_files']], fread, sep = '\t', header = F, select = a_col_indices, mc.cores = snakemake@threads) %>% rbindlist -> a_dat

names(a_dat) <- a_cols

setDTthreads(snakemake@threads)

fwrite(a_dat, file = snakemake@output[['combined_sum_stats_A']], sep = '\t')

setDTthreads(1)

z_column_name_B = sprintf("zsim.%s", snakemake@wildcards[['tag_B']])
beta_column_name_B = sprintf("betasim.%s", snakemake@wildcards[['tag_B']])
p_column_name_B = sprintf("p.%s", snakemake@wildcards[['tag_B']])

z_column_index_B <- 8 + as.integer(snakemake@wildcards[['tag_B']])
beta_column_index_B <- 8 + (as.integer(snakemake@wildcards[['no_reps']]) * 2) + as.integer(snakemake@wildcards[['tag_B']])
p_column_index_B <- 8 + (as.integer(snakemake@wildcards[['no_reps']]) * 3) + as.integer(snakemake@wildcards[['tag_B']])

b_col_indices <- c(1:7, z_column_index_B, beta_column_index_B, p_column_index_B, ncases_column_index, ncontrols_column_index, chr_column_index, block_effect_column_index)

b_cols <- c("position", "a0", "a1", "id", "block", "TYPE", "EUR", z_column_name_B, beta_column_name_B, p_column_name_B, "ncases", "ncontrols", "chr", "block_effect_size")

mclapply(snakemake@input[['b_block_files']], fread, sep = '\t', header = F, select = b_col_indices, mc.cores = snakemake@threads) %>% rbindlist -> b_dat

names(b_dat) <- b_cols

setDTthreads(snakemake@threads)

fwrite(b_dat, file = snakemake@output[['combined_sum_stats_B']], sep = '\t')
