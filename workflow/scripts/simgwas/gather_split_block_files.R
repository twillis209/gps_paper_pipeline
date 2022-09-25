library(data.table)
setDTthreads(snakemake@threads)
library(magrittr)

z_column_name_A <- sprintf("zsim.%s", snakemake@wildcards[['tag_A']])
beta_column_name_A <- sprintf("betasim.%s", snakemake@wildcards[['tag_A']])
p_column_name_A <- sprintf("p.%s", snakemake@wildcards[['tag_A']])

header_A <- c("position", "a0", "a1", "id", "block", "TYPE", "EUR", z_column_name_A, beta_column_name_A, p_column_name_A, "ncases", "ncontrols", "chr", "block_effect_size")

snakemake@input[['a_files']] %>%
  lapply(., fread, sep = '\t', header = F) %>%
  rbindlist -> a_dat

names(a_dat) <- header_A

fwrite(a_dat, sep = '\t', file = snakemake@output[['combined_sum_stats_A']], col.names = T)

rm(a_dat)

z_column_name_B <- sprintf("zsim.%s", snakemake@wildcards[['tag_B']])
beta_column_name_B <- sprintf("betasim.%s", snakemake@wildcards[['tag_B']])
p_column_name_B <- sprintf("p.%s", snakemake@wildcards[['tag_B']])

header_B <- c("position", "a0", "a1", "id", "block", "TYPE", "EUR", z_column_name_B, beta_column_name_B, p_column_name_B, "ncases", "ncontrols", "chr", "block_effect_size")

snakemake@input[['b_files']] %>%
  lapply(., fread, sep = '\t', header = F) %>%
  rbindlist -> b_dat

names(b_dat) <- header_B

fwrite(b_dat, sep = '\t', file = snakemake@output[['combined_sum_stats_B']], col.names = T)
