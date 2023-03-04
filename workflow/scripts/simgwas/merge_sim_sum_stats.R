library(data.table)

setDTthreads(snakemake@threads)

sum_stats_A_dat <- fread(snakemake@input[['sum_stats_file_A']], sep = '\t', header = T, select = c('id', 'chr', 'position', 'a0', 'a1', 'block', snakemake@params[['file_A_stat_cols']], 'block_effect_size'))
sum_stats_B_dat <- fread(snakemake@input[['sum_stats_file_B']], sep = '\t', header = T, select = c('id', 'chr', 'position', 'a0', 'a1', 'block', snakemake@params[['file_B_stat_cols']], 'block_effect_size'))

merged_dat <- merge(sum_stats_A_dat, sum_stats_B_dat, by = c('chr', 'position', 'a0', 'a1', 'block'), suffixes = c(".A", ".B"))

merged_dat[, id.B := NULL]

setnames(merged_dat, 'id.A', 'id')

if(is.null(snakemake@wildcards[['chr']]) & length(unique(merged_dat$chr)) != 22) {
  stop(sprintf("Incorrect number of chromosomes present in merged summary statistics: %s", length(unique(merged_dat$chr))))
} else {
  merged_dat <- unique(merged_dat, by = c('chr', 'position'))
  fwrite(merged_dat, file = snakemake@output[[1]], sep = '\t')
}
