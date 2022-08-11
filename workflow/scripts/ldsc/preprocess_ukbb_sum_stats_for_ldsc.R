library(data.table)

setDTthreads(snakemake@resources[['threads']])

sum_stats_dat <- fread(snakemake@input[['sum_stats_file']], sep = '\t', header = T, tmpdir = snakemake@resources[['tmpdir']], select = c('variant', snakemake@params[['pval_col']], snakemake@params[['tstat_col']], snakemake@params[['n_col']]))

sum_stats_dat[, c('chr', 'bp', 'ref', 'alt') := tstrsplit(variant, split = ':')]

sum_stats_dat[, chr := as.character(chr)]
sum_stats_dat[, bp := as.character(bp)]

snp_dat <- fread(snakemake@input[['snplist_file']], sep = '\t', header = T)

snp_dat[, CHR := as.character(CHR)]
snp_dat[, BP := as.character(BP)]

merged_dat <- merge(sum_stats_dat, snp_dat, by.x = c('chr', 'bp'), by.y = c('CHR', 'BP'))

merged_dat <- merged_dat[(ref == A1  & alt == A2) | (ref == A2 & alt == A1)]

merged_dat[, c('ref', 'alt') := NULL]

setnames(merged_dat, c('chr', 'bp'), c('CHR', 'BP'))

fwrite(merged_dat, file = snakemake@output[[1]], sep = '\t')
