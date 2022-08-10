library(data.table)

setDTthreads(snakemake@resources[['threads']])

sum_stats_dat <- fread(snakemake@input[['sum_stats_file']], sep = '\t', header = T, tmpdir = snakemake@params[['tmpdir']], select = c('variant', snakemake@params))

sum_stats_dat[, c('chr', 'bp', 'ref', 'alt') := tstrsplit(variant, split = ':')]

snp_dat <- fread(snakemake@input[['snplist_file']], sep = '\t', header = T)

merged_dat <- merge(sum_stats_dat, snp_dat, by.x = c('chr', 'bp'), by.y = c('CHR', 'BP'))

merged_dat <- merged_dat[(ref == A1  & alt == A2) | (ref == A2 & alt == A1)]

merged_dat[, c('ref', 'alt') := NULL]

fwrite(sum_stats_dat, file = snakemake@output[[1]], sep = '\t')
