library(data.table)
setDTthreads(8)



sum_stats_dat <- fread(snakemake@input[['sum_stats_file']], sep = '\t', header = T)

for(i in 1:22) {
  bim_dat <- fread(snakemake@input[['bim_files']][i], sep = '\t', header = F, col.names = c('chr', 'ID', 'Cm', 'bp', 'A1', 'A2'))

  bim_dat[, variant_12 := paste(chr, bp, A1, A2, sep = ':')]
  bim_dat[, variant_21 := paste(chr, bp, A2, A1, sep = ':')]

  bim_dat <- bim_dat[variant_12 %in% sum_stats_dat$variant | variant_21 %in% sum_stats_dat$variant]

  fwrite(bim_dat, file = snakemake@output[['bim_files']][i], row.names = F, sep = ' ', col.names = F, quote = F)
}

bim_dat <- fread(snakemake@input[['bim_files']][23], sep = '\t', header = F, col.names = c('chr', 'ID', 'Cm', 'bp', 'A1', 'A2'))

bim_dat[, variant_12 := paste(chr, bp, A1, A2, sep = ':')]
bim_dat[, variant_21 := paste(chr, bp, A2, A1, sep = ':')]

bim_dat <- bim_dat[variant_12 %in% sum_stats_dat$variant | variant_21 %in% sum_stats_dat$variant]

fwrite(bim_dat, file = snakemake@output[['bim_files']][23], row.names = F, sep = ' ', col.names = F, quote = F)
