library(data.table)
setDTthreads(8)

sum_stats_dat <- fread("resources/ukbb_sum_stats/merged_sum_stats.tsv", sep = '\t', header = T)

for(i in 1:22) {
  bim_dat <- fread(sprintf("resources/1000g/euro/qc/chr%d_qc.bim", i), sep = '\t', header = F, col.names = c('chr', 'ID', 'Cm', 'bp', 'A1', 'A2'))

  bim_dat[, variant_12 := paste(chr, bp, A1, A2, sep = ':')]
  bim_dat[, variant_21 := paste(chr, bp, A2, A1, sep = ':')]

  bim_dat <- bim_dat[variant_12 %in% sum_stats_dat$variant | variant_21 %in% sum_stats_dat$variant]

  fwrite(bim_dat, file = sprintf("resources/ukbb_sum_stats/plink_ranges/chr%d.txt", i), row.names = F, sep = ' ', col.names = F, quote = F)
}

bim_dat <- fread("resources/1000g/euro/qc/chrX_qc.bim", sep = '\t', header = F, col.names = c('chr', 'ID', 'Cm', 'bp', 'A1', 'A2'))

bim_dat[, variant_12 := paste(chr, bp, A1, A2, sep = ':')]
bim_dat[, variant_21 := paste(chr, bp, A2, A1, sep = ':')]

bim_dat <- bim_dat[variant_12 %in% sum_stats_dat$variant | variant_21 %in% sum_stats_dat$variant]

fwrite(bim_dat, file = "resources/ukbb_sum_stats/plink_ranges/chrX.txt", row.names = F, sep = ' ', col.names = F, quote = F)
