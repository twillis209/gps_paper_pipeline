library(data.table)

setDTthreads(snakemake@threads)

sum_stats_dat <- fread(snakemake@input[['sum_stats_file']], sep = '\t', header = T)

if(snakemake@params[['sans_mhc']]) {
  sum_stats_dat[, c('chr', 'bp') := tstrsplit(variant, split = ':', keep = 1:2)]

  sum_stats_dat <- sum_stats_dat[!(chr == 6 & bp %between% c(24e6, 45e6))]

  sum_stats_dat[, c('chr', 'bp') := NULL]
}

pruned_rsid_dat <- fread(snakemake@input[['pruned_range_file']], sep = ' ', header = F)

names(pruned_rsid_dat) <- 'ID'

bim_dat <- fread(snakemake@input[['bim_file']], sep = '\t', header = F, col.names = c('chr', 'ID', 'Cm', 'bp', 'A1', 'A2'))

# Prune the rsIDs
bim_dat <- bim_dat[ID %in% pruned_rsid_dat$ID]

bim_dat[, variant_12 := paste(chr, bp, A1, A2, sep = ':')]
bim_dat[, variant_21 := paste(chr, bp, A2, A1, sep = ':')]

# Prune the summary statistics; there are no rsIDs in this file so we need to construct the IDs from coordinates and alleles contained in the concatenated bim file
sum_stats_dat <- sum_stats_dat[variant %in% bim_dat$variant_12 | variant %in% bim_dat$variant_21]

fwrite(sum_stats_dat, file = snakemake@output[[1]], sep = '\t')

system(sprintf("sed -i 's/pval\\.//g' %s", snakemake@output[[1]]))
