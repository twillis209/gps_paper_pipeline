library(data.table)
setDTthreads(8)

sum_stats_dat <- fread('resources/ukbb_sum_stats/merged_sum_stats.tsv', sep = '\t', header = T)

pruned_rsid_dat <- fread('resources/ukbb_sum_stats/plink_subsets/pruned_ranges/all.prune.in', sep = ' ', header = F)

names(pruned_rsid_dat) <- 'ID'

bim_dat <- fread('resources/ukbb_sum_stats/plink_subsets/all.bim', sep = '\t', header = F, col.names = c('chr', 'ID', 'Cm', 'bp', 'A1', 'A2'))

# Prune the rsIDs
bim_dat <- bim_dat[ID %in% pruned_rsid_dat$ID]

bim_dat[, variant_12 := paste(chr, bp, A1, A2, sep = ':')]
bim_dat[, variant_21 := paste(chr, bp, A2, A1, sep = ':')]

# Prune the summary statistics; there are no rsIDs in this file so we need to construct the IDs from coordinates and alleles contained in the concatenated bim file
sum_stats_dat <- sum_stats_dat[variant %in% bim_dat$variant_12 | variant %in% bim_dat$variant_21]

fwrite(sum_stats_dat, file = 'resources/ukbb_sum_stats/pruned_merged_sum_stats.tsv', sep = '\t')
