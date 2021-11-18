library(data.table)
library(argparse)

parser <- ArgumentParser(description = 'Prune merged GWAS summary statistics file')
parser$add_argument('-is', '--sum_stats_file', type = 'character', help = 'Path to summary statistics file')
parser$add_argument('-ib', '--bim_file', type = 'character', help = 'Path to bim file')
parser$add_argument('-ir', '--range_file', type = 'character', help = 'Path to range file')
parser$add_argument('-o', '--output_path', type = 'character', help = 'Path to pruned summary statistics file', required = T)
parser$add_argument('-nt', '--no_of_threads', type = 'integer', help = 'Number of threads to use', default = 1)

args <- parser$parse_args()

setDTthreads(args$no_of_threads)

sum_stats_dat <- fread(args$sum_stats_file, sep = '\t', header = T)

pruned_rsid_dat <- fread(args$range_file, sep = ' ', header = F)

names(pruned_rsid_dat) <- 'ID'

bim_dat <- fread(args$bim_file, sep = '\t', header = F, col.names = c('chr', 'ID', 'Cm', 'bp', 'A1', 'A2'))

# Prune the rsIDs
bim_dat <- bim_dat[ID %in% pruned_rsid_dat$ID]

bim_dat[, variant_12 := paste(chr, bp, A1, A2, sep = ':')]
bim_dat[, variant_21 := paste(chr, bp, A2, A1, sep = ':')]

# Prune the summary statistics; there are no rsIDs in this file so we need to construct the IDs from coordinates and alleles contained in the concatenated bim file
sum_stats_dat <- sum_stats_dat[variant %in% bim_dat$variant_12 | variant %in% bim_dat$variant_21]

fwrite(sum_stats_dat, file = args$output_path, sep = '\t')
