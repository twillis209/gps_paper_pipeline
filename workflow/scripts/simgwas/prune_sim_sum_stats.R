library(data.table)
library(argparse)

parser <- ArgumentParser(description = 'Prune simulated summary statistics file')
parser$add_argument('--sum_stats_file', type = 'character', help = 'Path to summary statistics file')
parser$add_argument('--bim_file', type = 'character', help = 'Path to bim file')
parser$add_argument('--prune_file', type = 'character', help = 'Path to file containing pruned IDs')
parser$add_argument('-o', '--output_path', type = 'character', help = 'Path to pruned summary statistics file', required = T)
parser$add_argument('-nt', '--no_of_threads', type = 'integer', help = 'Number of threads to use', default = 1)

args <- c("--sum_stats_file", "results/simgwas/simulated_sum_stats/merged/10000_10000_10000_10000/1-m0:24_1-m10:34_sum_stats.tsv.gz", "--bim_file", "resources/1000g/euro/qc/chr1-22_qc.bim", "--prune_file", "resources/plink_ranges/simgwas/pruned_ranges/window_1000kb_step_50/all.prune.in", "-o",  "results/simgwas/simulated_sum_stats/pruned/window_1000kb_step_50/10000_10000_10000_10000/1-m0:24_1-m10:34_sum_stats.tsv", "-nt", 4)

args <- parser$parse_args()

setDTthreads(args$no_of_threads)

sum_stats_dat <- fread(args$sum_stats_file, sep = '\t', header = T)

pruned_rsid_dat <- fread(args$prune_file, sep = ' ', header = F, col.names = 'ID')

bim_dat <- fread(args$bim_file, sep = '\t', header = F, col.names = c('chr', 'ID', 'Cm', 'bp', 'A1', 'A2'))

pruned_rsid_dat <- merge(bim_dat, pruned_rsid_dat, by = 'ID')

rm(bim_dat)

merged_dat <- merge(sum_stats_dat, pruned_rsid_dat[ , .(chr, bp, A1, A2)], by.x = c('chr', 'position'), by.y = c('chr', 'bp'))

merged_dat <- merged_dat[(a0 == A2 & a1 == A1) | (a0 == A1 & a1 == A2)]

# Some indels are present as duplicates with flipped A1/A2 entries
merged_dat <- merged_dat[!duplicated(merged_dat, by = c('chr', 'position'))]

merged_dat[, c("A1", "A2") := NULL]

fwrite(merged_dat, file = args$output_path, sep = '\t')
