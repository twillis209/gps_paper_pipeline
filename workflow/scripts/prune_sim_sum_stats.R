library(data.table)
library(argparse)

parser <- ArgumentParser(description = 'Prune simulated summary statistics file')
parser$add_argument('--sum_stats_file', type = 'character', help = 'Path to summary statistics file')
parser$add_argument('--bim_file', type = 'character', help = 'Path to bim file')
parser$add_argument('--prune_file', type = 'character', help = 'Path to file containing pruned IDs')
parser$add_argument('-o', '--output_path', type = 'character', help = 'Path to pruned summary statistics file', required = T)
parser$add_argument('-nt', '--no_of_threads', type = 'integer', help = 'Number of threads to use', default = 1)

test_args <- c('--sum_stats_file', #'results/simgwas/simulated_sum_stats/20000_snps_1_cv_sum_stats.tsv.gz',
               'test.tsv', '--bim_file', 'resources/1000g/euro/qc/chr21_qc.bim', '--prune_file', 'resources/1000g/euro/qc/pruned_ranges/window_1000kb_step_50/chr21.prune.in', '-o', 'test.tsv.gz', '-nt', 8)

args <- parser$parse_args()

setDTthreads(args$no_of_threads)

sum_stats_dat <- fread(args$sum_stats_file, sep = '\t', header = T)

pruned_rsid_dat <- fread(args$prune_file, sep = ' ', header = F)

names(pruned_rsid_dat) <- 'ID'

bim_dat <- fread(args$bim_file, sep = '\t', header = F, col.names = c('chr', 'ID', 'Cm', 'bp', 'A1', 'A2'))

pruned_rsid_dat <- merge(bim_dat, pruned_rsid_dat, by = 'ID')

merged_dat <- merge(sum_stats_dat, pruned_rsid_dat[ , .(bp, A1, A2)], by.x = c('position', 'a0', 'a1'), by.y = c('bp', 'A1', 'A2'))

fwrite(merged_dat, file = args$output_path, sep = '\t')
