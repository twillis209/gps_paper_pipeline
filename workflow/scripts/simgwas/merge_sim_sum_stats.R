library(data.table)
library(argparse)

parser <- ArgumentParser(description = 'Prune simulated summary statistics file')
parser$add_argument('--sum_stats_file_A', type = 'character', help = 'Path to first summary statistics file')
parser$add_argument('--sum_stats_file_B', type = 'character', help = 'Path to second summary statistics file')
parser$add_argument('--no_reps', type = 'integer', help = 'No. of replicates')
parser$add_argument('-o', '--output_path', type = 'character', help = 'Path to merged summary statistics file', required = T)
parser$add_argument('-nt', '--no_of_threads', type = 'integer', help = 'Number of threads to use', default = 1)

test_args <- c('--sum_stats_file_A', 'results/simgwas/simulated_sum_stats/chr1/whole_chr_sum_stats/10000_10000/blocks_0_121_v2_sum_stats.tsv.gz',
               '--sum_stats_file_B', 'results/simgwas/simulated_sum_stats/chr1/whole_chr_sum_stats/10000_10000/blocks_0_121_s2_sum_stats.tsv.gz',
               '-o', 'results/simgwas/simulated_sum_stats/merged/chr1_1/10000_10000_10000_10000/blocks_0_121_v2_blocks_0_121_s2_sum_stats.tsv.gz',
               '-nt', 4)

args <- parser$parse_args()

setDTthreads(args$no_of_threads)

sum_stats_A_dat <- fread(args$sum_stats_file_A, sep = '\t', header = T)
sum_stats_B_dat <- fread(args$sum_stats_file_B, sep = '\t', header = T)

setnames(sum_stats_A_dat, paste0('zsim.', 1:args$no_reps), sprintf('zsim.%s.A', 1:args$no_reps))
setnames(sum_stats_B_dat, paste0('zsim.', 1:args$no_reps), sprintf('zsim.%s.B', 1:args$no_reps))
setnames(sum_stats_A_dat, paste0('vbetasim.', 1:args$no_reps), sprintf('vbetasim.%s.A', 1:args$no_reps))
setnames(sum_stats_B_dat, paste0('vbetasim.', 1:args$no_reps), sprintf('vbetasim.%s.B', 1:args$no_reps))
setnames(sum_stats_A_dat, paste0('betasim.', 1:args$no_reps), sprintf('betasim.%s.A', 1:args$no_reps))
setnames(sum_stats_B_dat, paste0('betasim.', 1:args$no_reps), sprintf('betasim.%s.B', 1:args$no_reps))
setnames(sum_stats_A_dat, paste0('p.', 1:args$no_reps), sprintf('p.%s.A', 1:args$no_reps))
setnames(sum_stats_B_dat, paste0('p.', 1:args$no_reps), sprintf('p.%s.B', 1:args$no_reps))

B_cols <- c('rsID', sprintf('zsim.%s.B', 1:args$no_reps), sprintf('vbetasim.%s.B', 1:args$no_reps), sprintf('betasim.%s.B', 1:args$no_reps), sprintf('p.%s.B', 1:args$no_reps))

merged_dat <- merge(sum_stats_A_dat, sum_stats_B_dat[, ..B_cols], by = 'rsID')

fwrite(merged_dat, file = args$output_path, sep = '\t')