library(data.table)
library(argparse)

parser <- ArgumentParser(description = 'Prune simulated summary statistics file')
parser$add_argument('--sum_stats_file_A', type = 'character', help = 'Path to first summary statistics file')
parser$add_argument('--sum_stats_file_B', type = 'character', help = 'Path to second summary statistics file')
parser$add_argument('--file_A_stat_cols', type = 'character', nargs = '+', help = 'Columns to retain from file A')
parser$add_argument('--file_B_stat_cols', type = 'character', nargs = '+', help = 'Columns to retain from file B')
parser$add_argument('--sum_stats_file_B', type = 'character', help = 'Path to second summary statistics file')
parser$add_argument('-o', '--output_path', type = 'character', help = 'Path to merged summary statistics file', required = T)
parser$add_argument('-nt', '--no_of_threads', type = 'integer', help = 'Number of threads to use', default = 1)

args <- parser$parse_args()

setDTthreads(args$no_of_threads)

sum_stats_A_dat <- fread(args$sum_stats_file_A, sep = '\t', header = T, select = c('id', 'chr', 'position', 'a0', 'a1', 'block', 'block_effect_size', args$file_A_stat_cols))
sum_stats_B_dat <- fread(args$sum_stats_file_B, sep = '\t', header = T, select = c('id', 'chr', 'position', 'a0', 'a1', 'block', 'block_effect_size', args$file_B_stat_cols))

merged_dat <- merge(sum_stats_A_dat, sum_stats_B_dat, by = c('id', 'chr', 'position', 'a0', 'a1', 'block'), suffixes = c(".A", ".B"))

if(length(unique(merged_dat$chr)) != 22) stop("Incorrect number of chromosomes present in merged summary statistics")

fwrite(merged_dat, file = args$output_path, sep = '\t')
