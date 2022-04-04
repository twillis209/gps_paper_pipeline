library(data.table)
library(argparse)

parser <- ArgumentParser(description = 'Prune simulated summary statistics file')
parser$add_argument('--sum_stats_file_A', type = 'character', help = 'Path to first summary statistics file')
parser$add_argument('--sum_stats_file_B', type = 'character', help = 'Path to second summary statistics file')
parser$add_argument('--no_reps', type = 'integer', help = 'No. of replicates')
parser$add_argument('-o', '--output_path', type = 'character', help = 'Path to merged summary statistics file', required = T)
parser$add_argument('-nt', '--no_of_threads', type = 'integer', help = 'Number of threads to use', default = 1)

args <- parser$parse_args()

setDTthreads(args$no_of_threads)

col_names <- c("position", "block", "a0", "a1", "TYPE", "EUR", "zexp", paste0("zsim.", 1:args$no_reps), paste0("vbetasim.", 1:args$no_reps), paste0("betasim.", 1:args$no_reps), paste0("p.", 1:args$no_reps), "chosen_or", "ncases", "ncontrols", "rsID", "chr")

sum_stats_A_dat <- fread(args$sum_stats_file_A, sep = '\t', header = F)
sum_stats_B_dat <- fread(args$sum_stats_file_B, sep = '\t', header = F)

names(sum_stats_A_dat) <- col_names
names(sum_stats_B_dat) <- col_names

setnames(sum_stats_A_dat, paste0('zsim.', 1:args$no_reps), sprintf('zsim.%s.A', 1:args$no_reps))
setnames(sum_stats_B_dat, paste0('zsim.', 1:args$no_reps), sprintf('zsim.%s.B', 1:args$no_reps))
setnames(sum_stats_A_dat, paste0('vbetasim.', 1:args$no_reps), sprintf('vbetasim.%s.A', 1:args$no_reps))
setnames(sum_stats_B_dat, paste0('vbetasim.', 1:args$no_reps), sprintf('vbetasim.%s.B', 1:args$no_reps))
setnames(sum_stats_A_dat, paste0('betasim.', 1:args$no_reps), sprintf('betasim.%s.A', 1:args$no_reps))
setnames(sum_stats_B_dat, paste0('betasim.', 1:args$no_reps), sprintf('betasim.%s.B', 1:args$no_reps))
setnames(sum_stats_A_dat, paste0('p.', 1:args$no_reps), sprintf('p.%s.A', 1:args$no_reps))
setnames(sum_stats_B_dat, paste0('p.', 1:args$no_reps), sprintf('p.%s.B', 1:args$no_reps))

B_cols <- c('rsID', sprintf('zsim.%s.B', 1:args$no_reps), sprintf('vbetasim.%s.B', 1:args$no_reps), sprintf('betasim.%s.B', 1:args$no_reps), sprintf('p.%s.B', 1:args$no_reps), 'chr', 'position', 'a0', 'a1')

merged_dat <- merge(sum_stats_A_dat, sum_stats_B_dat[, ..B_cols], by = c('chr', 'position', 'a0', 'a1'))

fwrite(merged_dat, file = args$output_path, sep = '\t')
