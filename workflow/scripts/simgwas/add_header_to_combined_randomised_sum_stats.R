library(data.table)
library(argparse)

parser <- ArgumentParser(description = 'Add header to combined randomised summary statistics')
parser$add_argument('--sum_stats_file', type = 'character', help = 'Path to first summary statistics file')
parser$add_argument('--no_reps', type = 'integer', help = 'No. of replicates')
parser$add_argument('-o', '--output_path', type = 'character', help = 'Path to merged summary statistics file', required = T)
parser$add_argument('-nt', '--no_of_threads', type = 'integer', help = 'Number of threads to use', default = 1)

args <- parser$parse_args()

setDTthreads(args$no_of_threads)

sum_stats_dat <- fread(args$sum_stats_file, sep = '\t', header = F)

col_names <- c("position", "a0", "a1", "id", "block", "TYPE", "EUR", "zexp", paste0("zsim.", 1:args$no_reps), paste0("vbetasim.", 1:args$no_reps), paste0("betasim.", 1:args$no_reps), paste0("p.", 1:args$no_reps), "chosen_or", "ncases", "ncontrols", "chr")

names(sum_stats_dat) <- col_names

fwrite(sum_stats_dat, file = args$output_path, sep = '\t')
