library(data.table)
library(argparse)

parser <- ArgumentParser(description = 'Replace replicates in one file with those from another')
parser$add_argument('--sum_stats_file_A', type = 'character', help = 'Path to first summary statistics file')
parser$add_argument('--sum_stats_file_B', type = 'character', help = 'Path to second summary statistics file')
parser$add_argument('--no_reps_A', type = 'integer', help = 'Number of replicates in file A')
parser$add_argument('--no_reps_B', type = 'integer', help = 'Number of replicates in file B')
parser$add_argument('--file_A_stat_cols', type = 'character', nargs = '+', help = 'Columns to replace in file A')
parser$add_argument('--file_B_stat_cols', type = 'character', nargs = '+', help = 'Columns to copy from file B to file A')
parser$add_argument('-o', '--output_path', type = 'character', help = 'Path to merged summary statistics file', required = T)
parser$add_argument('-nt', '--no_of_threads', type = 'integer', help = 'Number of threads to use', default = 1)

test_args <- c('--sum_stats_file_A', 'results/simgwas/simulated_sum_stats/chr12/block_sum_stats/null/500_10000/block_54_sum_stats.tsv.gz',
               '--sum_stats_file_B', 'results/simgwas/simulated_sum_stats/rep_blocks/chr12/block_sum_stats/null/500_10000/block_54_sum_stats.tsv.gz',
               '--no_reps_A', 20,
               '--no_reps_B', 2,
               '--file_A_stat_cols', 'p.13', 'p.14',
               '--file_B_stat_cols', 'p.1', 'p.2',
               '-o', 'test.tsv.gz',
               '-nt', 1)

args <- parser$parse_args(test_args)

args <- parser$parse_args()

setDTthreads(args$no_of_threads)

sum_stats_A_dat <- fread(args$sum_stats_file_A, sep = '\t', header = F)
sum_stats_B_dat <- fread(args$sum_stats_file_B, sep = '\t', header = F)

col_names_A <- c("position", "a0", "a1", "id", "block", "TYPE", "EUR", "zexp", paste0("zsim.", 1:args$no_reps_A), paste0("vbetasim.", 1:args$no_reps_A), paste0("betasim.", 1:args$no_reps_A), paste0("p.", 1:args$no_reps_A), "chosen_or", "ncases", "ncontrols", "chr", "block_effect_size")

col_names_B <- c("position", "a0", "a1", "id", "block", "TYPE", "EUR", "zexp", paste0("zsim.", 1:args$no_reps_B), paste0("vbetasim.", 1:args$no_reps_B), paste0("betasim.", 1:args$no_reps_B), paste0("p.", 1:args$no_reps_B), "chosen_or", "ncases", "ncontrols", "chr", "block_effect_size")

names(sum_stats_A_dat) <- col_names_A
names(sum_stats_B_dat) <- col_names_B

if(any(sum_stats_A_dat$id != sum_stats_B_dat$id)) {
  stop("sum_stats_A_dat$id != sum_stats_B_dat$id")
}

sum_stats_A_dat[, args$file_A_stat_cols := NULL]

# Need to rename these columns to the names we wish them to have in the A data.table
setnames(sum_stats_B_dat, args$file_B_stat_cols, args$file_A_stat_cols)

cols_to_replace <- args$file_A_stat_cols

sum_stats_A_dat <- cbind(sum_stats_A_dat, sum_stats_B_dat[ , ..cols_to_replace])

# Set correct column order
sum_stats_A_dat <- sum_stats_A_dat[, ..col_names_A]

fwrite(sum_stats_A_dat, file = args$output_path, sep = '\t', col.names = F)
