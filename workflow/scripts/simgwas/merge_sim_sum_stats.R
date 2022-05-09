library(data.table)
library(argparse)

parser <- ArgumentParser(description = 'Merge simulated summary statistics files')
parser$add_argument('--sum_stats_file_A', type = 'character', help = 'Path to first summary statistics file')
parser$add_argument('--sum_stats_file_B', type = 'character', help = 'Path to second summary statistics file')
parser$add_argument('--file_A_stat_cols', type = 'character', nargs = '+', help = 'Columns to retain from file A')
parser$add_argument('--file_B_stat_cols', type = 'character', nargs = '+', help = 'Columns to retain from file B')
parser$add_argument('-o', '--output_path', type = 'character', help = 'Path to merged summary statistics file', required = T)
parser$add_argument('-nt', '--no_of_threads', type = 'integer', help = 'Number of threads to use', default = 1)

#test_args <- c("--sum_stats_file_A", "results/simgwas/simulated_sum_stats/whole_genome_sum_stats/randomised/5000_10000_5000_10000/m50_m50_m0/m50_seed_2124_sum_stats_A_tag_e_of_ef.tsv.gz", "--sum_stats_file_B", "results/simgwas/simulated_sum_stats/whole_genome_sum_stats/randomised/5000_10000_5000_10000/m50_m50_m0/m50_seed_2124_sum_stats_B_tag_f_of_ef.tsv.gz", "--file_A_stat_cols", "p.5", "--file_B_stat_cols", "p.6", "-o", "results/simgwas/simulated_sum_stats/whole_genome_sum_stats/randomised/5000_10000_5000_10000/m50_m50_m0/seed_2124_merged_sum_stats_tags_ef.tsv.gz", "-nt", 1)

#args <- parser$parse_args(test_args)

args <- parser$parse_args()

setDTthreads(args$no_of_threads)

sum_stats_A_dat <- fread(args$sum_stats_file_A, sep = '\t', header = T, select = c('id', 'chr', 'position', 'a0', 'a1', 'block', args$file_A_stat_cols, 'block_effect_size'))
sum_stats_B_dat <- fread(args$sum_stats_file_B, sep = '\t', header = T, select = c('id', 'chr', 'position', 'a0', 'a1', 'block', args$file_B_stat_cols, 'block_effect_size'))

merged_dat <- merge(sum_stats_A_dat, sum_stats_B_dat, by = c('chr', 'position', 'a0', 'a1', 'block'), suffixes = c(".A", ".B"))

merged_dat[, id.B := NULL]

setnames(merged_dat, 'id.A', 'id')

if(length(unique(merged_dat$chr)) != 22) {
  fwrite(merged_dat, file = args$output_path, sep = '\t')
  stop(sprintf("Incorrect number of chromosomes present in merged summary statistics: %s", length(unique(merged_dat$chr))))
} else {
  fwrite(merged_dat, file = args$output_path, sep = '\t')
}

