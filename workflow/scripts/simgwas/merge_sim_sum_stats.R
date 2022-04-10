library(data.table)
library(argparse)

parser <- ArgumentParser(description = 'Prune simulated summary statistics file')
parser$add_argument('--sum_stats_file_A', type = 'character', help = 'Path to first summary statistics file')
parser$add_argument('--sum_stats_file_B', type = 'character', help = 'Path to second summary statistics file')
parser$add_argument('--file_A_stat_cols', type = 'character', nargs = '+', help = 'Columns to retain from file A')
parser$add_argument('--file_B_stat_cols', type = 'character', nargs = '+', help = 'Columns to retain from file B')
parser$add_argument('-o', '--output_path', type = 'character', help = 'Path to merged summary statistics file', required = T)
parser$add_argument('-nt', '--no_of_threads', type = 'integer', help = 'Number of threads to use', default = 1)

test_args <- c("--sum_stats_file_A", "results/simgwas/simulated_sum_stats/whole_genome_sum_stats/randomised/10000_10000_10000_10000/m50_m50_m10/m50_seed_353_sum_stats_A_tag_c_of_cd.tsv.gz", "--sum_stats_file_B", "results/simgwas/simulated_sum_stats/whole_genome_sum_stats/randomised/10000_10000_10000_10000/m50_m50_m10/m50_seed_353_sum_stats_B_tag_d_of_cd.tsv.gz", "--file_A_stat_cols", "p.3", "--file_B_stat_cols", "p.4", "-o", "results/simgwas/simulated_sum_stats/whole_genome_sum_stats/randomised/10000_10000_10000_10000/m50_m50_m10/seed_353_merged_sum_stats_tags_cd.tsv.gz", "-nt", 4)

test_args <- c("--sum_stats_file_A", "results/simgwas/simulated_sum_stats/whole_genome_sum_stats/randomised/10000_10000_10000_10000/m50_m50_m0/m50_seed_400_sum_stats_A_tag_k_of_kl.tsv.gz", "--sum_stats_file_B", "results/simgwas/simulated_sum_stats/whole_genome_sum_stats/randomised/10000_10000_10000_10000/m50_m50_m0/m50_seed_400_sum_stats_B_tag_l_of_kl.tsv.gz", "--file_A_stat_cols", "p.11", "--file_B_stat_cols", "p.12", "-o", "results/simgwas/simulated_sum_stats/whole_genome_sum_stats/randomised/10000_10000_10000_10000/m50_m50_m0/seed_400_merged_sum_stats_tags_kl.tsv.gz", "-nt", 4)

args <- parser$parse_args(test_args)

args <- parser$parse_args()

setDTthreads(args$no_of_threads)

sum_stats_A_dat <- fread(args$sum_stats_file_A, sep = '\t', header = T, select = c('id', 'chr', 'position', 'a0', 'a1', 'block', args$file_A_stat_cols, 'block_effect_size'))
sum_stats_B_dat <- fread(args$sum_stats_file_B, sep = '\t', header = T, select = c('id', 'chr', 'position', 'a0', 'a1', 'block', args$file_B_stat_cols, 'block_effect_size'))

merged_dat <- merge(sum_stats_A_dat, sum_stats_B_dat, by = c('id', 'chr', 'position', 'a0', 'a1', 'block'), suffixes = c(".A", ".B"))

if(length(unique(merged_dat$chr)) != 22) stop("Incorrect number of chromosomes present in merged summary statistics")

fwrite(merged_dat, file = args$output_path, sep = '\t')
