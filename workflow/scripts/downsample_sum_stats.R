library(data.table)
library(argparse)

parser <- ArgumentParser(description = 'Prune merged GWAS summary statistics file')
parser$add_argument('-i', '--sum_stats_file', type = 'character', help = 'Path to summary statistics file')
parser$add_argument('-n', '--no_snps', type = 'integer', help = 'Number of SNPs to retain in output')
parser$add_argument('-o', '--output_path', type = 'character', help = 'Path to output file', required = T)
parser$add_argument('-nt', '--no_of_threads', type = 'integer', help = 'Number of threads to use', default = 1)

args <- parser$parse_args()

setDTthreads(args$no_of_threads)

sum_stats_dat <- fread(args$sum_stats_file, sep = '\t', header = T)

args$no_snps <- min(nrow(sum_stats_dat), args$no_snps)

sum_stats_dat <- sum_stats_dat[sample(1:args$no_snps)]

fwrite(sum_stats_dat, file = args$output_path, sep = '\t')
