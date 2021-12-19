library(data.table)
library(argparse)

parser <- ArgumentParser(description = 'Excise MHC from pruned, merged GWAS summary statistics file')
parser$add_argument('-i', '--sum_stats_file', type = 'character', help = 'Path to summary statistics file')
parser$add_argument('-n', '--no_snps', type = 'integer', help = 'Number of SNPs to retain in output')
parser$add_argument('-o', '--output_path', type = 'character', help = 'Path to output file', required = T)
parser$add_argument('-nt', '--no_of_threads', type = 'integer', help = 'Number of threads to use', default = 1)

args <- parser$parse_args()

setDTthreads(args$no_of_threads)

sum_stats_dat <- fread(args$sum_stats_file, sep = '\t', header = T)

sum_stats_dat <- sum_stats_dat[!(CHR19 == 6 & BP19 %between% c(24e6, 45e6))]

fwrite(sum_stats_dat, file = args$output_path, sep = '\t')
