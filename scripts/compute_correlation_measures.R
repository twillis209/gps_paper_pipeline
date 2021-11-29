library(data.table)
library(argparse)

parser <- ArgumentParser(description = 'Calculate correlation measures for a pair of traits')
parser$add_argument('-i', '--sum_stats_file', type = 'character', help = 'Path to summary statistics file')
parser$add_argument('-a', '--trait_A', type = 'character', help = 'Trait A code')
parser$add_argument('-b', '--trait_B', type = 'character', help = 'Trait B code')
parser$add_argument('-m', '--methods', type = 'character', help = 'Comma-delimited list of correlation measures', default = 'pearson,spearman')
parser$add_argument('-o', '--output_path', type = 'character', help = 'Path to output file', required = T)
parser$add_argument('-nt', '--no_of_threads', type = 'integer', help = 'Number of threads to use', default = 1)

args <- parser$parse_args()

setDTthreads(args$no_of_threads)

args$methods <- unlist(strsplit(args$methods, ','))

sum_stats_dat <- fread(args$sum_stats_file, sep = '\t', header = T, select = c(args$trait_A, args$trait_B))

res_dat <- data.table(trait_A = args$trait_A, trait_B = args$trait_B, n = nrow(sum_stats_dat))

for(x in args$methods) {
  corr <- cor(sum_stats_dat[[args$trait_A]], sum_stats_dat[[args$trait_B]], method = x)
  res_dat[, (x) := corr]
}

fwrite(res_dat, file = args$output_path, sep = '\t', col.names = T, row.names = F)
