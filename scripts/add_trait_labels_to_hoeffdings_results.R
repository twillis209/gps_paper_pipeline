library(data.table)
library(argparse)

parser <- ArgumentParser(description = 'Add trait labels to Hoeffding\'s test results')
parser$add_argument('-r', '--results_file', type = 'character', help = 'Results file path')
parser$add_argument('-l', '--lookup_file', type = 'character', help = 'Path to file containing trait-label lookup')
parser$add_argument('-o', '--output_path', type = 'character', help = 'Path to output file', required = T)

args <- parser$parse_args()

results_dat <- fread(args$results_file, sep = '\t', header = T)

lookup_dat <- fread(args$lookup_file, sep = '\t', header = T)

results_dat <- merge(results_dat, lookup_dat[,.(code, abbrv)], all.x = T, by.x = 'trait_A', by.y = 'code')
setnames(results_dat, "abbrv", "abbrv_A")
results_dat <- merge(results_dat, lookup_dat[,.(code, abbrv)], all.x = T, by.x = 'trait_B', by.y = 'code')
setnames(results_dat, "abbrv", "abbrv_B")

results_dat <- results_dat[, c("trait_A", "trait_B", "abbrv_A", "abbrv_B", "n", "Dn", "scaled", "p.value")]

fwrite(results_dat, file = args$output_path, sep = '\t', col.names = T, row.names = F)
