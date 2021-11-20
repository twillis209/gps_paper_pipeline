library(data.table)
library(argparse)

parser <- ArgumentParser(description = 'Add trait labels to pvalue file')
parser$add_argument('-p', '--pvalue_file', type = 'character', help = 'P-value file path')
parser$add_argument('-l', '--lookup_file', type = 'character', help = 'Path to file containing trait-label lookup')
parser$add_argument('-o', '--output_path', type = 'character', help = 'Path to output file', required = T)

args <- parser$parse_args()

pvalue_dat <- fread(args$pvalue_file, sep = '\t', header = T)

lookup_dat <- fread(args$lookup_file, sep = '\t', header = T)

pvalue_dat <- merge(pvalue_dat, lookup_dat[,.(code, abbrv)], all.x = T, by.x = 'trait_A', by.y = 'code')
setnames(pvalue_dat, "abbrv", "abbrv_A")
pvalue_dat <- merge(pvalue_dat, lookup_dat[,.(code, abbrv)], all.x = T, by.x = 'trait_B', by.y = 'code')
setnames(pvalue_dat, "abbrv", "abbrv_B")

pvalue_dat <- pvalue_dat[, .(trait_A, trait_B, abbrv_A, abbrv_B, gps, n, loc, loc.sd, scale, scale.sd, shape, shape.sd, pval)]

fwrite(pvalue_dat, file = args$output_path, sep = '\t', col.names = T, row.names = F)
