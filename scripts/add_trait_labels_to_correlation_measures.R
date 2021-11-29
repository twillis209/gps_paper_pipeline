library(data.table)
library(argparse)

parser <- ArgumentParser(description = 'Add trait labels to correlation measures')
parser$add_argument('-p', '--corr_file', type = 'character', help = 'Correlation measures file path')
parser$add_argument('-l', '--lookup_file', type = 'character', help = 'Path to file containing trait-label lookup')
parser$add_argument('-o', '--output_path', type = 'character', help = 'Path to output file', required = T)

args <- parser$parse_args()

corr_dat <- fread(args$corr_file, sep = '\t', header = T)

lookup_dat <- fread(args$lookup_file, sep = '\t', header = T)

corr_dat <- merge(corr_dat, lookup_dat[,.(code, abbrv)], all.x = T, by.x = 'trait_A', by.y = 'code')
setnames(corr_dat, "abbrv", "abbrv_A")
corr_dat <- merge(corr_dat, lookup_dat[,.(code, abbrv)], all.x = T, by.x = 'trait_B', by.y = 'code')
setnames(corr_dat, "abbrv", "abbrv_B")

fwrite(corr_dat, file = args$output_path, sep = '\t', col.names = T, row.names = F)
