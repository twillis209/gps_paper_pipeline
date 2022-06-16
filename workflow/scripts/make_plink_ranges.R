library(data.table)
setDTthreads(8)
library(argparse)

parser <- ArgumentParser(description = 'Produce PLINK range files')
parser$add_argument('-is', '--sum_stats_file', type = 'character', help = 'Path to summary statistics file')
parser$add_argument('-ib', '--input_bim_files', type = 'character', nargs = '+', help = 'List of bim file paths')
parser$add_argument('-nc', '--non_na_cols', type = 'character', nargs = '+', help = 'Comma-delimited list of columns to assert as non-NA')
parser$add_argument('-o', '--output_bim_files', type = 'character', nargs = '+', help = 'List of output file paths', required = T)
parser$add_argument('-nt', '--no_of_threads', type = 'integer', help = 'Number of threads to use', default = 1)

args <- parser$parse_args()

setDTthreads(args$no_of_threads)

if(length(args$output_bim_files) != length(args$input_bim_files)) {
  stop("Length of input and output bim file vectors must match")
}

sum_stats_dat <- fread(args$sum_stats_file, sep = '\t', header = T)

# Retain only those SNPs for which we have p-values in all relevant traits
for(x in args$non_na_cols) {
  sum_stats_dat <- sum_stats_dat[!is.na(get(x))]
}

for(i in seq_along(args$input_bim_files)) {
  bim_dat <- fread(args$input_bim_files[i], sep = '\t', header = F, col.names = c('chr', 'ID', 'Cm', 'bp', 'A1', 'A2'))

  bim_dat[, variant_12 := paste(chr, bp, A1, A2, sep = ':')]
  bim_dat[, variant_21 := paste(chr, bp, A2, A1, sep = ':')]

  bim_dat <- bim_dat[variant_12 %in% sum_stats_dat$variant | variant_21 %in% sum_stats_dat$variant]

  fwrite(bim_dat, file = args$output_bim_files[i], row.names = F, sep = ' ', col.names = F, quote = F)
}
