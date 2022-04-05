library(data.table)
setDTthreads(8)
library(argparse)

parser <- ArgumentParser(description = 'Produce PLINK range files')
parser$add_argument('--sum_stats_file', type = 'character', help = 'Path to summary statistics file')
parser$add_argument('--input_bim_file', type = 'character', help = 'Path to combined bim file')
parser$add_argument('--bp_pos', type = 'integer', help = 'Column no. of bp/position column')
parser$add_argument('--chr_pos', type = 'integer', help = 'Column no. of chr column')
parser$add_argument('--a0_pos', type = 'integer', help = 'Column no. of a0 column')
parser$add_argument('--a1_pos', type = 'integer', help = 'Column no. of a1 column')
parser$add_argument('-o', '--output_range_files', type = 'character', nargs = '+', help = 'List of output file paths', required = T)
parser$add_argument('-nt', '--no_of_threads', type = 'integer', help = 'Number of threads to use', default = 1)

args <- parser$parse_args()

setDTthreads(args$no_of_threads)

sum_stats_dat <- fread(args$sum_stats_file, sep = '\t', header = F, select = c(args$bp_pos, args$chr_pos, args$a0_pos, args$a1_pos))

names(sum_stats_dat) <- c('position', 'chr', 'a0', 'a1')

bim_dat <- fread(args$input_bim_file, sep = '\t', header = F, col.names = c('chr', 'id', 'Cm', 'bp', 'A1', 'A2'))

merged_dat <- merge(sum_stats_dat, bim_dat, by.x = c('chr', 'position'), by.y = c('chr', 'bp'))

merged_dat <- merged_dat[(a0 == A1 & a1 == A2) | (a0 == A2 & a1 == A1)]

merged_dat[, a0 := NULL]
merged_dat[, a1 := NULL]

for(i in 1:22) {
  fwrite(merged_dat[chr == i], file = args$output_range_files[i], row.names = F, sep = ' ', col.names = F, quote = F)
}
