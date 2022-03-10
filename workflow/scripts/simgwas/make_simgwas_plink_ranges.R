library(data.table)
setDTthreads(8)
library(argparse)

parser <- ArgumentParser(description = 'Produce PLINK range files')
parser$add_argument('--sum_stats_file', type = 'character', help = 'Path to summary statistics file')
parser$add_argument('--input_bim_file', type = 'character', help = 'Path to combined bim file')
parser$add_argument('-o', '--output_range_files', type = 'character', nargs = '+', help = 'List of output file paths', required = T)
parser$add_argument('-nt', '--no_of_threads', type = 'integer', help = 'Number of threads to use', default = 1)

test_args <- c("--sum_stats_file", "results/simgwas/simulated_sum_stats/whole_genome_sum_stats/10000_10000/null_sum_stats.tsv.gz", "--input_bim_file", "resources/1000g/euro/qc/chr1-22_qc.bim", "-o", "test.bim", "-nt", 4)

args <- parser$parse_args()

setDTthreads(args$no_of_threads)

sum_stats_dat <- fread(args$sum_stats_file, sep = '\t', header = T, select = c('position', 'chr', 'a0', 'a1'))

bim_dat <- fread(args$input_bim_file, sep = '\t', header = F, col.names = c('chr', 'ID', 'Cm', 'bp', 'A1', 'A2'))

merged_dat <- merge(sum_stats_dat, bim_dat, by.x = c('chr', 'position'), by.y = c('chr', 'bp'))

merged_dat <- merged_dat[(a0 == A1 & a1 == A2) | (a0 == A2 & a1 == A1)]

merged_dat[, a0 := NULL]
merged_dat[, a1 := NULL]

merged_dat[, ID := paste(chr, position, A1, A2, sep = ':')]

for(i in 1:22) {
  fwrite(merged_dat[chr == i], file = args$output_range_files[i], row.names = F, sep = ' ', col.names = F, quote = F)
}
