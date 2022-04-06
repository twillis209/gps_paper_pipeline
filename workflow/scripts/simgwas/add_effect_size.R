library(data.table)
library(argparse)

parser <- ArgumentParser(description = 'Add effect size to block file')
parser$add_argument('--block_file', type = 'character', help = 'Path to block file')
parser$add_argument('--effect_size', type = 'character', help = 'Effect size')
parser$add_argument('-o', '--output_file', type = 'character', help = 'Path to output file', required = T)
parser$add_argument('-nt', '--no_of_threads', type = 'integer', help = 'Number of threads to use', default = 1)

args <- parser$parse_args()

setDTthreads(args$no_of_threads)

block_dat <- fread(file = args$block_file, sep = '\t', header = F)

block_dat[, block_effect_size := args$effect_size]

fwrite(block_dat, file = args$output_file, sep = '\t', col.names = F)
