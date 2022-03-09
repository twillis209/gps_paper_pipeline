library(data.table)
library(argparse)

parser <- ArgumentParser(description = 'Generate summary statistics file containing only causal variants')
parser$add_argument('--sum_stats_file', type = 'character', help = 'Path to summary statistics file')
parser$add_argument('--block_file', type = 'character', help = 'Path to block file')
parser$add_argument('-o', '--output_file', type = 'character', help = 'Path to output file', required = T)
parser$add_argument('-nt', '--no_of_threads', type = 'integer', help = 'Number of threads to use', default = 1)

#test_args <- c("--sum_stats_file", "results/simgwas/simulated_sum_stats/whole_genome_sum_stats/10000_10000/null_sum_stats.tsv.gz", "--block_file", "resources/ldetect/blocks.txt", "-o", "cv_dat.tsv.gz", "-nt", 4)

args <- parser$parse_args()

setDTthreads(args$no_of_threads)

sum_stats_dat <- fread(args$sum_stats_file, sep = '\t', header = T, select = c('position', 'a0', 'a1', 'chr', 'rsID', 'EUR'))

blocks_dat <- fread(args$block_file, sep = " ", header = F, col.names = c("chr_block", "block_chr", "start", "stop"))

blocks_dat[, block := tstrsplit(chr_block, "_block")[[2]]]

for(i in 1:nrow(sum_stats_dat)) {
  block_no <- blocks_dat[block_chr == sum_stats_dat[i, chr] & sum_stats_dat[i, position] >= start & sum_stats_dat[i, position] <= stop, block]

  if(length(block_no) > 1) {
    print(i)
  }

  sum_stats_dat[i, block := block_no[1]]
}

fwrite(sum_stats_dat, file = args$output_file, sep = '\t')
