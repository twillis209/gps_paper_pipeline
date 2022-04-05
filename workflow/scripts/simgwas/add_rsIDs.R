library(data.table)
library(argparse)
library(magrittr)

parser <- ArgumentParser(description = 'Add rsID column')
parser$add_argument('--leg_file', type = 'character', help = 'Path to legend file')
parser$add_argument('--block_file', type = 'character', help = 'Path to block file')
parser$add_argument('--no_reps', type = 'integer', help = 'No. of replicates', default = 1)
parser$add_argument('-o', '--output_file', type = 'character', help = 'Path to output file', required = T)
parser$add_argument('-nt', '--no_of_threads', type = 'integer', help = 'Number of threads to use', default = 1)

test_args <- c("--leg_file", "resources/simgwas/1000g/1000GP_Phase3_chr1_eur_common_maf.legend.gz", "--block_file", "results/simgwas/simulated_sum_stats/chr1/block_sum_stats/null/10000_10000/block_7_sum_stats.tsv.gz", "--no_reps", 20, "-o", "test.tsv.gz", "-nt", 4)

args <- parser$parse_args()

setDTthreads(args$no_of_threads)

block_dat <- fread(file = args$block_file, sep = '\t', header = F)

# Old, incorrect columns
#names(block_dat) <- c("position", "block", "a0", "a1", "TYPE", "EUR", "zexp", paste0("zsim.", 1:args$no_reps), paste0("vbetasim.", 1:args$no_reps), paste0("betasim.", 1:args$no_reps), paste0("p.", 1:args$no_reps), "chosen_or", "ncases", "ncontrols", "id", "chr")

names(block_dat) <- c("position", "a0", "a1", "id", "block", "TYPE", "EUR", "zexp", paste0("zsim.", 1:20), paste0("vbetasim.", 1:20), paste0("betasim.", 1:20), paste0("p.", 1:20), "chosen_or", "ncases", "ncontrols", "chr")

block_dat[, id := NULL]

leg_dat <- fread(file = args$leg_file, sep = ' ', header = T, select = c('id', 'position', 'a0', 'a1'))

merged_dat <- merge(leg_dat, block_dat, by = c('position', 'a0', 'a1'))

merged_dat[, id := tstrsplit(id, split = ':', keep = 1)]

fwrite(merged_dat, file = args$output_file, sep = '\t', col.names = F)
