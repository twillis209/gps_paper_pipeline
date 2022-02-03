library(data.table)
library(simGWAS)
library(argparse)
library(parallel)

parser <- ArgumentParser(description = 'Computes blockwise LD matrices for chromosome 21')
parser$add_argument('--hap_file', type = 'character', help = 'Path to haplotype file')
parser$add_argument('--leg_file', type = 'character', help = 'Path to legend file')
parser$add_argument('--bim_file', type = 'character', help = 'Path to bim file')
parser$add_argument('--block_file', type = 'character', help = 'Path to block file')
parser$add_argument('-o', '--output_file', type = 'character', help = 'Path to output file', required = T)
parser$add_argument('-nt', '--no_of_threads', type = 'integer', help = 'Number of threads to use', default = 1)

test_args <- c("--hap_file", "resources/simgwas/1000g/1000GP_Phase3_chr21_with_meta_eur.hap.gz", "--leg_file", "resources/simgwas/1000g/1000GP_Phase3_chr21.legend.gz", "--bim_file", "resources/1000g/chr21.bim", "-o", "results/simgwas/chr21_block_ld_matrices_alt.RData", "-nt", 8, "--block_file", "resources/ldetect/blocks.txt")

args <- parser$parse_args(test_args)

args <- parser$parse_args()

setDTthreads(args$no_of_threads)

leg_dat <- fread(file = args$leg_file, sep = ' ', header = T)
hap_dat <- fread(file = args$hap_file, sep = ' ', header = F)
bim_dat <- fread(file = args$bim_file, sep = '\t', header = F, col.names = c('chr', 'rsID', 'Cm', 'bp', 'A1', 'A2'))
block_dat <- fread(file = args$block_file, sep = ' ', header = F, col.names = c('block', 'chr', 'start', 'stop'))

# Drop metadata rows
hap_dat <- hap_dat[5:nrow(hap_dat)]

# Only using data from chr21
block_dat <- block_dat[chr == 21]

leg_dat[, rs := make.names(leg_dat$id)]

hap_dat[, rs := leg_dat$rs]

leg_dat <- leg_dat[EUR %between% c(0.01, 0.99)]

for(i in 1:nrow(block_dat)) {
   leg_dat[position %between% c(block_dat[i, start], block_dat[i, stop]), block := (i-1)]
}

# Drop SNPs we can't assign to a block
leg_dat <- leg_dat[!is.na(block)]

hap_dat <- hap_dat[rs %in% leg_dat$rs]

hap_dat[, block := leg_dat$block]

hap_dat <- hap_dat[, c(seq(1, ncol(hap_dat)-3, by = 2), seq(2, ncol(hap_dat)-2, by = 2), ncol(hap_dat)-1, ncol(hap_dat)), with = F]

for(j in 1:(ncol(hap_dat)-2)) {
  set(hap_dat, j = j, value = as.numeric(hap_dat[[j]]))
}

hap_mat <- as.matrix(hap_dat[,1:(ncol(hap_dat)-2)])

freq_dat <- data.table(t(hap_mat)+1)

colnames(freq_dat) <- hap_dat$rs

freq_dat[, Probability := 1/.N]

rm(hap_dat, hap_mat)

block_snps <- table(leg_dat$block)

ld_mat_rsIDs <- lapply(sort(unique(leg_dat$block)), function(i) leg_dat[block == i, rs])

compute_ld_matrix <- function(freq, rsIDs) {
  corpcor::make.positive.definite(simGWAS:::wcor2(as.matrix(freq[,setdiff(colnames(freq),"Probability"), with = F][, rsIDs, with = F]), freq$Probability))
}

ld_mats <- mclapply(ld_mat_rsIDs, compute_ld_matrix, freq = freq_dat, mc.cores = args$no_of_threads)

names(ld_mats) <- as.character(sort(unique(leg_dat$block)))

save(ld_mats, file = args$output_file)
