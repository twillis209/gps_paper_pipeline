library(data.table)
library(simGWAS)
library(argparse)
library(parallel)

parser <- ArgumentParser(description = 'Computes blockwise LD matrix for specified chromosome and LD block.')
parser$add_argument('--hap_file', type = 'character', help = 'Path to haplotype file')
parser$add_argument('--leg_file', type = 'character', help = 'Path to legend file')
parser$add_argument('--block_file', type = 'character', help = 'Path to block file')
parser$add_argument('--block_no', type = 'integer', help = 'LD block number')
parser$add_argument('--chr_no', type = 'integer', help = 'Number of chromosome')
parser$add_argument('--output_file', type = 'character', help = 'Path to output file', required = T)
parser$add_argument('-nt', '--no_of_threads', type = 'integer', help = 'Number of threads to use', default = 1)

test_args <- c('--hap_file', "resources/simgwas/1000g/1000GP_Phase3_chr1_with_meta_eur.hap.gz", '--leg_file', "resources/simgwas/1000g/1000GP_Phase3_chr1.legend.gz", '--output_file',            "results/simgwas/chr1_ld_matrices/chr1_block_21_ld_matrix.RData", '-nt', 4, '--block_file', 'resources/ldetect/blocks.txt', '--chr_no', 1, '--block_no', 21)
args <- parser$parse_args(test_args)

args <- parser$parse_args()

setDTthreads(args$no_of_threads)

leg_dat <- fread(file = args$leg_file, sep = ' ', header = T)
# Skip leading metadata rows
hap_dat <- fread(file = args$hap_file, sep = ' ', header = F, skip = 4)
block_dat <- fread(file = args$block_file, sep = ' ', header = F, col.names = c('block', 'chr', 'start', 'stop'))

block_dat <- block_dat[chr == args$chr_no]

leg_dat[, rs := make.names(leg_dat$id)]

hap_dat[, rs := leg_dat$rs]

for(i in 1:nrow(block_dat)) {
   leg_dat[position %between% c(block_dat[i, start], block_dat[i, stop]), block := (i-1)]
}

# Drop SNPs we can't assign to a block
leg_dat <- leg_dat[!is.na(block)]

hap_dat <- hap_dat[rs %in% leg_dat$rs]

hap_dat <- hap_dat[block == args$block_no]

leg_dat <- leg_dat[block == args$block_no]

hap_dat <- hap_dat[, c(seq(1, ncol(hap_dat)-3, by = 2), seq(2, ncol(hap_dat)-2, by = 2), ncol(hap_dat)-1, ncol(hap_dat)), with = F]

for(j in 1:(ncol(hap_dat)-2)) {
  set(hap_dat, j = j, value = as.numeric(hap_dat[[j]]))
}

hap_mat <- as.matrix(hap_dat[,1:(ncol(hap_dat)-2)])

freq_dat <- data.table(t(hap_mat)+1)

colnames(freq_dat) <- hap_dat$rs

freq_dat[, Probability := 1/.N]

rm(hap_dat, hap_mat)

ld_mat <- corpcor::make.positive.definite(simGWAS:::wcor2(as.matrix(freq_dat[,setdiff(colnames(freq_dat),"Probability"), with = F][, leg_dat$rs, with = F]), freq_dat$Probability))

save(ld_mat, file = args$output_file)
