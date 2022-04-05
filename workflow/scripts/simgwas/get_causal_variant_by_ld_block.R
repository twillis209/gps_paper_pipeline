library(data.table)
library(argparse)

parser <- ArgumentParser(description = 'Get causal variant in LD block.')
parser$add_argument('--hap_file', type = 'character', help = 'Path to haplotype file')
parser$add_argument('--leg_file', type = 'character', help = 'Path to legend file')
parser$add_argument('--bim_file', type = 'character', help = 'Path to bim file')
parser$add_argument('--ld_mat_file', type = 'character', help = 'Path to LD matrix file')
parser$add_argument('--chr_no', type = 'integer', help = 'Number of chromosome')
parser$add_argument('--causal_variant_ind', type = 'integer', nargs = '*', help = 'Indices of causal variants within LD block.')
parser$add_argument('-o', '--output_file', type = 'character', help = 'Path to output file', required = T)
parser$add_argument('-nt', '--no_of_threads', type = 'integer', help = 'Number of threads to use', default = 1)

#test_args <- c('--hap_file', 'resources/simgwas/1000g/blockwise/chr1/block_0.hap.gz', '--leg_file', 'resources/simgwas/1000g/blockwise/chr1/block_0.legend.gz', '--bim_file', 'resources/1000g/chr1.bim', '--ld_mat_file', 'results/simgwas/chr1_ld_matrices/chr1_block_0_ld_matrix.RData', '--chr_no', 1, '--causal_variant_ind', 2000, '--effect_size', 'null', '--output_file', 'chr1_block_0_sum_stats.tsv.gz', '--no_of_threads', 1, '--no_reps', 2, '--no_controls', 40000, '--no_cases', 40000)
#args <- parser$parse_args(test_args)

args <- parser$parse_args()

setDTthreads(args$no_of_threads)

leg_dat <- fread(file = args$leg_file, sep = ' ', header = T)
hap_dat <- fread(file = args$hap_file, sep = ' ', header = F)
bim_dat <- fread(file = args$bim_file, sep = '\t', header = F, col.names = c('chr', 'bim.id', 'Cm', 'bp', 'A1', 'A2'))
load(file = args$ld_mat_file)

if(is.null(args$causal_variant_ind)) {
  args$causal_variant_ind <- sample(1:ncol(ld_mat), size = 1)
}

if(args$causal_variant_ind > ncol(ld_mat)) {
  args$causal_variant_ind <- ncol(ld_mat)/2
}

for(j in 1:(ncol(hap_dat)-2)) {
  set(hap_dat, j = j, value = as.numeric(hap_dat[[j]]))
}

hap_mat <- as.matrix(hap_dat[,1:(ncol(hap_dat)-2)])

hap_meta_dat <- hap_dat[, (ncol(hap_dat)-1):ncol(hap_dat)]

names(hap_meta_dat) <- c('rs', 'block')

rm(hap_dat)

freq_dat <- data.table(t(hap_mat)+1)

rm(hap_mat)

colnames(freq_dat) <- hap_meta_dat$rs

# Required by the make_GenoProbList function below (at least)
freq_dat[, Probability := 1/.N]

chosen_snps <- colnames(ld_mat)

cv_ind <- args$causal_variant_ind

cv_snp <- chosen_snps[cv_ind]

result_dat <- data.table(leg_dat[rs == cv_snp, .(id, position, block, a0, a1, TYPE, EUR)])

result_dat <- merge(result_dat, bim_dat[, .(bim.id, bp, A1, A2)], by.x = 'position', by.y = 'bp', all.x = T)

result_dat <- result_dat[(a0 == A2 & a1 == A1) | (a0 == A1 & a1 == A2)]

result_dat[, c("A1", "A2", "bim.id") := NULL]

result_dat[, chr := args$chr_no]

result_dat[, rsID := tstrsplit(id, split = ':')[[1]]]

fwrite(result_dat, file = args$output_file, sep = '\t')
