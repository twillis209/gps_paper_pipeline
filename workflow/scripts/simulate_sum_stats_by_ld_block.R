library(data.table)
library(simGWAS)
library(argparse)
library(parallel)

simulated_z_score_par <- function(exp_z_score, ld_mat, nrep=1, ncores = 1){
    sim_z_score <- mvnfast::rmvn(n = nrep, mu = exp_z_score, sigma = ld_mat, ncores 
= ncores)
    if(nrep==1)
        return(c(sim_z_score))
    sim_z_score
}

parser <- ArgumentParser(description = 'Uses simGWAS to simulate GWAS summary statistics for a single LD block.')
parser$add_argument('--hap_file', type = 'character', help = 'Path to haplotype file')
parser$add_argument('--leg_file', type = 'character', help = 'Path to legend file')
parser$add_argument('--bim_file', type = 'character', help = 'Path to bim file')
parser$add_argument('--block_no', type = 'integer', help = 'Block number')
parser$add_argument('-b', '--block_file', type = 'character', help = 'Path to block file')
parser$add_argument('--ld_mat_file', type = 'character', help = 'Path to LD matrix file')
parser$add_argument('--chr_no', type = 'integer', help = 'Number of chromosome')
parser$add_argument('--causal_variant_ind', type = 'integer', nargs = '*', help = 'Indices of causal variants within LD block.')
parser$add_argument('--effect_size', type = 'character', help = 'Determines odds ratios for causal variants', required = T)
parser$add_argument('--no_controls', type = 'integer', help = 'No. of controls', default = 10000)
parser$add_argument('--no_cases', type = 'integer', help = 'No. of cases', default = 10000)
parser$add_argument('--no_reps', type = 'integer', help = 'No. of replicates', default = 1)
parser$add_argument('-o', '--output_file', type = 'character', help = 'Path to output file', required = T)
parser$add_argument('-nt', '--no_of_threads', type = 'integer', help = 'Number of threads to use', default = 1)

test_args <- c('--hap_file', 'resources/simgwas/1000g/1000GP_Phase3_chr1_with_meta_eur_common_maf.hap.gz', '--leg_file', 'resources/simgwas/1000g/1000GP_Phase3_chr1_eur_common_maf.legend.gz', '--bim_file', 'resources/1000g/chr1.bim', '-b', 'resources/ldetect/blocks.txt', '--block_no', 40, '--ld_mat_file', 'results/simgwas/chr1_ld_matrices/chr1_block_40_ld_matrix.RData', '--chr_no', 1, '--causal_variant_ind', 2000, '--effect_size', 'null', '--output_file', 'chr1_block_40.tsv.gz', '--no_of_threads', 8, '--no_reps', 1)
args <- parser$parse_args(test_args)

args <- parser$parse_args()

setDTthreads(args$no_of_threads)

leg_dat <- fread(file = args$leg_file, sep = ' ', header = T)
# Skip leading metadata rows
hap_dat <- fread(file = args$hap_file, sep = ' ', header = F, skip = 4)
bim_dat <- fread(file = args$bim_file, sep = '\t', header = F, col.names = c('chr', 'rsID', 'Cm', 'bp', 'A1', 'A2'))
block_dat <- fread(file = args$block_file, sep = ' ', header = F, col.names = c('block', 'chr', 'start', 'stop'))
load(file = args$ld_mat_file)

if(is.null(args$causal_variant_ind)) {
  args$causal_variant_ind <- sample(1:ncol(ld_mat), size = 1)
}

odds_ratios <- list('null' = 1,
                    'small' = 1.05,
                    'medium' = 1.2,
                    'large' = 1.4)

if(!(args$effect_size %in% names(odds_ratios))) {
  stop(sprintf("Effect size %s must be one of the following: '%s.", args$effect_size, paste(names(odds_ratios), collapse = ',')))
}

args$odds_ratios <- sample(odds_ratios[[args$effect_size]], size = length(args$causal_variant_ind), replace = T)

if(args$causal_variant_ind > ncol(ld_mat)) {
  args$causal_variant_ind <- ncol(ld_mat)/2
}

block_dat <- block_dat[chr == args$chr_no]

leg_dat[, rs := make.names(id)]

hap_dat[, rs := leg_dat$rs]

for(i in 1:nrow(block_dat)) {
   leg_dat[position %between% c(block_dat[i, start], block_dat[i, stop]), block := (i-1)]
}

leg_dat <- leg_dat[block == args$block_no]

hap_dat <- hap_dat[rs %in% leg_dat$rs]

hap_dat[, block := leg_dat$block]

hap_dat <- hap_dat[, c(seq(1, ncol(hap_dat)-3, by = 2), seq(2, ncol(hap_dat)-2, by = 2), ncol(hap_dat)-1, ncol(hap_dat)), with = F]

for(j in 1:(ncol(hap_dat)-2)) {
  set(hap_dat, j = j, value = as.numeric(hap_dat[[j]]))
}

hap_mat <- as.matrix(hap_dat[,1:(ncol(hap_dat)-2)])

hap_meta_dat <- hap_dat[, (ncol(hap_dat)-1):ncol(hap_dat)]

rm(hap_dat)

freq_dat <- data.table(t(hap_mat)+1)

rm(hap_mat)

colnames(freq_dat) <- hap_meta_dat$rs

# Required by the make_GenoProbList function below (at least)
freq_dat[, Probability := 1/.N]

chosen_snps <- colnames(ld_mat)

cv_ind <- args$causal_variant_ind

cv_snp <- chosen_snps[cv_ind]

sub_freq_dat <- freq_dat[, c(colnames(ld_mat), 'Probability'), with = F]

geno_probs <- make_GenoProbList(snps = colnames(ld_mat),
                                W = colnames(ld_mat)[cv_ind],
                                freq = sub_freq_dat)

zexp <- expected_z_score(N0 = args$no_controls,
                          N1 = args$no_cases,
                          snps = chosen_snps,
                          W = cv_snp,
                          gamma.W = log(args$odds_ratios),
                          freq = sub_freq_dat,
                          GenoProbList = geno_probs)

zsim <- simulated_z_score_par(exp_z_score = zexp, ld_mat = ld_mat, nrep = args$no_reps, ncores = args$no_of_threads)

vbetasim <- simulated_vbeta(N0 = args$no_controls,
                            N1 = args$no_cases,
                            snps = chosen_snps,
                            W = cv_snp,
                            gamma.W = log(args$odds_ratios),
                            freq = sub_freq_dat,
                            nrep = args$no_reps)

betasim <- zsim * sqrt(vbetasim)

if(args$no_reps == 1 ) {
  result_dat <- data.table(leg_dat[rs %in% chosen_snps, .(id, position, a0, a1, TYPE, EUR)], zexp, zsim, vbetasim = t(vbetasim), betasim = t(betasim))
} else {
  # Add p-values for multiple reps
  result_dat <- data.table(leg_dat[rs %in% chosen_snps, .(id, position, a0, a1, TYPE, EUR)], zexp, t(zsim), t(vbetasim), t(betasim))
}

for(j in 1:args$no_reps) {
  result_dat[, c(paste0('p.', j)) := 2*pnorm(abs(result_dat[[7+j]]), lower.tail = F)]

  result_dat[, block := args$block_no]
  result_dat[, or := 0]
  result_dat[cv_ind, or := args$odds_ratios]
}

names(result_dat) <- c('id', 'position', 'a0', 'a1', 'TYPE', 'EUR', 'zexp',
                    paste0('zsim.', 1:args$no_reps),
                    paste0('vbetasim.', 1:args$no_reps),
                    paste0('betasim.', 1:args$no_reps),
                    paste0('p.', 1:args$no_reps),
                    'block',
                    'chosen_or')

result_dat[, `:=` (ncases = args$no_cases, ncontrols = args$no_controls)]

result_dat <- merge(result_dat, bim_dat[, .(rsID, bp, A1, A2)], by.x = 'position', by.y = 'bp', all.x = T)

result_dat <- result_dat[(a0 == A2 & a1 == A1) | (a0 == A1 & a1 == A2)]

result_dat[, c("A1", "A2", "id") := NULL]

result_dat[, chr := args$chr_no]

fwrite(result_dat, file = args$output_file, sep = '\t')
