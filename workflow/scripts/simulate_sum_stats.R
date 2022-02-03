library(data.table)
library(simGWAS)
library(argparse)
library(parallel)

parser <- ArgumentParser(description = 'Uses simGWAS to simulate GWAS summary statistics for SNPs on chromosome 21')
parser$add_argument('--hap_file', type = 'character', help = 'Path to haplotype file')
parser$add_argument('--leg_file', type = 'character', help = 'Path to legend file')
parser$add_argument('--bim_file', type = 'character', help = 'Path to bim file')
parser$add_argument('--ld_mats_file', type = 'character', help = 'Path to file containing LD matrices')
parser$add_argument('-b', '--block_file', type = 'character', help = 'Path to block file')
parser$add_argument('--no_blocks', type = 'integer', help = 'Number of LD blocks on chromosome 21 to simulate')
parser$add_argument('--causal_variant_ind', type = 'integer', nargs = '*', help = 'Indices of causal variants with each LD block.')
parser$add_argument('--odds_ratios', type = 'double', nargs = '+', help = 'Odds ratios for specified causal variants', default = c(1.1, 1.3, 1.5))
parser$add_argument('--no_controls', type = 'integer', help = 'No. of controls', default = 1e4)
parser$add_argument('--no_cases', type = 'integer', help = 'No. of cases', default = 1e4)
parser$add_argument('--no_reps', type = 'integer', help = 'No. of replicates', default = 1)
parser$add_argument('-o', '--output_file', type = 'character', help = 'Path to output file', required = T)
parser$add_argument('-nt', '--no_of_threads', type = 'integer', help = 'Number of threads to use', default = 1)

#test_args <- c('--hap_file', 'resources/simgwas/1000g/1000GP_Phase3_chr21_with_meta_eur.hap.gz', '--leg_file', 'resources/simgwas/1000g/1000GP_Phase3_chr21.legend.gz', '--bim_file', 'resources/1000g/chr21.bim', '--block_file', 'resources/ldetect/blocks.txt', '--ld_mats_file', 'results/simgwas/chr21_block_ld_matrices.RData', '--no_blocks', 1, '--causal_variant_ind', 2000, 2000, '--odds_ratios', 3, '--output_file', 'test.tsv.gz', '--no_of_threads', 8, '--no_reps', 2)
#args <- parser$parse_args(test_args)

args <- parser$parse_args()

setDTthreads(args$no_of_threads)

leg_dat <- fread(file = args$leg_file, sep = ' ', header = T)
hap_dat <- fread(file = args$hap_file, sep = ' ', header = F)
bim_dat <- fread(file = args$bim_file, sep = '\t', header = F, col.names = c('chr', 'rsID', 'Cm', 'bp', 'A1', 'A2'))
block_dat <- fread(file = args$block_file, sep = ' ', header = F, col.names = c('block', 'chr', 'start', 'stop'))
load(file = args$ld_mats_file)

if(!is.null(args$causal_variant_ind)) {
  args$causal_variant_ind <- integer()

  for(i in 1:args$no_blocks) {
    args$causal_variant_ind[i] <- sample(1:ncol(ld_mats[[i]]), size = 1)
  }
} else if(length(args$causal_variant_ind) < length(args$no_blocks)) {

  for(i in (length(args$causal_variant_ind)+1):args$no_blocks) {
    args$causal_variant_ind[i] <- sample(1:ncol(ld_mats[[i]]), size = 1)
  }
}

if(length(args$odds_ratios) < length(args$causal_variant_ind)) {
  args$odds_ratios <- c(args$odds_ratios, sample(1:ncol(ld_mats[[i]]), size = length(args$causal_variant_ind)-length(args$odds_ratios)))
}

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

hap_meta_dat <- hap_dat[, (ncol(hap_dat)-1):ncol(hap_dat)]

rm(hap_dat)

freq_dat <- data.table(t(hap_mat)+1)

rm(hap_mat)

colnames(freq_dat) <- hap_meta_dat$rs

# Required by the make_GenoProbList function below (at least)
freq_dat[, Probability := 1/.N]

ld_mats <- ld_mats[1:args$no_blocks]

result_dat_list <- mclapply(seq_along(ld_mats), function(i) {
  ld_mat <- ld_mats[[i]]

  chosen_snps <- colnames(ld_mat)

  cv_ind <- args$causal_variant_ind[i]

  cv_snp <- chosen_snps[cv_ind]

  or <- args$odds_ratios[i]

  sub_freq_dat <- freq_dat[, c(colnames(ld_mat), 'Probability'), with = F]

  geno_probs <- make_GenoProbList(snps = colnames(ld_mat),
                                  W = colnames(ld_mat)[cv_ind],
                                  freq = sub_freq_dat)

  zexp <- expected_z_score(N0 = args$no_controls,
                           N1 = args$no_cases,
                           snps = chosen_snps,
                           W = cv_snp,
                           gamma.W = log(or),
                           freq = sub_freq_dat,
                           GenoProbList = geno_probs)

  zsim <- simulated_z_score_par(exp_z_score = zexp, ld_mat = ld_mat, nrep = args$no_reps, ncores = min(args$no_reps, args$no_of_threads %/% args$no_blocks))

  vbetasim <- simulated_vbeta(N0 = args$no_controls,
                              N1 = args$no_cases,
                              snps = chosen_snps,
                              W = cv_snp,
                              gamma.W = log(or),
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
  }

  result_dat[, block := names(ld_mats)[[i]]]
  result_dat[id == cv_snp, or := or]

  result_dat
}, mc.cores = args$no_blocks)

res_dat <- rbindlist(result_dat_list)

names(res_dat) <- c('id', 'position', 'a0', 'a1', 'TYPE', 'EUR', 'zexp',
                    paste0('zsim.', 1:args$no_reps),
                    paste0('vbetasim.', 1:args$no_reps),
                    paste0('betasim.', 1:args$no_reps),
                    paste0('p.', 1:args$no_reps),
                    'block',
                    'chosen_or')

res_dat[, `:=` (ncases = args$no_cases, ncontrols = args$no_controls)]

res_dat <- merge(res_dat, bim_dat[, .(rsID, bp, A1, A2)], by.x = 'position', by.y = 'bp', all.x = T)

res_dat <- res_dat[(a0 == A2 & a1 == A1) | (a0 == A1 & a1 == A2)]

res_dat[, c("A1", "A2", "id") := NULL]

res_dat[, chr := 21]

fwrite(res_dat, file = args$output_file, sep = '\t')
