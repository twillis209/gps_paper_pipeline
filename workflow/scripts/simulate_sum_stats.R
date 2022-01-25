library(data.table)
library(simGWAS)
library(argparse)

parser <- ArgumentParser(description = 'Computes p-value for GPS statistic using permuted data')
parser$add_argument('-a', '--hap_file', type = 'character', help = 'Path to haplotype file')
parser$add_argument('-l', '--leg_file', type = 'character', help = 'Path to legend file')
parser$add_argument('-s', '--no_snps', type = 'integer', help = 'Number of SNPs on chromosome 21 to simulate', default = 2e4)
parser$add_argument('-c', '--no_causal_variants', type = 'integer', help = 'Number of causal variants', required = T)
parser$add_argument('-n', '--causal_variant_ind', type = 'integer', nargs = '*', help = 'Indices of causal variants. If omitted, random indices will be selected.', required = F)
parser$add_argument('-r', '--odds_ratios', type = 'double', nargs = '+', help = 'Odds ratios for specified causal variants', required = T)
parser$add_argument('--no_controls', type = 'integer', help = 'No. of controls', default = 1e4)
parser$add_argument('--no_cases', type = 'integer', help = 'No. of cases', default = 1e4)
parser$add_argument('--no_reps', type = 'integer', help = 'No. of replicates', default = 1)
parser$add_argument('-o', '--output_file', type = 'character', help = 'Path to output file', required = T)
parser$add_argument('-nt', '--no_of_threads', type = 'integer', help = 'Number of threads to use', default = 1)

# TODO remove me
test_args <- c('--hap_file', 'resources/simgwas/1000g/1000GP_Phase3_chr21_with_meta_eur.hap.gz', '--leg_file', 'resources/simgwas/1000g/1000GP_Phase3_chr21.legend.gz', '--no_causal_variants', 1, '--odds_ratios', 3, '-o', 'test.tsv', '--no_of_threads', 8)

args <- parser$parse_args()

if(!is.null(args$causal_variant_ind)) {
  if(args$no_causal_variants != length(args$causal_variant_ind)) {
    stop("Must specify indices for all causal variants")
  }
}

if(args$no_causal_variants != length(args$odds_ratios)) {
  stop("Must specify odds ratios for all causal variants")
}

setDTthreads(args$no_of_threads)

leg_dat <- fread(file = args$leg_file, sep = ' ', header = T)
hap_dat <- fread(file = args$hap_file, sep = ' ', header = F)

# Drop metadata rows
hap_dat <- hap_dat[5:nrow(hap_dat)]

leg_dat <- leg_dat[1:args$no_snps]
hap_dat <- hap_dat[1:args$no_snps]

hap_dat <- hap_dat[, c(seq(1, ncol(hap_dat)-1, by = 2), seq(2, ncol(hap_dat), by = 2)), with = F]

hap_dat <- hap_dat[, lapply(.SD, as.numeric)]

hap_mat <- as.matrix(hap_dat)

freq_dat <- data.table(t(hap_mat)+1)

leg_dat[, rs := make.names(leg_dat$id)]
colnames(freq_dat) <- leg_dat$rs

mask <- leg_dat$EUR %between% c(0.01, 0.99) & apply(freq_dat, 2, var) > 0

chosen_snp_dat <- leg_dat[mask]

freq_dat <- freq_dat[, mask, with = F]

chosen_snp_dat[, maf := colMeans(freq_dat-1)]

ld_mat <- cor(freq_dat)
ld_mat <- as.matrix(corpcor::make.positive.definite(ld_mat))

freq_dat[, Probability := 1/.N]

if(is.null(args$causal_variant_ind)) {
  args$causal_variant_ind <- sample(1:nrow(chosen_snp_dat[EUR %between% c(0.2, 0.8)]), size = args$no_causal_variants)
}

# Genotype probabilities at each SNP conditional on genotypes at causal variants
geno_probs <- make_GenoProbList(snps = chosen_snp_dat$rs,
                                W = chosen_snp_dat$rs[args$causal_variant_ind],
                                freq = freq_dat)

zexp <- expected_z_score(N0 = args$no_controls, # number of controls
                    N1 = args$no_cases, # number of cases
                    snps = chosen_snp_dat$rs, # column names in freq of SNPs for which Z scores should be generated
                    W = chosen_snp_dat$rs[args$causal_variant_ind], # causal variants, subset of snps
                    gamma.W = log(args$odds_ratios), # odds ratios
                    freq = freq_dat, # reference haplotypes
                    GenoProbList = geno_probs)

# NB: Taken from the simGWAS intro vignette; some of the functions have changed, I know Chris has developed simGWAS since the paper
zsim <- simulated_z_score(N0 = args$no_controls,
                          N1 = args$no_cases,
                          snps = chosen_snp_dat$rs,
                          W = chosen_snp_dat$rs[args$causal_variant_ind],
                          gamma.W = log(args$odds_ratio),
                          freq = freq_dat,
                          nrep = args$no_reps)

vbetasim <- simulated_vbeta(N0= args$no_controls,
                            N1= args$no_cases,
                            snps = chosen_snp_dat$rs,
                            W = chosen_snp_dat$rs[args$causal_variant_ind],
                            gamma.W = log(args$odds_ratio),
                            freq = freq_dat,
                            nrep = args$no_reps)

betasim <- zsim / sqrt(vbetasim)

if(args$no_reps == 1 ) {
  #res_dat <- data.table(chosen_snp_dat[, .(id, position, a0, a1, TYPE, EUR)], zexp, zsim, vbetasim, betasim)
  res_dat <- data.table(chosen_snp_dat[, .(id, position, a0, a1, TYPE, EUR)], zexp, zsim, vbetasim = t(vbetasim), betasim = t(betasim))
} else {
  res_dat <- data.table(chosen_snp_dat[, .(id, position, a0, a1, TYPE, EUR)], zexp, t(zsim), t(vbetasim), t(betasim))
}

res_dat[args$causal_variant_ind, chosen_or := args$odds_ratios]

names(res_dat) <- c('id', 'position', 'a0', 'a1', 'TYPE', 'EUR', 'zexp',
                    paste0('zsim.', 1:args$no_reps),
                    paste0('vbetasim.', 1:args$no_reps),
                    paste0('betasim.', 1:args$no_reps), 'chosen_or')

fwrite(res_dat, file = args$output_file, sep = '\t')
