library(data.table)
library(simGWAS)
library(tictoc)
library(parallel)

no_snps <- 20000
no_causal_variants <- 1
odds_ratios <- 1.3
output_file <- 'benchmarking_test.tsv'
no_of_threads <- 8
no_cases <- no_controls <- 10000
no_reps <- 1

setDTthreads(no_of_threads)

leg_dat <- fread(file = 'resources/simgwas/1000g/1000GP_Phase3_chr21.legend.gz', sep = ' ', header = T)
hap_dat <- fread(file = 'resources/simgwas/1000g/1000GP_Phase3_chr21_with_meta_eur.hap.gz', sep = ' ', header = F)

# Drop metadata rows
hap_dat <- hap_dat[5:nrow(hap_dat)]

leg_dat <- leg_dat[1:no_snps]
hap_dat <- hap_dat[1:no_snps]

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

causal_variant_ind <- sample(1:nrow(chosen_snp_dat[EUR %between% c(0.2, 0.8)]), size = no_causal_variants)

# Genotype probabilities at each SNP conditional on genotypes at causal variants
geno_probs <- make_GenoProbList(snps = chosen_snp_dat$rs,
                                W = chosen_snp_dat$rs[causal_variant_ind],
                                freq = freq_dat)

zexp <- expected_z_score(N0 = no_controls, # number of controls
                    N1 = no_cases, # number of cases
                    snps = chosen_snp_dat$rs, # column names in freq of SNPs for which Z scores should be generated
                    W = chosen_snp_dat$rs[causal_variant_ind], # causal variants, subset of snps
                    gamma.W = log(odds_ratios), # odds ratios
                    freq = freq_dat, # reference haplotypes
                    GenoProbList = geno_probs)

# NB: Taken from the simGWAS intro vignette; some of the functions have changed, I know Chris has developed simGWAS since the paper
tic()
zsim <- simulated_z_score(N0 = no_controls,
                          N1 = no_cases,
                          snps = chosen_snp_dat$rs,
                          W = chosen_snp_dat$rs[causal_variant_ind],
                          gamma.W = log(odds_ratios),
                          freq = freq_dat,
                          nrep = no_reps)
toc()


#tic()
#zsim <- simulated_z_score_par(N0 = no_controls,
#                          N1 = no_cases,
#                          snps = chosen_snp_dat$rs,
#                          W = chosen_snp_dat$rs[causal_variant_ind],
#                          gamma.W = log(odds_ratios),
#                          freq = freq_dat,
#                          nrep = no_reps,
#                          ncores = no_of_threads)
#toc()

vbetasim <- simulated_vbeta(N0= no_controls,
                            N1= no_cases,
                            snps = chosen_snp_dat$rs,
                            W = chosen_snp_dat$rs[causal_variant_ind],
                            gamma.W = log(odds_ratios),
                            freq = freq_dat,
                            nrep = no_reps)

betasim <- zsim * sqrt(vbetasim)

if(no_reps == 1 ) {
  #res_dat <- data.table(chosen_snp_dat[, .(id, position, a0, a1, TYPE, EUR)], zexp, zsim, vbetasim, betasim)
  res_dat <- data.table(chosen_snp_dat[, .(id, position, a0, a1, TYPE, EUR)], zexp, zsim, vbetasim = t(vbetasim), betasim = t(betasim))
  res_dat[, p := 2*pnorm(abs(zsim), lower.tail = F)]
} else {
  res_dat <- data.table(chosen_snp_dat[, .(id, position, a0, a1, TYPE, EUR)], zexp, t(zsim), t(vbetasim), t(betasim))
}

res_dat[causal_variant_ind, chosen_or := odds_ratios]

names(res_dat) <- c('id', 'position', 'a0', 'a1', 'TYPE', 'EUR', 'zexp',
                    paste0('zsim.', 1:no_reps),
                    paste0('vbetasim.', 1:no_reps),
                    paste0('betasim.', 1:no_reps), 'chosen_or')

res_dat[, `:=` (ncases = no_cases, ncontrols = no_controls)]

fwrite(res_dat, file = output_file, sep = '\t')
