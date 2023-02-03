library(data.table)
library(simGWAS)
library(magrittr)

simulated_z_score_par <- function(exp_z_score, ld_mat, nrep=1, ncores = 1){
    sim_z_score <- mvnfast::rmvn(n = nrep, mu = exp_z_score, sigma = ld_mat, ncores 
= ncores)
    if(nrep==1)
        return(c(sim_z_score))
    sim_z_score
}

no_reps <- snakemake@params[['no_reps']]

setDTthreads(snakemake@threads)

set.seed(snakemake@wildcards[['seed']])

leg_dat <- fread(file = snakemake@input[['block_legend_file']], sep = ' ', header = T)
hap_dat <- fread(file = snakemake@input[['block_haplotype_file']], sep = ' ', header = F)
bim_dat <- fread(file = snakemake@input[['bim_file']], sep = '\t', header = F, col.names = c('chr', 'bim.id', 'Cm', 'bp', 'A1', 'A2'))
load(file = snakemake@input[['ld_mat_file']])

log_odds_ratio <- log(snakemake@params[['odds_ratio']])

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

cv_ind <- sample(1:length(chosen_snps), size = snakemake@wildcards[['no_cvariants']])

cv_snps <- chosen_snps[cv_ind]

sub_freq_dat <- freq_dat[, c(colnames(ld_mat), 'Probability'), with = F]

geno_probs <- make_GenoProbList(snps = colnames(ld_mat),
                                W = colnames(ld_mat)[cv_ind],
                                freq = sub_freq_dat)

zexp <- expected_z_score(N0 = snakemake@params[['no_controls']],
                          N1 = snakemake@params[['no_cases']],
                          snps = chosen_snps,
                          W = cv_snps,
                          gamma.W = rep(log_odds_ratio, length(cv_snps)),
                          freq = sub_freq_dat,
                          GenoProbList = geno_probs)

zsim <- simulated_z_score_par(exp_z_score = zexp, ld_mat = ld_mat, nrep = no_reps, ncores = snakemake@threads)

# Both vbetasim and betasim overflow if passed as integers
vbetasim <- simulated_vbeta(N0 = as.numeric(snakemake@params[['no_controls']]),
                            N1 = as.numeric(snakemake@params[['no_cases']]),
                            snps = chosen_snps,
                            W = cv_snps,
                            gamma.W = rep(log_odds_ratio, length(cv_snps)),
                            freq = sub_freq_dat,
                            nrep = no_reps)

betasim <- zsim * sqrt(vbetasim)

if(no_reps == 1 ) {
  result_dat <- data.table(leg_dat[rs %in% chosen_snps, .(id, position, block, a0, a1, TYPE, EUR)], zexp, zsim, vbetasim = t(vbetasim), betasim = t(betasim))
} else {
  # Add p-values for multiple reps
  result_dat <- data.table(leg_dat[rs %in% chosen_snps, .(id, position, block, a0, a1, TYPE, EUR)], zexp, t(zsim), t(vbetasim), t(betasim))
}

names(result_dat) <- c('id', 'position', 'block', 'a0', 'a1', 'TYPE', 'EUR', 'zexp',
                       paste0('zsim.', 1:no_reps),
                       paste0('vbetasim.', 1:no_reps),
                       paste0('betasim.', 1:no_reps))

for(j in 1:no_reps) {
  result_dat[, c(paste0('p.', j)) := 2*pnorm(abs(get(paste0('zsim.', j))), lower.tail = F)]
}

result_dat[, or := 1]
result_dat[cv_ind, or := exp(log_odds_ratio)]
setnames(result_dat, 'or', 'chosen_or')

result_dat[, `:=` (ncases = snakemake@params[['no_cases']], ncontrols = snakemake@params[['no_controls']])]

# This merging step confines our SNPs to those present in the bim file; we keep the same order of a0 and a1 alleles, even if this differs in the bim file
result_dat <- merge(result_dat, bim_dat[, .(bim.id, bp, A1, A2)], by.x = 'position', by.y = 'bp', all.x = T)

result_dat <- result_dat[(a0 == A2 & a1 == A1) | (a0 == A1 & a1 == A2)]

result_dat[, c("A1", "A2", "bim.id") := NULL]

result_dat[, chr := snakemake@params[['chr_no']]]

result_dat[, id := tstrsplit(id, split = ':')[[1]]]

# NB: We add 'block_effect_size' to the end of this
cols <- c("position", "a0", "a1", "id", "block", "TYPE", "EUR", "zexp", paste0("zsim.", 1:no_reps), paste0("vbetasim.", 1:no_reps), paste0("betasim.", 1:no_reps), paste0("p.", 1:no_reps), "chosen_or", "ncases", "ncontrols", "chr")

result_dat <- result_dat[, ..cols]

result_dat[, block_effect_size := snakemake@wildcards[['effect_size']]]

fwrite(result_dat, file = snakemake@output[[1]], sep = '\t', col.names = F)
