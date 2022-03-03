library(data.table)

dat <- fread('chr1_s0:119_maf.tsv', sep = '\t', header = T)

dat[, geno_var := 2*EUR*(1-EUR)]

odds_ratios <- c(1.05, 1.2, 1.4, 2)

log_odds_ratios <- log(odds_ratios)

beta_2 <- log_odds_ratios^2

names(beta_2) <- c('small', 'medium', 'large', 'huge')

no_blocks <- c(10, 20, 40, 80, 100, 120)

h2_dat <- data.table(odds_ratio = numeric(), no_blocks = integer(), h2.theo.obs = numeric())

obs_liab_trans <- function(h2.obs, P, K) {
  z_2 <- dnorm(qnorm(1-K))^2

  h2.obs * ((K*(1-K))^2/(P*(1-P)))/z_2
}

for(i in seq_along(beta_2)) {
  for(j in seq_along(no_blocks)) {
    h2_dat <- rbind(h2_dat,  data.table(odds_ratio = odds_ratios[i], no_blocks = no_blocks[j], h2.theo.obs = sum(dat[1:no_blocks[j], geno_var]*beta_2[i])/0.25))
  }
}

h2_dat[, h2.theo.lia := obs_liab_trans(h2.theo.obs, P = 0.5, K = 0.02)]

chr1_sweep_dat <- fread('results/ldsc/h2/compiled/chr1_sweep.tsv', sep = '\t', header = T)
wg_sweep_dat <- fread('results/ldsc/h2/compiled/whole_genome_sweep.tsv', sep = '\t', header = T)

chr1_sweep_dat[, no_blocks := as.integer(tstrsplit(blocks, '-')[[2]])+1]
wg_sweep_dat[, no_blocks := as.integer(tstrsplit(blocks, '-')[[2]])+1]
chr1_sweep_dat[, dataset := 'chr1']
wg_sweep_dat[, dataset := 'wg']

setnames(chr1_sweep_dat, c('h2', 'tag'), c('h2.est.lia', 'replicate'))
setnames(wg_sweep_dat, c('h2', 'tag'), c('h2.est.lia', 'replicate'))

merged_dat <- rbind(merge(h2_dat, chr1_sweep_dat[, .(ncases, ncontrols, h2.est.lia, se, odds_ratio, no_blocks, dataset, replicate)], by = c('odds_ratio', 'no_blocks')),
                    merge(h2_dat, wg_sweep_dat[, .(ncases, ncontrols, h2.est.lia, se, odds_ratio, no_blocks, dataset, replicate)], by = c('odds_ratio', 'no_blocks')))

merged_dat <- merged_dat[, c('dataset', 'odds_ratio', 'no_blocks', 'ncases', 'ncontrols', 'replicate', 'h2.theo.obs', 'h2.theo.lia', 'h2.est.lia', 'se')]

fwrite(merged_dat, file = 'merged_theo_est_h2.tsv', sep = '\t')
