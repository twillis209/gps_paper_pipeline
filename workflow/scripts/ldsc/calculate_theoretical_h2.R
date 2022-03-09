library(data.table)
library(argparse)

obs_liab_trans <- function(h2.obs, P, K) {
  z_2 <- dnorm(qnorm(1-K))^2

  h2.obs * ((K*(1-K))^2/(P*(1-P)))/z_2
}

parser <- ArgumentParser(description = 'Calculate theoretical h2 for simulation from effect sizes and MAFs')
parser$add_argument('--sum_stats_file', type = 'character', help = 'Path to summary statistics file')
parser$add_argument('-P', type = 'double', help = 'Sample case proportion')
parser$add_argument('-K', type = 'double', help = 'Population case prevalence')
parser$add_argument('-o', '--output_path', type = 'character', help = 'Path to pruned summary statistics file', required = T)
parser$add_argument('-nt', '--no_of_threads', type = 'integer', help = 'Number of threads to use', default = 1)

args <- parser$parse_args()

setDTthreads(args$no_of_threads)

sum_stats_dat <- fread(args$sum_stats_file, sep = '\t', header = T, select = c("chosen_or", "EUR"))

# NB: I've been inconsistent in the past with specifying chosen_or as either 0 or 1, so this handles both cases
cv_dat <- sum_stats_dat[chosen_or > 1]

cv_dat[, geno_var := 2*EUR*(1-EUR)]

cv_dat[, beta_2 := log(chosen_or)^2]

cv_dat[, h2.theo.obs := cumsum(beta_2*geno_var)/0.25]
cv_dat[, h2.theo.liab:= obs_liab_trans(h2.theo.obs, P = args$P, K = args$K)]

cv_dat[, beta_2.tiny := log(1.02)^2]
cv_dat[, beta_2.small := log(1.05)^2]
cv_dat[, beta_2.medium := log(1.2)^2]
cv_dat[, beta_2.large := log(1.4)^2]
cv_dat[, beta_2.huge := log(2)^2]

cv_dat[, h2.theo.obs.tiny := cumsum(beta_2.tiny*geno_var)/0.25]
cv_dat[, h2.theo.obs.small := cumsum(beta_2.small*geno_var)/0.25]
cv_dat[, h2.theo.obs.medium := cumsum(beta_2.medium*geno_var)/0.25]
cv_dat[, h2.theo.obs.large := cumsum(beta_2.large*geno_var)/0.25]
cv_dat[, h2.theo.obs.huge := cumsum(beta_2.huge*geno_var)/0.25]

cv_dat[, h2.theo.liab.tiny := obs_liab_trans(h2.theo.obs.tiny, P = args$P, K = args$K)]
cv_dat[, h2.theo.liab.small := obs_liab_trans(h2.theo.obs.small, P = args$P, K = args$K)]
cv_dat[, h2.theo.liab.medium := obs_liab_trans(h2.theo.obs.medium, P = args$P, K = args$K)]
cv_dat[, h2.theo.liab.large := obs_liab_trans(h2.theo.obs.large, P = args$P, K = args$K)]
cv_dat[, h2.theo.liab.huge := obs_liab_trans(h2.theo.obs.huge, P = args$P, K = args$K)]

fwrite(cv_dat, file = args$output_path, sep = '\t')
