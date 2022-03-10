library(data.table)
library(argparse)
library(stringr)
library(magrittr)

obs_liab_trans <- function(h2.obs, P, K) {
  z_2 <- dnorm(qnorm(1-K))^2

  h2.obs * ((K*(1-K))^2/(P*(1-P)))/z_2
}

get_chr_block_for_variants <- function(effect_blocks_string) {
  str_split(effect_blocks_string, "\\+") %>%
    unlist %>%
    str_match(., "(\\d+)-([simlh])(\\d+):(\\d+)") %>%
    apply(., 1, FUN = function(x) paste(as.integer(x[2]), as.integer(x[4]):as.integer(x[5]), sep = '_')) %>% unlist
}

parser <- ArgumentParser(description = 'Calculate theoretical h2 for simulation from effect sizes and MAFs')
parser$add_argument('--cv_file', type = 'character', help = 'Path to summary statistics file containing causal variants')
parser$add_argument('--effect_blocks', type = 'character', help = 'Effect blocks for trait')
parser$add_argument('--odds_ratio', type = 'double', help = 'Odds ratio')
parser$add_argument('-P', type = 'double', help = 'Sample case proportion')
parser$add_argument('-K', type = 'double', help = 'Population case prevalence')
parser$add_argument('-o', '--output_path', type = 'character', help = 'Path to pruned summary statistics file', required = T)
parser$add_argument('-nt', '--no_of_threads', type = 'integer', help = 'Number of threads to use', default = 1)

args <- parser$parse_args()

setDTthreads(args$no_of_threads)

cv_dat <- fread(args$cv_file, sep = '\t', header = T)

cv_dat[, geno_var := 2*EUR*(1-EUR)]

cv_dat[, chr_block := paste(chr, block, sep = '_')]

cv_dat[, in_effect_blocks := chr_block %in% get_chr_block_for_variants(args$effect_blocks)]

cv_dat[in_effect_blocks == T, beta := log(args$odds_ratio)]

cv_dat[in_effect_blocks == T, beta_2 := beta^2]

V_A <- with(cv_dat[in_effect_blocks == T], sum(beta_2*geno_var))

h2.theo.obs <- V_A/(args$P*(1-args$P))

h2.theo.liab <- obs_liab_trans(h2.theo.obs, P = args$P, K = args$K)

res_dat <- data.table(odds_ratio = args$odds_ratio,
                      effect_blocks = args$effect_blocks,
                      no_blocks = cv_dat[in_effect_blocks == T, .N],
                      h2.theo.obs,
                      h2.theo.liab)

fwrite(res_dat, file = args$output_path, sep = '\t', na = 'NA')
