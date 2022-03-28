library(data.table)
library(argparse)
library(stringr)
library(magrittr)

obs_liab_trans <- function(h2.obs, P, K) {
  z_2 <- dnorm(qnorm(1-K))^2

  h2.obs * ((K*(1-K))^2/(P*(1-P)))/z_2
}

parser <- ArgumentParser(description = 'Calculate theoretical rg')
parser$add_argument('--cv_file', type = 'character', help = 'Path to summary statistics file containing causal variants')
parser$add_argument('--blocks_file', type = 'character', help = 'Path to available blocks file')
parser$add_argument('--effect_blocks_a', type = 'character', nargs = '+', help = 'Effect blocks for trait A')
parser$add_argument('--effect_blocks_b', type = 'character', nargs = '+', help = 'Effect blocks for trait B')
parser$add_argument('--odds_ratio_a', type = 'double', help = 'Odds ratio A')
parser$add_argument('--odds_ratio_b', type = 'double', help = 'Odds ratio B')
parser$add_argument('--P_a', type = 'double', help = 'Sample case proportion for trait A')
parser$add_argument('--K_a', type = 'double', help = 'Population case prevalence for trait A')
parser$add_argument('--P_b', type = 'double', help = 'Sample case proportion for trait B')
parser$add_argument('--K_b', type = 'double', help = 'Population case prevalence for trait B')
parser$add_argument('-o', '--output_path', type = 'character', help = 'Path to output file', required = T)
parser$add_argument('-nt', '--no_of_threads', type = 'integer', help = 'Number of threads to use', default = 1)

test_args <- c("--cv_file", "results/simgwas/combined_causal_variants.tsv",
               "--blocks_file", "resources/simgwas/available_blocks.tsv",
               "--effect_blocks_a", "1-m0:104",
               "--effect_blocks_b", "1-m0:9",
               "--odds_ratio_a", 1.2,
               "--odds_ratio_b", 1.2,
               "--P_a", 0.02,
               "--P_b", 0.02,
               "--K_a", 0.5,
               "--K_b", 0.5,
               "-o", "test.tsv",
               "-nt", 4)

args <- parser$parse_args()

setDTthreads(args$no_of_threads)

cv_dat <- fread(args$cv_file, sep = '\t', header = T)

blocks_dat <- fread(args$blocks_file, sep = '\t', header = T)

blocks_dat[, chr_block := paste(chr, block, sep = '_')]

get_chr_block_for_variants <- function(effect_blocks_string) {
  str_split(effect_blocks_string, "\\+") %>%
    unlist %>%
    str_match(., "(\\d+)-([simlh])(\\d+):(\\d+)") %>%
    apply(., 1, FUN = function(x) paste(as.integer(x[2]), as.integer(x[4]):as.integer(x[5]), sep = '_')) %>%
    unlist 
}

cv_dat[, geno_var := 2*EUR*(1-EUR)]

cv_dat[, chr_block := paste(chr, block, sep = '_')]

cv_dat[, in_effect_blocks_A := chr_block %in% get_chr_block_for_variants(args$effect_blocks_a) & chr_block %in% blocks_dat$chr_block]
cv_dat[, in_effect_blocks_B := chr_block %in% get_chr_block_for_variants(args$effect_blocks_b) & chr_block %in% blocks_dat$chr_block]

cv_dat[in_effect_blocks_A == T, beta.A := log(args$odds_ratio_a)]
cv_dat[in_effect_blocks_B == T, beta.B := log(args$odds_ratio_b)]

cv_dat[in_effect_blocks_A == T, beta_2.A := beta.A^2]
cv_dat[in_effect_blocks_B == T, beta_2.B := beta.B^2]

V_A.A <- with(cv_dat[in_effect_blocks_A == T], sum(beta_2.A*geno_var))
V_A.B <- with(cv_dat[in_effect_blocks_B == T], sum(beta_2.B*geno_var))

h2.theo.obs.A <- V_A.A/(args$P_a*(1-args$P_a))
h2.theo.obs.B <- V_A.B/(args$P_b*(1-args$P_b))

h2.theo.liab.A <- obs_liab_trans(h2.theo.obs.A, P = args$P_a, K = args$K_a)
h2.theo.liab.B <- obs_liab_trans(h2.theo.obs.B, P = args$P_b, K = args$K_b)

C_A.AB <- with(cv_dat[in_effect_blocks_A == T & in_effect_blocks_B == T], sum(beta.A * beta.B * geno_var))

r_A.AB <- C_A.AB/(sqrt(V_A.A)*sqrt(V_A.B))

res_dat <- data.table(odds_ratio.A = args$odds_ratio_a,
                      odds_ratio.B = args$odds_ratio_b,
                      no_blocks.A = cv_dat[in_effect_blocks_A == T, .N],
                      no_blocks.B = cv_dat[in_effect_blocks_B == T, .N],
                      no_shared_blocks = cv_dat[in_effect_blocks_A == T & in_effect_blocks_B == T, .N],
                      h2.theo.obs.A,
                      h2.theo.obs.B,
                      h2.theo.liab.A,
                      h2.theo.liab.B,
                      V_A.A,
                      V_A.B,
                      C_A.AB,
                      r_A.AB)

fwrite(res_dat, file = args$output_path, sep = '\t', na = 'NA')
