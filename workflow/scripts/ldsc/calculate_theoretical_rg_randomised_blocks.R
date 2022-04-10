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
parser$add_argument('--a_blocks_file', type = 'character')
parser$add_argument('--b_blocks_file', type = 'character')
parser$add_argument('--odds_ratio_a', type = 'double', help = 'Odds ratio A')
parser$add_argument('--odds_ratio_b', type = 'double', help = 'Odds ratio B')
parser$add_argument('--P_a', type = 'double', help = 'Sample case proportion for trait A')
parser$add_argument('--K_a', type = 'double', help = 'Population case prevalence for trait A')
parser$add_argument('--P_b', type = 'double', help = 'Sample case proportion for trait B')
parser$add_argument('--K_b', type = 'double', help = 'Population case prevalence for trait B')
parser$add_argument('-o', '--output_path', type = 'character', help = 'Path to output file', required = T)
parser$add_argument('-nt', '--no_of_threads', type = 'integer', help = 'Number of threads to use', default = 1)

#test_args <- c("--cv_file", "results/simgwas/combined_causal_variants.tsv", "--a_blocks_file", "results/ldsc/rg/whole_genome/randomised/theoretical_rg/block_files/m50_m50_m10_seed_111_a.tsv", "--b_blocks_file", "results/ldsc/rg/whole_genome/randomised/theoretical_rg/block_files/m50_m50_m10_seed_111_b.tsv", "--odds_ratio_a", "1.2", "--odds_ratio_b", "1.2", "--P_a", "0.5", "--P_b", "0.5", "--K_a", "0.02", "--K_b", "0.02", "-o", "results/ldsc/rg/whole_genome/randomised/theoretical_rg/m50_m50_m10_seed_111_theo_rg.tsv", "-nt", 1)

test_args <- c("--cv_file", "results/simgwas/combined_causal_variants.tsv", "--a_blocks_file", "results/ldsc/rg/whole_genome/randomised/theoretical_rg/1000_1000_1000_1000/block_files/m50_m50_m0_seed_196_q_qr.tsv", "--b_blocks_file", "results/ldsc/rg/whole_genome/randomised/theoretical_rg/1000_1000_1000_1000/block_files/m50_m50_m0_seed_196_r_qr.tsv", "--odds_ratio_a", "1.2", "--odds_ratio_b", "1.2", "--P_a", "0.5", "--P_b", "0.5", "--K_a", "0.02", "--K_b", "0.02", "-o", "results/ldsc/rg/whole_genome/randomised/theoretical_rg/1000_1000_1000_1000/m50_m50_m0_seed_196_qr_theo_rg.tsv", "-nt", 1)
#
args <- parser$parse_args(test_args)

args <- parser$parse_args()

setDTthreads(args$no_of_threads)

cv_dat <- fread(args$cv_file, sep = '\t', header = T)

a_blocks_dat <- fread(args$a_blocks_file, sep = '\t', header = T)
b_blocks_dat <- fread(args$b_blocks_file, sep = '\t', header = T)

cv_dat[, geno_var := 2*EUR*(1-EUR)]

cv_dat[, c('in_a_blocks', 'in_b_blocks') := F]

for(i in seq_along(1:nrow(a_blocks_dat))) {
  if(cv_dat[a_blocks_dat[i, chr] == chr & a_blocks_dat[i, block] == block, .N] == 0) {
    print(sprintf("Missing chr%d:block%d in a file", a_blocks_dat[i, chr], a_blocks_dat[i, block]))
  } else {
    cv_dat[a_blocks_dat[i, chr] == chr & a_blocks_dat[i, block] == block, in_a_blocks := T]
  }
}

for(i in seq_along(1:nrow(b_blocks_dat))) {
  if(cv_dat[b_blocks_dat[i, chr] == chr & b_blocks_dat[i, block] == block, .N] == 0) {
    print(sprintf("Missing chr%d:block%d in b file", b_blocks_dat[i, chr], b_blocks_dat[i, block]))
  } else {
    cv_dat[b_blocks_dat[i, chr] == chr & b_blocks_dat[i, block] == block, in_b_blocks := T]
  }
}

if(cv_dat[in_a_blocks == T, .N] != a_blocks_dat[, .N]) stop(sprintf('Missing %d causal variants from A set', a_blocks_dat[, .N] - cv_dat[in_a_blocks == T, .N] ))

if(cv_dat[in_b_blocks == T, .N] != b_blocks_dat[, .N]) stop(sprintf('Missing %d causal variants from B set', b_blocks_dat[, .N] - cv_dat[in_b_blocks == T, .N]))

if(cv_dat[in_a_blocks == T & in_b_blocks == T, .N] != merge(a_blocks_dat, b_blocks_dat, by = c('chr', 'block'))[, .N]) stop('Missing shared causal variants')


cv_dat[in_a_blocks == T, beta.A := log(args$odds_ratio_a)]
cv_dat[in_b_blocks == T, beta.B := log(args$odds_ratio_b)]

cv_dat[in_a_blocks == T, beta_2.A := beta.A^2]
cv_dat[in_b_blocks == T, beta_2.B := beta.B^2]

V_A.A <- with(cv_dat[in_a_blocks == T], sum(beta_2.A*geno_var))
V_A.B <- with(cv_dat[in_b_blocks == T], sum(beta_2.B*geno_var))

h2.theo.obs.A <- V_A.A/(args$P_a*(1-args$P_a))
h2.theo.obs.B <- V_A.B/(args$P_b*(1-args$P_b))

h2.theo.liab.A <- obs_liab_trans(h2.theo.obs.A, P = args$P_a, K = args$K_a)
h2.theo.liab.B <- obs_liab_trans(h2.theo.obs.B, P = args$P_b, K = args$K_b)

C_A.AB <- with(cv_dat[in_a_blocks == T & in_b_blocks == T], sum(beta.A * beta.B * geno_var))

r_A.AB <- C_A.AB/(sqrt(V_A.A)*sqrt(V_A.B))

res_dat <- data.table(odds_ratio.A = args$odds_ratio_a,
                      odds_ratio.B = args$odds_ratio_b,
                      no_blocks.A = cv_dat[in_a_blocks == T, .N],
                      no_blocks.B = cv_dat[in_b_blocks == T, .N],
                      no_shared_blocks = cv_dat[in_a_blocks == T & in_b_blocks == T, .N],
                      h2.theo.obs.A,
                      h2.theo.obs.B,
                      h2.theo.liab.A,
                      h2.theo.liab.B,
                      V_A.A,
                      V_A.B,
                      C_A.AB,
                      r_A.AB)

fwrite(res_dat, file = args$output_path, sep = '\t', na = 'NA')
