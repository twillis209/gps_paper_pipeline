library(data.table)
library(stringr)
library(magrittr)

save.image('calculate_theoretical_rg_randomised_blocks.RData')

obs_liab_trans <- function(h2.obs, P, K) {
  z_2 <- dnorm(qnorm(1-K))^2

  h2.obs * ((K*(1-K))^2/(P*(1-P)))/z_2
}

odds_ratios <- list('null' = 1,
                    'tiny' = 1.02,
                    'small' = 1.05,
                    'infinitesimal' = 1.1,
                    'medium' = 1.2,
                    'large' = 1.4,
                    'vlarge' = 2)

setDTthreads(snakemake@threads)
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3059431/
# According to the paper in which the ascertainment-corrected liability scale transformation is set out, the convention is that P is the sample prevalence and K the population prevalence
P_a <- snakemake@params[['sample_prevalence_A']]
P_b <- snakemake@params[['sample_prevalence_B']]
K_a <- snakemake@params[['population_prevalence_A']]
K_b <- snakemake@params[['population_prevalence_B']]

cv_dat <- fread(snakemake@input[['combined_causal_variants_file']], sep = '\t', header = T)

a_blocks_dat <- fread(snakemake@input[['a_block_file']], sep = '\t', header = T)
b_blocks_dat <- fread(snakemake@input[['b_block_file']], sep = '\t', header = T)

cv_dat[, geno_var := 2*EUR*(1-EUR)]

cv_dat[, c('in_a_blocks', 'in_b_blocks') := F]
cv_dat[, c('odds_ratio_a', 'odds_ratio_b') := 1]

for(i in 1:nrow(a_blocks_dat)) {
  if(a_blocks_dat[i, no_cvs] > 0) {
    if(cv_dat[a_blocks_dat[i, chr] == chr & a_blocks_dat[i, block] == block, .N] == 0) {
      print(sprintf("Missing chr%d:block%d in a file", a_blocks_dat[i, chr], a_blocks_dat[i, block]))
      } else {
        for(j in 1:a_blocks_dat[i, no_cvs]) {
        # Need to iterate over j here otherwise we'll get two cvs where we sometimes only want one
          cv_dat[a_blocks_dat[i, chr] == chr & a_blocks_dat[i, block] == block & a_blocks_dat[i, effect] != 'null'][j] <- cv_dat[a_blocks_dat[i, chr] == chr & a_blocks_dat[i, block] == block & a_blocks_dat[i, effect] != 'null'][j][, `:=` (odds_ratio_a = unlist(odds_ratios[[a_blocks_dat[i, effect]]]), in_a_blocks = T)]
        }
    }
  }
}

for(i in 1:nrow(b_blocks_dat)) {
  if(b_blocks_dat[i, no_cvs] > 0) {
    if(cv_dat[b_blocks_dat[i, chr] == chr & b_blocks_dat[i, block] == block, .N] == 0) {
      print(sprintf("Missing chr%d:block%d in a file", b_blocks_dat[i, chr], b_blocks_dat[i, block]))
    } else {
      for(j in 1:b_blocks_dat[i, no_cvs]) {
      # Need to iterate over j here otherwise we'll get two cvs where we sometimes only want one
        cv_dat[b_blocks_dat[i, chr] == chr & b_blocks_dat[i, block] == block & b_blocks_dat[i, effect] != 'null'][j] <- cv_dat[b_blocks_dat[i, chr] == chr & b_blocks_dat[i, block] == block & b_blocks_dat[i, effect] != 'null'][j][, `:=` (odds_ratio_b = unlist(odds_ratios[[b_blocks_dat[i, effect]]]), in_b_blocks = T)]
      }
    }
  }
}

if(cv_dat[in_a_blocks == T, .N] != a_blocks_dat[, sum(no_cvs)]) stop(sprintf('Missing %d causal variants from A set', a_blocks_dat[, sum(no_cvs)] - cv_dat[in_a_blocks == T, .N] ))

if(cv_dat[in_b_blocks == T, .N] != b_blocks_dat[, sum(no_cvs)]) stop(sprintf('Missing %d causal variants from B set', b_blocks_dat[, sum(no_cvs)] - cv_dat[in_b_blocks == T, .N]))

if(cv_dat[in_a_blocks == T & in_b_blocks == T, .N] != merge(a_blocks_dat, b_blocks_dat, by = c('chr', 'block'))[effect.x != 'null' & effect.y != 'null', sum(no_cvs.x)]) stop('Missing shared causal variants')

cv_dat[in_a_blocks == T, beta.A := log(odds_ratio_a)]
cv_dat[in_b_blocks == T, beta.B := log(odds_ratio_b)]

cv_dat[in_a_blocks == T, beta_2.A := beta.A^2]
cv_dat[in_b_blocks == T, beta_2.B := beta.B^2]

V_A.A <- with(cv_dat[in_a_blocks == T], sum(beta_2.A*geno_var))
V_A.B <- with(cv_dat[in_b_blocks == T], sum(beta_2.B*geno_var))

h2.theo.obs.A <- V_A.A/(P_a*(1-P_a))
h2.theo.obs.B <- V_A.B/(P_b*(1-P_b))

h2.theo.liab.A <- obs_liab_trans(h2.theo.obs.A, P = P_a, K = K_a)
h2.theo.liab.B <- obs_liab_trans(h2.theo.obs.B, P = P_b, K = K_b)

C_A.AB <- with(cv_dat[in_a_blocks == T & in_b_blocks == T], sum(beta.A * beta.B * geno_var))

r_A.AB <- C_A.AB/(sqrt(V_A.A)*sqrt(V_A.B))

res_dat <- data.table(odds_ratio.A = paste(unique(cv_dat[in_a_blocks == T, odds_ratio_a]), collapse = ','),
                      odds_ratio.B = paste(unique(cv_dat[in_b_blocks == T, odds_ratio_b]), collapse = ','),
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

fwrite(res_dat, file = snakemake@output[['theo_rg_file']], sep = '\t', na = 'NA')
