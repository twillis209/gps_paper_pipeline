library(data.table)
library(stringr)
library(magrittr)

save.image('calculate_theoretical_rg_randomised_blocks.RData')

obs_liab_trans <- function(h2.obs, P, K) {
  z_2 <- dnorm(qnorm(1-K))^2

  h2.obs * ((K*(1-K))^2/(P*(1-P)))/z_2
}

odds_ratios <- list('null' = 1,
                    'tiny' = 1.01,
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

for(k in 1:2) {
  if(k == 1) {
    x <- a_blocks_dat
  } else {
    x <- b_blocks_dat
  }

  in_blocks_label <- ifelse(k == 1, 'in_a_blocks', 'in_b_blocks')
  odds_ratio_label <- ifelse(k == 1, 'odds_ratio_a', 'odds_ratio_b')

  for(i in 1:nrow(x)) {
    if(x[i, no_cvs] > 0) {
        if(cv_dat[x[i, chr] == chr & x[i, block] == block, .N] == 0) {
            print(sprintf("Missing chr%d:block%d in a file", x[i, chr], x[i, block]))
        } else {
          for(j in 1:x[i, no_cvs]) {
          # Need to iterate over j here otherwise we'll get two cvs where we sometimes only want one
            cv_dat[x[i, chr] == chr & x[i, block] == block & x[i, effect] != 'null'][j, `:=` (in_blocks = T, odds_ratio = unlist(odds_ratios[[x[i, effect]]])), env = list(in_blocks = in_blocks_label, odds_ratio = odds_ratio_label)]
          }
      }
    }
  }
}

if(cv_dat[in_a_blocks == T, .N] != a_blocks_dat[, sum(no_cvs)]) stop(sprintf('Missing %d causal variants from A set', a_blocks_dat[, sum(no_cvs)] - cv_dat[in_a_blocks == T, .N] ))

if(cv_dat[in_b_blocks == T, .N] != b_blocks_dat[, sum(no_cvs)]) stop(sprintf('Missing %d causal variants from B set', b_blocks_dat[, sum(no_cvs)] - cv_dat[in_b_blocks == T, .N]))

if(cv_dat[in_a_blocks == T & in_b_blocks == T, .N] != merge(a_blocks_dat, b_blocks_dat, by = c('chr', 'block'))[effect.x != 'null' & effect.y != 'null', sum(no_cvs)]) stop('Missing shared causal variants')

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
