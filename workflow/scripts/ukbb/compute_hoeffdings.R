library(data.table)
library(Rcpp)
library(independence)

setDTthreads(snakemake@threads)

sourceCpp(code = '
#include <Rcpp.h>
#include <map>

using namespace Rcpp;

// [[Rcpp::export]]
NumericVector perturbDuplicates(NumericVector values) {
  std::map<double, int> freqMap;

  for(size_t i = 0; i < values.size(); ++i){
    freqMap[values[i]]++;

    if(freqMap[values[i]] > 1) {
      for(int j = 1; j < freqMap[values[i]]; ++j) {
        values[i] = values[i] + (freqMap[values[i]] * std::numeric_limits<double>::epsilon());
      }
    }
  }

  return values;
}
')

run_hoeffding <- function(dat, trait_A_code, trait_B_code) {
  dat[, trait_A_code := perturbDuplicates(get(trait_A_code))]
  dat[, trait_B_code := perturbDuplicates(get(trait_B_code))]

  dat <- unique(unique(dat, by = trait_A_code), by = trait_B_code)

  hoeffding.D.test(xs = dat[[trait_A_code]], ys = dat[[trait_B_code]])
}

sum_stats_dat <- fread(snakemake@input[['sum_stats_file']], sep = '\t', header = T, select = c(snakemake@wildcards[['trait_A']], snakemake@wildcards[['trait_B']]))

hoeffding_res <- run_hoeffding(sum_stats_dat, snakemake@wildcards[['trait_A']], snakemake@wildcards[['trait_B']])

print(hoeffding_res)

res_dat <- data.table(t(unlist(hoeffding_res[c('n', 'Dn', 'scaled', 'p.value')])))

res_dat[, `:=` (trait_A = snakemake@wildcards[['trait_A']], trait_B = snakemake@wildcards[['trait_B']])]

res_dat <- res_dat[, c('trait_A', 'trait_B', 'n', 'Dn', 'scaled', 'p.value')]

fwrite(res_dat, file = snakemake@output[[1]], sep = '\t', col.names = T, row.names = F)
