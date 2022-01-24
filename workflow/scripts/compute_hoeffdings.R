library(data.table)
library(Rcpp)
library(argparse)
library(independence)

parser <- ArgumentParser(description = 'Run Hoeffding\'s test on a pair of traits')
parser$add_argument('-i', '--sum_stats_file', type = 'character', help = 'Path to summary statistics file')
parser$add_argument('-a', '--trait_A', type = 'character', help = 'Trait A code')
parser$add_argument('-b', '--trait_B', type = 'character', help = 'Trait B code')
parser$add_argument('-o', '--output_path', type = 'character', help = 'Path to output file', required = T)
parser$add_argument('-nt', '--no_of_threads', type = 'integer', help = 'Number of threads to use', default = 1)

args <- parser$parse_args()

setDTthreads(args$no_of_threads)

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

sum_stats_dat <- fread(args$sum_stats_file, sep = '\t', header = T, select = c(args$trait_A, args$trait_B))

hoeffding_res <- run_hoeffding(sum_stats_dat, args$trait_A, args$trait_B)

print(hoeffding_res)

res_dat <- data.table(t(unlist(hoeffding_res[c('n', 'Dn', 'scaled', 'p.value')])))

res_dat[, `:=` (trait_A = args$trait_A, trait_B = args$trait_B)]

res_dat <- res_dat[, c('trait_A', 'trait_B', 'n', 'Dn', 'scaled', 'p.value')]

fwrite(res_dat, file = args$output_path, sep = '\t', col.names = T, row.names = F)
