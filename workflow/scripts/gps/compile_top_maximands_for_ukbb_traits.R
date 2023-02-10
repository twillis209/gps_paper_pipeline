library(data.table)
library(magrittr)
library(stringr)

setDTthreads(snakemake@threads)

for(i in seq_along(snakemake@input[['annot_files']])) {
  base <- basename(snakemake@input[['annot_files']][i])

  str_replace(base, pattern = '_intermediates_annot.tsv', replacement = '')[[1]] %>%
    str_split(., pattern = '-') %>%
    unlist -> trait_pair

  dat <- fread(snakemake@input[['annot_files']][i], sep = '\t', header = T)

  dat <- dat[order(maximand, decreasing = T)][1]

  dat[, `:=` (trait_A = trait_pair[[1]], trait_B = trait_pair[[2]])]

  pval_dat <- fread(snakemake@input[['pvalue_files']][i], sep = '\t', header = T)

  if(i == 1) {
    out_dat <- cbind(dat, pval_dat[, .(gps, pval)])
  } else {
    out_dat <- rbind(out_dat, cbind(dat, pval_dat[, .(gps, pval)]))
  }
}

fwrite(out_dat, file = snakemake@output[[1]], sep = '\t', col.names = T)
