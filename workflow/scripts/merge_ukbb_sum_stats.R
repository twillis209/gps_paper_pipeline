library(data.table)
setDTthreads(snakemake@threads)

ukbb_trait_codes <- snakemake@params[['ukbb_trait_codes']]
sum_stats_files <- snakemake@input[['ukbb_files']]

left_dat <- fread(sum_stats_files[1], sep = '\t', header = T, select = c('variant', 'pval', 'tstat', 'n_complete_samples'))

setnames(left_dat, c('variant', 'pval', 'tstat', 'n_complete_samples'), c('variant', paste0('pval.', ukbb_trait_codes[1]), paste0('tstat.', ukbb_trait_codes[1]), paste0('n_complete_samples.', ukbb_trait_codes[1])))

for(i in 2:length(sum_stats_files)) {
  right_dat <- fread(sum_stats_files[i], sep = '\t', header = T, select = c('variant', 'pval', 'tstat', 'n_complete_samples'))

  setnames(right_dat, c('variant', 'pval', 'tstat', 'n_complete_samples'), c('variant', paste0('pval.', ukbb_trait_codes[i]), paste0('tstat.', ukbb_trait_codes[i]), paste0('n_complete_samples.', ukbb_trait_codes[i])))

  left_dat <- merge(left_dat, right_dat, by = 'variant')
}

if(snakemake@params[['sans_mhc']]) {
  left_dat[, c('chr', 'bp') := tstrsplit(variant, split = ':', keep = 1:2)]

  left_dat <- left_dat[!(chr == 6 & bp %between% c(24e6, 45e6))]

  left_dat[, c('chr', 'bp') := NULL]
}

fwrite(left_dat, file = snakemake@output[[1]], sep = '\t')
