library(magrittr)
library(data.table)

setDTthreads(snakemake@threads)

lapply(snakemake@input[['ldsc_files']], fread, sep = '\t', select = c('CHR', 'SNP', 'BP')) %>% rbindlist -> merged_dat

allele_dat <- fread(snakemake@input[['snplist_file']], sep = '\t')

merged_dat <- merge(merged_dat, allele_dat, by = 'SNP')

fwrite(merged_dat, file = snakemake@output[[1]], sep = '\t')
