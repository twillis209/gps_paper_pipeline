library(data.table)

setDTthreads(snakemake@threads)

intermediates_dat <- fread(snakemake@input[['intermediates_file']], sep = '\t', header = T)
sum_stats_dat <- fread(snakemake@input[['sum_stats_file']], sep = '\t', header = T, select = 'variant')

sum_stats_dat[, c('CHR19', 'BP19', 'ALT', 'REF') := tstrsplit(variant, split = ':', keep = 1:4)]

if(sum_stats_dat[, .N] != intermediates_dat[, .N]) stop("No. of rows in sum stats and intermediate output files differs")

out_dat <- cbind(sum_stats_dat, intermediates_dat)

fwrite(out_dat, file = snakemake@output[[1]], sep = '\t', col.names = T)
