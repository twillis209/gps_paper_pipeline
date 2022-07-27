library(data.table)

setDTthreads(snakemake@resources[['threads']])

sum_stats_dat <- fread(snakemake@input[[1]], sep = '\t', header = T, tmpdir = snakemake@params[['tmpdir']])

sum_stats_dat[, c('chr', 'bp') := tstrsplit(variant, split = ':', keep = 1:2)]

sum_stats_dat <- sum_stats_dat[!(chr == 6 & bp %between% c(24e6, 45e6))]

sum_stats_dat[, c('chr', 'bp') := NULL]

fwrite(sum_stats_dat, file = snakemake@output[[1]], sep = '\t')
