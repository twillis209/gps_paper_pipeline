library(data.table)

setDTthreads(snakemake@threads)

sum_stats_dat <- fread(snakemake@input[[1]], sep = '\t', header = T)

no_snps <- min(nrow(sum_stats_dat), as.integer(snakemake@wildcards[['no_snps']]))

sum_stats_dat <- sum_stats_dat[sample(1:no_snps)]

fwrite(sum_stats_dat, file = snakemake@output[[1]], sep = '\t')
