library(data.table)
setDTthreads(snakemake@threads)

dat <- fread(snakemake@input[[1]], sep = '\t', header = T)

dat[, c('chr', 'bp') := tstrsplit(variant, split = ':', keep = 1:2)]

dat <- dat[!(chr == 6 & bp %between% c(24e6, 45e6))]

dat[, c('chr', 'bp') := NULL]

fwrite(dat, file = snakemake@output[[1]], sep = '\t')
