library(data.table)

setDTthreads(snakemake@threads)

chr_colname <- snakemake@params[['chr_colname']]
bp_colname <- snakemake@params[['bp_colname']]
a1_colname <- snakemake@params[['a1_colname']]
a2_colname <- snakemake@params[['a2_colname']]
z_colname <- snakemake@params[['z_colname']]
n <- snakemake@params[['sample_size']]

dat <- fread(snakemake@input[[1]], sep = '\t', header = T, select = c(chr_colname, bp_colname, a1_colname, a2_colname, z_colname))

setnames(dat, c(a1_colname, a2_colname, z_colname), c('A2', 'A1', 'Z'))

dat[, Predictor := paste(get(chr_colname), get(bp_colname), A2, A1, sep = ':')]

dat[, n := n]

dat[, `:=` (len.A1 = nchar(A1), len.A2 = nchar(A2))]

dat <- dat[len.A1 == 1 & len.A2 == 1]

dat <- dat[!duplicated(dat, by = 'Predictor')]

dat <- dat[, .(Predictor, A1, A2, n, Z)]

fwrite(dat, file = snakemake@output[[1]], sep = '\t')
