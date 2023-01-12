library(data.table)
setDTthreads(snakemake@threads)

dat <- fread(snakemake@input[[1]])

for(i in 1:22) {
  fwrite(dat[chr == i, env = list(chr = snakemake@params[['chr_col']])], file = snakemake@output[[i]], sep = '\t')
}
