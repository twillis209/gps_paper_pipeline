library(data.table)
setDTthreads(snakemake@threads)

dat <- fread(snakemake@input[[1]], sep = '\t', header = T)

dat[f.6148.0.0 == 2 | f.6148.0.1 == 2 | f.6148.0.2 == 2 | f.6148.0.3 == 2 | f.6148.0.4 == 2, glaucoma := "6148_2"]
dat[f.6148.0.0 == 5 | f.6148.0.1 == 5 | f.6148.0.2 == 5 | f.6148.0.3 == 5 | f.6148.0.4 == 5, md := "6148_5"]

dat[f.22126.0.0 == 0, f.22126.0.0 := NA]
dat[f.22126.0.0 == 1, f.22126.0.0 := 22126]

fwrite(dat, file = snakemake@output[[1]], sep = '\t')

