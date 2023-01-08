daf <- read.table(snakemake@input[[1]], sep = '\t', header = T)

out_daf <- data.frame(gps = daf$GPS, pval = pexp(daf$GPS^-2))

write.table(out_daf, file = snakemake@output[[1]], sep = '\t')

