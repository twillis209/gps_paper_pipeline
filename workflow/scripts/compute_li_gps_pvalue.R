daf <- read.table(snakemake@input[[1]], sep = '\t')

write.table(data.frame(gps = daf$GPS, pval = pexp(daf$GPS^-2)), file = snakemake@output[[1]], sep = '\t')
