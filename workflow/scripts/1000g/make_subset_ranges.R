library(data.table)
setDTthreads(snakemake@threads)

bim <- fread(snakemake@input[['bim']], sep = '\t', header = F, col.names = c('chr', 'ID', 'Cm', 'bp', 'A1', 'A2'))

if(snakemake@wildcards[['snp_set']] == 'ukbb_sans_mhc') {
  bim <- bim[!(chr == 6 & bp %between% c(24e6, 45e6))]
}

if(snakemake@wildcards[['snp_set']] == 'ukbb_sans_mhc' | snakemake@wildcards[['snp_set']] == 'ukbb') {
  ukbb <- fread(snakemake@input[['ukbb']], sep = '\t', header = T, select = 'variant')
  
  bim[, variant_12 := paste(chr, bp, A1, A2, sep = ':')]
  bim[, variant_21 := paste(chr, bp, A2, A1, sep = ':')]

  bim <- bim[variant_12 %in% ukbb$variant | variant_21 %in% ukbb$variant]
}

fwrite(bim, sep = '\t', col.names = F, file = snakemake@output[[1]])


