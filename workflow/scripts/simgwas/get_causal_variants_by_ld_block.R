library(data.table)

setDTthreads(snakemake@threads)

leg_dat <- fread(file = snakemake@input[['block_legend_file']], sep = ' ', header = T)
hap_dat <- fread(file = snakemake@input[['block_haplotype_file']], sep = ' ', header = F)
bim_dat <- fread(file = snakemake@input[['bim_file']], sep = '\t', header = F, col.names = c('chr', 'bim.id', 'Cm', 'bp', 'A1', 'A2'))
load(file = snakemake@input[['ld_mat_file']])

for(j in 1:(ncol(hap_dat)-2)) {
  set(hap_dat, j = j, value = as.numeric(hap_dat[[j]]))
}

hap_mat <- as.matrix(hap_dat[,1:(ncol(hap_dat)-2)])

hap_meta_dat <- hap_dat[, (ncol(hap_dat)-1):ncol(hap_dat)]

names(hap_meta_dat) <- c('rs', 'block')

rm(hap_dat)

freq_dat <- data.table(t(hap_mat)+1)

rm(hap_mat)

colnames(freq_dat) <- hap_meta_dat$rs

# Required by the make_GenoProbList function below (at least)
freq_dat[, Probability := 1/.N]

cv_ind <- snakemake@params[['causal_variant_indices']]

result_dat <- data.table(leg_dat[cv_ind, .(id, position, block, a0, a1, TYPE, EUR)])

result_dat <- merge(result_dat, bim_dat[, .(bim.id, bp, A1, A2)], by.x = 'position', by.y = 'bp', all.x = T)

result_dat[, c("A1", "A2", "bim.id") := NULL]

result_dat[, chr := snakemake@params[['chr_no']]]

result_dat[, rsID := tstrsplit(id, split = ':', keep = 1)]

fwrite(result_dat, file = snakemake@output[[1]], sep = '\t')
