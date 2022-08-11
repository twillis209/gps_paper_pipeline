library(data.table)
library(simGWAS)
library(argparse)
library(parallel)

parser <- ArgumentParser(description = 'Computes blockwise LD matrix')
parser$add_argument('--hap_file', type = 'character', help = 'Path to haplotype file')
parser$add_argument('--leg_file', type = 'character', help = 'Path to legend file')
parser$add_argument('--output_file', type = 'character', help = 'Path to output file', required = T)
parser$add_argument('-nt', '--no_of_threads', type = 'integer', help = 'Number of threads to use', default = 1)

args <- parser$parse_args()

setDTthreads(args$no_of_threads)

leg_dat <- fread(file = args$leg_file, sep = ' ', header = T)
hap_dat <- fread(file = args$hap_file, sep = ' ', header = F)

for(j in 1:(ncol(hap_dat)-2)) {
  set(hap_dat, j = j, value = as.numeric(hap_dat[[j]]))
}

hap_mat <- as.matrix(hap_dat[,1:(ncol(hap_dat)-2)])

freq_dat <- data.table(t(hap_mat)+1)

colnames(freq_dat) <- hap_dat[[(ncol(hap_dat)-1)]]

freq_dat[, Probability := 1/.N]

rm(hap_dat, hap_mat)

ld_mat <- corpcor::make.positive.definite(simGWAS:::wcor2(as.matrix(freq_dat[,setdiff(colnames(freq_dat),"Probability"), with = F][, leg_dat$rs, with = F]), freq_dat$Probability))

save(ld_mat, file = args$output_file)
