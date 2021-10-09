library(data.table)
library(argparse)

parser <- ArgumentParser(description = 'Extract LD-pruned subset of SNPs from joined GWAS file.')
parser$add_argument('-i', '--gwas_file', type = 'character', help = 'Path to merged GWAS summary statistics file')
parser$add_argument('-p', '--prune_file', type = 'character', help = 'Path to file containing IDs of SNPs in LD-pruned subset')
parser$add_argument('-chr', type = 'character', help = 'Label of chromosome column in GWAS file', default = 'CHR38')
parser$add_argument('-bp', type = 'character', help = 'Label of BP column in GWAS file', default = 'BP38')
parser$add_argument('-ref', type = 'character', help = 'Label of reference allele column in GWAS file', default = 'REF')
parser$add_argument('-alt', type = 'character', help = 'Label of alternative allele column in GWAS file', default = 'ALT')
parser$add_argument('-prin', type = 'character', help = 'Label of principal p-value column', default = 'P.A')
parser$add_argument('-aux', type = 'character', help = 'Label of auxiliary p-value column', default = 'P.B')
parser$add_argument('-o', '--output_file', type = 'character', help = 'Path to output file', required = T)
parser$add_argument('-nt', '--no_of_threads', type = 'integer', help = 'Number of threads to use', default = 1)

#args_vec <- c('-i', 'gwas/pid_aster.tsv.gz', '-p', 'gwas/pid_aster/prune/all.prune.in', '-o', 'gwas/pid_aster/pruned/pid_aster.tsv.gz', '-nt', 8)
#args<-parser$parse_args(args_vec)

args<-parser$parse_args()

setDTthreads(threads=args$no_of_threads)

gwas_dat <- fread(args$gwas_file, sep = '\t', header = T, select = c(args$chr, args$bp, args$ref, args$alt, args$prin, args$aux))

prune_dat <- fread(args$prune_file, sep = ' ', header = F)

#prune_dat[, c('CHR38', 'BP38', 'ALT', 'REF') := tstrsplit(V1, ':')]
#prune_dat[, V1 := NULL]

gwas_dat[, ref_short := ifelse(nchar(REF) > 10, substr(REF, 1, 10), REF)]
gwas_dat[, alt_short := ifelse(nchar(ALT) > 10, substr(ALT, 1, 10), ALT)]

gwas_dat[, 'ref_alt' := paste(get(args$chr), get(args$bp), ref_short, alt_short, sep = ':')]
gwas_dat[, 'alt_ref' := paste(get(args$chr), get(args$bp), alt_short, ref_short, sep = ':')]

gwas_dat[, c('ref_short', 'alt_short') := NULL]

gwas_dat <- gwas_dat[alt_ref %in% prune_dat$V1 | ref_alt %in% prune_dat$V1]

gwas_dat[, c('ref_alt', 'alt_ref') := NULL]

fwrite(gwas_dat, file = args$output_file, sep = '\t', col.names = T, row.names = F, quote = F)
