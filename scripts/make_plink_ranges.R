library(data.table)
library(argparse)

parser <- ArgumentParser(description = 'Produce PLINK SNP ID files to extract designated SNPs from PLINK files')
parser$add_argument('-i', '--gwas_file', type = 'character', help = 'Path to merged GWAS summary statistics file')
parser$add_argument('-b', '--bim_dir', type = 'character', help = 'Path to directory containing bim files')
parser$add_argument('-r', '--bim_regex', type = 'character', help = 'Regex for bim files. Must contain \'%d\' for autosome numbers.', default = 'chr%d.bim')
parser$add_argument('-chr', type = 'character', help = 'Label of chromosome column in GWAS file', default = 'CHR38')
parser$add_argument('-bp', type = 'character', help = 'Label of BP column in GWAS file', default = 'BP38')
parser$add_argument('-ref', type = 'character', help = 'Label of reference allele column in GWAS file', default = 'REF')
parser$add_argument('-alt', type = 'character', help = 'Label of alternative allele column in GWAS file', default = 'ALT')
parser$add_argument('-prin', type = 'character', help = 'Label of principal p-value column', default = 'P.A')
parser$add_argument('-aux', type = 'character', help = 'Label of auxiliary p-value column', default = 'P.B')
parser$add_argument('-o', '--output_dir', type = 'character', help = 'Path to output directory', required = T)
parser$add_argument('-nt', '--no_of_threads', type = 'integer', help = 'Number of threads to use', default = 1)

# args_vec <- c('-i', 'gwas/pid_aster.tsv.gz', '-b', '1000g/euro/qc', '-r', 'chr%d_qc.bim', '-aux', 'P.B', '-o' , 'test.bim', '-nt', 8)

args<-parser$parse_args()

setDTthreads(threads=args$no_of_threads)

gwas_dat <- fread(args$gwas_file, sep = '\t', header = T, select = c(args$chr, args$bp, args$ref, args$alt, args$prin, args$aux))

gwas_dat[ , args$chr := as.character(get(args$chr))]

#gwas_dat[, ref_short := ifelse(nchar(REF) > 10, substr(REF, 1, 10), REF)]
#gwas_dat[, alt_short := ifelse(nchar(ALT) > 10, substr(ALT, 1, 10), ALT)]
#
#gwas_dat[, 'ref_alt' := paste(get(args$chr), get(args$bp), ref_short, alt_short, sep = ':')]
#gwas_dat[, 'alt_ref' := paste(get(args$chr), get(args$bp), alt_short, ref_short, sep = ':')]
#
#gwas_dat[, c('ref_short', 'alt_short') := NULL]

for(i in 1:22) {
  bim_dat <- fread(file.path(args$bim_dir, sprintf(args$bim_regex, i)), sep = '\t', header = F, col.names = c('CHR38', 'ID', 'Cm', 'BP38', 'A1', 'A2'))

  bim_dat[, CHR38 := as.character(CHR38)]

  bim_join <- merge(bim_dat, gwas_dat, by.x = c(args$chr, args$bp), by.y = c('CHR38', 'BP38'), sort = F)

#  bim_join[, freq := table(BP38)[as.character(BP38)]]

#  bim_join <- bim_join[freq == 1]

  # Make sure alleles match
  bim_join <- bim_join[(REF == A1 & ALT == A2) | (REF == A2 & ALT == A1)]

  #  bim_join <- bim_join[, .(CHR38, BP38, BP38, ID)]

  bim_join <- bim_join[, .(ID)]

  fwrite(bim_join, file = file.path(args$output_dir, sprintf("chr%d.txt", i)), row.names = F, sep = ' ', col.names = F, quote = F)
}
