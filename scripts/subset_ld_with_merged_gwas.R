library(data.table)
library(argparse)

PID_ROOT <- Sys.getenv('pidRoot')

parser <- ArgumentParser(description = 'Subset 1kGP SNPs with LD values based on SNPs present in merged GWAS summary statistics file')
parser$add_argument('-i', '--gwas_file', type = 'character', help = 'Path to merged GWAS summary statistics file')
parser$add_argument('-ld', '--ld_root', type = 'character', help = 'Path to directory containing LD statistics files')
parser$add_argument('-o', '--output_path', type = 'character', help = 'Path to output directory to which subsetted LD files will be written', required = T)
parser$add_argument('-chr', type = 'character', help = 'Label of chromosome column in GWAS file A', default = 'CHR38')
parser$add_argument('-bp', type = 'character', help = 'Label of BP column in GWAS file', default = 'BP38')
parser$add_argument('-ref', type = 'character', help = 'Label of reference allele column in GWAS file', default = 'REF')
parser$add_argument('-alt', type = 'character', help = 'Label of alternative allele column in GWAS file', default = 'ALT')
parser$add_argument('-nt', '--no_of_threads', type = 'integer', help = 'Number of threads to use', default = 1)

args_vec <- c('-i', 'gwas/pid_aster.tsv.gz', '-ld', '1000g/euro/qc/ld', '-o', 'gwas/pid_aster/ld', '-nt', '8')

args <- parser$parse_args()

setDTthreads(args$no_of_threads)

gwas_dat <- fread(args$gwas_file, sep = '\t', header = T, select = c(args$chr, args$bp, args$ref, args$alt))

gwas_dat[ , (args$chr) := as.character(get(args$chr))]

# TODO doesn't work for chrX as LD file has no IDs in it
# TODO for indels, max 10 characters of alt allele included in the ID made for the LD file

gwas_dat[, 'alt_ref' := paste(get(args$chr), get(args$bp), get(args$alt), get(args$ref), sep = ':')]
gwas_dat[, 'ref_alt' := paste(get(args$chr), get(args$bp), get(args$ref), get(args$alt), sep = ':')]

# Only handling 1-22 at the moment
gwas_dat <- gwas_dat[CHR38 %in% as.character(1:22)]

for(i in 1:22) {
  ld_dat <- fread(file.path(args$ld_root, sprintf('chr%d.ld', i)), sep = ' ', header = T)

  chr_dat <- gwas_dat[CHR38 == as.character(i)]

  ld_dat <- ld_dat[SNP_A %in% chr_dat$alt_ref | SNP_A %in% chr_dat$ref_alt]

  fwrite(ld_dat, file = file.path(args$output_path, sprintf('chr%d.ld', i)), sep = ' ')
}
