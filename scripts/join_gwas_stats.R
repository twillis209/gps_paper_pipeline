library(data.table)
library(argparse)

PID_ROOT <- Sys.getenv('pidRoot')

parser <- ArgumentParser(description = 'Join summary statistics in pair of GWAS files')
parser$add_argument('-a', '--input_file_a', type = 'character', help = 'Path to GWAS summary statistics file for trait A')
parser$add_argument('-b', '--input_file_b', type = 'character', help = 'Path to GWAS summary statistics file for trait B', default = file.path(PID_ROOT, 'sourceData/pid.tsv.gz'))
parser$add_argument('-chr_a', type = 'character', help = 'Label of chromosome column in input file A', default = 'CHR38')
parser$add_argument('-chr_b', type = 'character', help = 'Label of chromosome column in input file B', default = 'CHR38')
parser$add_argument('-bp_a', type = 'character', help = 'Label of BP column in input file A', default = 'BP38')
parser$add_argument('-bp_b', type = 'character', help = 'Label of BP column in input file B', default = 'BP38')
parser$add_argument('-ref_a', type = 'character', help = 'Label of reference allele column in input file A', default = 'REF')
parser$add_argument('-ref_b', type = 'character', help = 'Label of reference allele column in input file B', default = 'REF')
parser$add_argument('-alt_a', type = 'character', help = 'Label of alternative allele column in input file A', default = 'ALT')
parser$add_argument('-alt_b', type = 'character', help = 'Label of alternative allele column in input file B', default = 'ALT')
parser$add_argument('-p_a', type = 'character', help = 'Label of p-value column in input file A', default = 'P')
parser$add_argument('-p_b', type = 'character', help = 'Label of p-value column in input file B', default = 'P')
parser$add_argument('-o', '--output_path', type = 'character', help = 'Path to output file', required = T)
parser$add_argument('-nt', '--no_of_threads', type = 'integer', help = 'Number of threads to use', default = 1)

args <- parser$parse_args()

setDTthreads(args$no_of_threads)

dat_a <- fread(args$input_file_a, sep = '\t', header = T, select = c(args$chr_a, args$bp_a, args$ref_a, args$alt_a, args$p_a))
dat_b <- fread(args$input_file_b, sep = '\t', header = T, select = c(args$chr_b, args$bp_b, args$ref_b, args$alt_b, args$p_b))

dat_a[ , c(args$ref_a, args$alt_a) := list(toupper(get(args$ref_a)), toupper(get(args$alt_a)))]
dat_a <- dat_a[get(args$ref_a) %in% c('A','T','C','G') & get(args$alt_a) %in% c('A','T','C','G')]
dat_a[, (args$chr_a) := as.character(get(args$chr_a))]
dat_a <- na.omit(dat_a)

dat_b[ , c(args$ref_b, args$alt_b) := list(toupper(get(args$ref_b)), toupper(get(args$alt_b)))]
dat_b <- dat_b[get(args$ref_b) %in% c('A','T','C','G') & get(args$alt_b) %in% c('A','T','C','G')]
dat_b[, (args$chr_b) := as.character(get(args$chr_b))]
dat_b <- na.omit(dat_b)

merged_dat <- merge(dat_a, dat_b, by.x = c(args$chr_a, args$bp_a), by.y = c(args$chr_b, args$bp_b), suffixes = c('', '.y'))

merged_dat <- merged_dat[!(get(args$chr_a) == 6 & get(args$bp_a) %in% c(24e6, 45e6))]

ref_a <- args$ref_a
alt_a <- args$alt_a
ref_b <- paste0(args$ref_b, '.y')
alt_b <- paste0(args$alt_b, '.y')
p_a <- args$p_a
p_b <- paste0(args$p_b, '.y')

merged_dat <- merged_dat[
  (get(ref_a) == get(ref_b) & get(alt_a) == get(alt_b)) |
      (get(ref_a) == get(alt_b) & get(alt_a) == get(ref_b))
]

names(merged_dat)[names(merged_dat) == p_a] <- 'P.A'
names(merged_dat)[names(merged_dat) == p_b] <- 'P.B'

merged_dat <- unique(merged_dat, cols = c('CHR38', 'BP38'))

merged_dat[, c(ref_b, alt_b) := NULL]

fwrite(merged_dat, file = args$output_path, sep = '\t', row.names = F)
