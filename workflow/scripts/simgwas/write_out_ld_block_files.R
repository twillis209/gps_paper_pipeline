library(data.table)
library(argparse)

parser <- ArgumentParser(description = 'Chops legend and haplotype files into LD block-specific files')
parser$add_argument('--hap_file', type = 'character', help = 'Path to haplotype file')
parser$add_argument('--leg_file', type = 'character', help = 'Path to legend file')
parser$add_argument('-b', '--block_file', type = 'character', help = 'Path to block file')
parser$add_argument('--chr_no', type = 'integer', help = 'Number of chromosome')
parser$add_argument('-o', '--output_root', type = 'character', help = 'Path to output directory', required = T)
parser$add_argument('-nt', '--no_of_threads', type = 'integer', help = 'Number of threads to use', default = 1)

args <- parser$parse_args()

setDTthreads(args$no_of_threads)

leg_dat <- fread(file = args$leg_file, sep = ' ', header = T)
# Skip leading metadata rows
hap_dat <- fread(file = args$hap_file, sep = ' ', header = F, skip = 4)
block_dat <- fread(file = args$block_file, sep = ' ', header = F, col.names = c('block', 'chr', 'start', 'stop'))

block_dat <- block_dat[chr == args$chr_no]

leg_dat[, rs := make.names(id)]

hap_dat[, rs := leg_dat$rs]

for(i in 1:nrow(block_dat)) {
   leg_dat[position %between% c(block_dat[i, start], block_dat[i, stop]), block := (i-1)]
}

hap_dat <- hap_dat[rs %in% leg_dat$rs]

hap_dat[, block := leg_dat$block]

hap_dat <- hap_dat[, c(seq(1, ncol(hap_dat)-3, by = 2), seq(2, ncol(hap_dat)-2, by = 2), ncol(hap_dat)-1, ncol(hap_dat)), with = F]

for(i in sort(unique(hap_dat$block))) {
  fwrite(leg_dat[block == i], file = file.path(args$output_root, sprintf('block_%s.legend.gz', i)), sep = ' ')
  fwrite(hap_dat[block == i], file = file.path(args$output_root, sprintf('block_%s.hap.gz', i)), sep = ' ', col.names = F)
}
