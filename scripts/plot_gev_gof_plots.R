library(argparse)
library(data.table)
library(fitdistrplus)
library(evd)

parser <- ArgumentParser(description = 'Plot goodness-of-fit plots for GEV fit to permuted data')
parser$add_argument('-f', '--fitdist_file', type = 'character', help = 'Path to fitted parameters file')
parser$add_argument('-p', '--perm_file', type = 'character', help = 'Path to file contained permuted null GPS statistics')
parser$add_argument('-l', '--trait_pair_label', type = 'character', help = 'Trait pair label')
parser$add_argument('-o', '--output_path', type = 'character', help = 'Path to output plot file', required = T)

args <- parser$parse_args()

perm_sample <- scan(args$perm_file, skip = 1)

fit_dat <- fread(args$fitdist_file, sep = '\t', header = T)

fgev.fitdist <- fitdist(perm_sample, 'gev', start = list(loc = fit_dat$loc, scale = fit_dat$scale, shape = fit_dat$shape))

fgev.gof <- gofstat(fgev.fitdist)

png(args$output_path)
plot(fgev.fitdist)
title(main = sprintf("%s GOF plots", args$trait_pair_label), outer = T, line = -1)
dev.off()
