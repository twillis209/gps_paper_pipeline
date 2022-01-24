library(argparse)
library(data.table)
library(evd)
library(fitdistrplus)

parser <- ArgumentParser(description = 'Fits GEV to permuted data using increasingly large subset.')

parser$add_argument('-a', '--trait_A', type = 'character', help = 'Trait A label', required = T)
parser$add_argument('-b', '--trait_B', type = 'character', help = 'Trait B label', required = T)
parser$add_argument('-i', '--perm_file', type = 'character', help = 'Path to file containing GPS values generated from permuted data', required = T)
parser$add_argument('-o', '--output_file', type = 'character', help = 'Path to output file', required = T)

args <- parser$parse_args()

perm_dat <- fread(args$perm_file, sep = '\t', header = T)

fgev.fit <- fgev(perm_dat$GPS)

fgev.fitdist <- fitdist(perm_dat$GPS, 'gev', start = list(loc = fgev.fit$estimate[['loc']], scale = fgev.fit$estimate[['scale']], shape = fgev.fit$estimate[['shape']]))

loc <- fgev.fitdist$estimate[['loc']]
loc.sd <- fgev.fitdist$sd[['loc']]
scale <- fgev.fitdist$estimate[['scale']]
scale.sd <- fgev.fitdist$sd[['scale']]
shape <- fgev.fitdist$estimate[['shape']]
shape.sd <- fgev.fitdist$sd[['shape']]

estimate_dat <- data.table(trait_A = args$trait_A, trait_B = args$trait_B, n = nrow(perm_dat), loc = loc, loc.sd = loc.sd, scale = scale, scale.sd = scale.sd, shape = shape, shape.sd = shape.sd)

fwrite(estimate_dat, sep = '\t', file = args$output_file)
