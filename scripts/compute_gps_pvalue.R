library(argparse)
library(data.table)
library(evd)
library(fitdistrplus)
setDTthreads(8)

parser <- ArgumentParser(description = 'Computes p-value for GPS statistic using permuted data')
parser$add_argument('-g', '--gps_file', type = 'character', help = 'Path to file containing GPS value', required = T)
parser$add_argument('-p', '--perm_file', type = 'character', help = 'Path to file containing GPS value generated from permuted data', required = T)
parser$add_argument('-a', '--trait_a', type = 'character', help = 'Trait A', required = T)
parser$add_argument('-b', '--trait_b', type = 'character', help = 'Trait B', required = T)
parser$add_argument('-o', '--output_file', type = 'character', help = 'Path to output file', required = T)

args <- parser$parse_args()

gps_dat <- fread(args$gps_file, sep = '\t', header = T)

gps <- gps_dat[Trait_A == args$trait_a & Trait_B == args$trait_b, GPS]

perm_dat <- fread(args$perm_file, sep = '\t', header = T)

fgev.fit <- fgev(perm_dat$GPS)

fgev.fitdist <- fitdist(perm_dat$GPS, 'gev', start = list(loc = fgev.fit$estimate[['loc']], scale = fgev.fit$estimate[['scale']], shape = fgev.fit$estimate[['shape']]))

loc <- fgev.fitdist$estimate[['loc']]
loc.sd <- fgev.fitdist$sd[['loc']]
scale <- fgev.fitdist$estimate[['scale']]
scale.sd <- fgev.fitdist$sd[['scale']]
shape <- fgev.fitdist$estimate[['shape']]
shape.sd <- fgev.fitdist$sd[['shape']]

pvalue <- pgev(gps, loc = loc, scale = scale, shape = shape, lower.tail = F)

fwrite(data.table(gps = gps, n = nrow(perm_dat), loc = loc, loc.sd = loc.sd, scale = scale, scale.sd = scale.sd, shape = shape, shape.sd = shape.sd, pval = pvalue), sep = '\t', file = args$output_file)
