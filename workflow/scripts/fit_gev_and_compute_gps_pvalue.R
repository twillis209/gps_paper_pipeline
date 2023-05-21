library(data.table)
library(evd)
library(fitdistrplus)

gps_dat <- fread(snakemake@input[['gps_file']], sep = '\t', header = T)

gps <- gps_dat[Trait_A == snakemake@params[['trait_A']] & Trait_B == snakemake@params[['trait_B']], GPS]

perm_dat <- fread(snakemake@input[['perm_file']], sep = '\t', header = T)

if(perm_dat[is.infinite(GPS), .N] > 0) {
  stop(sprintf("%d infinite-valued permuted GPS realisations", perm_dat[is.infinite(GPS), .N]))
}

fgev.fit <- tryCatch(
   fgev(perm_dat$GPS),
  error = function(c) {
    msg <- conditionMessage(c)
    if(msg == "observed information matrix is singular; use std.err = FALSE"){
      fgev(perm_dat$GPS, std.err = F)
    } else {
      stop(msg)
      }
    }
)

fgev.fitdist <- fitdist(perm_dat$GPS, 'gev', start = list(loc = fgev.fit$estimate[['loc']], scale = fgev.fit$estimate[['scale']], shape = fgev.fit$estimate[['shape']]))

loc <- fgev.fitdist$estimate[['loc']]
loc.sd <- fgev.fitdist$sd[['loc']]
scale <- fgev.fitdist$estimate[['scale']]
scale.sd <- fgev.fitdist$sd[['scale']]
shape <- fgev.fitdist$estimate[['shape']]
shape.sd <- fgev.fitdist$sd[['shape']]

pvalue <- pgev(gps, loc = loc, scale = scale, shape = shape, lower.tail = F)

fwrite(data.table(gps = gps, n = nrow(perm_dat), loc = loc, loc.sd = loc.sd, scale = scale, scale.sd = scale.sd, shape = shape, shape.sd = shape.sd, pval = pvalue), sep = '\t', file = snakemake@output[[1]])
