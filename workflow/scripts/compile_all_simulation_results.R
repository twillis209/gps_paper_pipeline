library(data.table)
library(magrittr)

effect_blocks <- c('s400', 'm25', 'm50', 's200-m25')

merged_dats <- list()

for(i in seq_along(effect_blocks)) {
  x <- effect_blocks[i]
  
  ldsc_dat <- fread(sprintf('results/ldsc/simgwas/400_reps/randomised/compiled_%s_results.tsv', x), sep = '\t', header = T)

  ldsc_dat <- ldsc_dat[!is.na(h2.A)]

  sumher_dat <- fread(sprintf('results/ldak/ldak-thin/simgwas/400_reps/randomised/rg/compiled_%s_results.tsv', x), sep = '\t', header = T)

  hoeffdings_dat <- fread(sprintf("results/hoeffdings/simgwas/400_reps/randomised/compiled_%s_results.tsv", x), sep = '\t', header = T)
  
  gps_dat <- fread(sprintf("results/gps/simgwas/400_reps/randomised/compiled_%s_results.tsv", x), sep = '\t', header = T)
  
  theo_dat <- fread(sprintf("results/ldsc/simgwas/400_reps/randomised/compiled_%s_theo_rg.tsv", x), sep = '\t', header = T)

  setnames(ldsc_dat, c('h2.A', 'h2.A.se', 'h2.B', 'h2.B.se', 'gcov', 'gcov.se', 'rg', 'rg.se', 'rg.z', 'rg.p'),
          c('h2.A.ldsc', 'h2.A.se.ldsc', 'h2.B.ldsc', 'h2.B.se.ldsc', 'gcov.ldsc', 'gcov.se.ldsc', 'rg.ldsc', 'rg.se.ldsc', 'rg.z.ldsc', 'rg.p.ldsc'))

  setnames(sumher_dat, c('h2.A', 'h2.A.se', 'h2.B', 'h2.B.se', 'gcov', 'gcov.se', 'rg', 'rg.se', 'rg.z', 'rg.p'),
          c('h2.A.sr', 'h2.A.se.sr', 'h2.B.sr', 'h2.B.se.sr', 'gcov.sr', 'gcov.se.sr', 'rg.sr', 'rg.se.sr', 'rg.z.sr', 'rg.p.sr'))

  setnames(gps_dat, 'pval', 'gps.p')

  theo_dat[, r_A.AB.mean := mean(r_A.AB), by = .(ncases.A, ncontrols.A, ncases.B, ncontrols.B, shared_blocks)]

  merge_keys <- c('ncases.A', 'ncontrols.A', 'ncases.B', 'ncontrols.B', 'shared_blocks', 'tag_pair', 'seed')

  merge(ldsc_dat[, .(ncases.A, ncontrols.A, ncases.B, ncontrols.B, shared_blocks, blocks.A, tag_pair, seed, h2.A.ldsc, h2.A.se.ldsc, h2.B.ldsc, h2.B.se.ldsc, gcov.ldsc, gcov.se.ldsc, rg.ldsc, rg.se.ldsc, rg.z.ldsc, rg.p.ldsc)],
        sumher_dat[, .(ncases.A, ncontrols.A, ncases.B, ncontrols.B, shared_blocks, tag_pair, seed, h2.A.sr, h2.A.se.sr, h2.B.sr, h2.B.se.sr, gcov.sr, gcov.se.sr, rg.sr, rg.se.sr, rg.z.sr, rg.p.sr)],
        by = merge_keys, all = T) %>%
    merge(.,
          hoeffdings_dat[, .(ncases.A, ncontrols.A, ncases.B, ncontrols.B, shared_blocks, tag_pair, seed, hoeff.p)],
          by = merge_keys, all = T) %>%
    merge(.,
          gps_dat[, .(ncases.A, ncontrols.A, ncases.B, ncontrols.B, shared_blocks, tag_pair, seed, gps, gps.p)],
          by = merge_keys, all = T) %>%
    merge(.,
          theo_dat[, .(ncases.A, ncontrols.A, ncases.B, ncontrols.B, shared_blocks, tag_pair, seed, r_A.AB, r_A.AB.mean)],
          by = merge_keys, all = T) -> merged_dats[[i]]
}

merged_dat <- rbindlist(merged_dats)

merged_dat <- merged_dat[!(blocks.A == 's400' & !(shared_blocks %in% c('s0', 's50', 's100', 's150', 's200', 's250', 's300', 's350', 's400')))]

merged_dat[blocks.A == 's400', shared_blocks := factor(shared_blocks, levels = c('s0', 's50', 's100', 's150', 's200', 's250', 's300', 's350', 's400'), ordered = T)]

merged_dat <- merged_dat[!(blocks.A == 'm25' & !(shared_blocks %in% c('m0', 'm5', 'm10', 'm15', 'm20', 'm25')))]

merged_dat[blocks.A == 'm25', shared_blocks := factor(shared_blocks, levels = c('m0', 'm2', 'm5', 'm7', 'm10', 'm12', 'm15', 'm17', 'm20', 'm25'), ordered = T)]

merged_dat <- merged_dat[!(blocks.A == 'm50' & !(shared_blocks %in% c('m0', 'm10', 'm20', 'm30', 'm40', 'm50')))]

merged_dat[blocks.A == 'm50', shared_blocks := factor(shared_blocks, levels = c('m0', 'm2', 'm5', 'm7', 'm10', 'm12', 'm15', 'm17', 'm20', 'm30', 'm40', 'm50'), ordered = T)]

merged_dat <- merged_dat[!(blocks.A == 's200-m25' & !(shared_blocks %in% c('s0-m0', 's100-m0', 's100-m15', 's100-m25', 's200-m0', 's200-m15', 's200-m25')))]

merged_dat[blocks.A == 's200-m25', shared_blocks := factor(shared_blocks, levels = c('s0-m0', 's100-m0', 's100-m15', 's100-m25', 's200-m0', 's200-m15', 's200-m25'), ordered = T)]

for(x in c('rg.p.ldsc', 'rg.p.sr', 'hoeff.p', 'gps.p')) {
  merged_dat[, pstat.trunc := ifelse(pstat < 1e-8, 1e-8, pstat), env = list(pstat = x, pstat.trunc = paste0(x, '.trunc'))]
}

fwrite(merged_dat, file = "results/all_sim.tsv", sep = '\t')
