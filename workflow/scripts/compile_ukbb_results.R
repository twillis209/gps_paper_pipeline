library(data.table)

trait_metadata_dat <- fread('resources/ukbb_sum_stats/trait_metadata.tsv', sep = '\t', header = T)

ukbb_all_ldsc_dat <- fread('results/ldsc/rg/ukbb/all/fixed_h2_free_rg_intercept/compiled_results.tsv', sep = '\t', header = T)

ukbb_sans_mhc_ldsc_dat <- fread('results/ldsc/rg/ukbb/sans_mhc/fixed_h2_free_rg_intercept/compiled_results.tsv', sep = '\t', header = T)

ukbb_all_sumher_dat <- fread("results/ldak/ldak-thin/ukbb/all/rg/compiled_ukbb_sumher_results.tsv", sep = '\t', header = T)

ukbb_sans_mhc_sumher_dat <- fread("results/ldak/ldak-thin/ukbb/sans_mhc/rg/compiled_ukbb_sumher_results.tsv", sep = '\t', header = T)

ukbb_all_hoeffdings_dat <- fread("results/combined/all/window_1000kb_step_50/hoeffdings_results.tsv", sep = '\t', header = T)

ukbb_sans_mhc_hoeffdings_dat <- fread("results/combined/sans_mhc/window_1000kb_step_50/hoeffdings_results.tsv", sep = '\t', header = T)

ukbb_all_gps_dat <- fread("results/gps/combined/all/window_1000kb_step_50/gps_pvalues_3000_permutations.tsv", sep = '\t', header = T)

ukbb_sans_mhc_gps_dat <- fread("results/gps/combined/sans_mhc/window_1000kb_step_50/gps_pvalues_3000_permutations.tsv", sep = '\t', header = T)

setnames(ukbb_all_ldsc_dat, 
         c('h2.A.obs.ldsc', 'h2.A.obs.se.ldsc', 'h2.B.obs.ldsc', 'h2.B.obs.se.ldsc', 'gcov.obs.ldsc', 'gcov.obs.se.ldsc', 'rg.ldsc', 'rg.se.ldsc', 'rg.z.ldsc', 'rg.p.ldsc'),
         c('h2.A.obs.ldsc.all', 'h2.A.obs.se.ldsc.all', 'h2.B.obs.ldsc.all', 'h2.B.obs.se.ldsc.all', 'gcov.obs.ldsc.all', 'gcov.obs.se.ldsc.all', 'rg.ldsc.all', 'rg.se.ldsc.all', 'rg.z.ldsc.all', 'rg.p.ldsc.all')
)

setnames(ukbb_sans_mhc_ldsc_dat, 
         c('h2.A.obs.ldsc', 'h2.A.obs.se.ldsc', 'h2.B.obs.ldsc', 'h2.B.obs.se.ldsc', 'gcov.obs.ldsc', 'gcov.obs.se.ldsc', 'rg.ldsc', 'rg.se.ldsc', 'rg.z.ldsc', 'rg.p.ldsc'),
         c('h2.A.obs.ldsc.sans_mhc', 'h2.A.obs.se.ldsc.sans_mhc', 'h2.B.obs.ldsc.sans_mhc', 'h2.B.obs.se.ldsc.sans_mhc', 'gcov.obs.ldsc.sans_mhc', 'gcov.obs.se.ldsc.sans_mhc', 'rg.ldsc.sans_mhc', 'rg.se.ldsc.sans_mhc', 'rg.z.ldsc.sans_mhc', 'rg.p.ldsc.sans_mhc')
)

setnames(ukbb_all_sumher_dat, 
         c('h2.A.obs.sr', 'h2.A.obs.se.sr', 'h2.B.obs.sr', 'h2.B.obs.se.sr', 'gcov.obs.sr', 'gcov.obs.se.sr', 'rg.sr', 'rg.se.sr', 'rg.z.sr', 'rg.p.sr'),
         c('h2.A.obs.sr.all', 'h2.A.obs.se.sr.all', 'h2.B.obs.sr.all', 'h2.B.obs.se.sr.all', 'gcov.obs.sr.all', 'gcov.obs.se.sr.all', 'rg.sr.all', 'rg.se.sr.all', 'rg.z.sr.all', 'rg.p.sr.all')
         )

setnames(ukbb_sans_mhc_sumher_dat, 
         c('h2.A.obs.sr', 'h2.A.obs.se.sr', 'h2.B.obs.sr', 'h2.B.obs.se.sr', 'gcov.obs.sr', 'gcov.obs.se.sr', 'rg.sr', 'rg.se.sr', 'rg.z.sr', 'rg.p.sr'),
         c('h2.A.obs.sr.sans_mhc', 'h2.A.obs.se.sr.sans_mhc', 'h2.B.obs.sr.sans_mhc', 'h2.B.obs.se.sr.sans_mhc', 'gcov.obs.sr.sans_mhc', 'gcov.obs.se.sr.sans_mhc', 'rg.sr.sans_mhc', 'rg.se.sr.sans_mhc', 'rg.z.sr.sans_mhc', 'rg.p.sr.sans_mhc')
         )

setnames(ukbb_all_hoeffdings_dat, c('trait_A', 'trait_B', 'p.value'), c('trait.A', 'trait.B', 'hoeff.p.all'))

setnames(ukbb_sans_mhc_hoeffdings_dat, c('trait_A', 'trait_B', 'p.value'), c('trait.A', 'trait.B', 'hoeff.p.sans_mhc'))

setnames(ukbb_all_gps_dat, c('trait_A', 'trait_B', 'pval', 'gps'), c('trait.A', 'trait.B', 'gps.p.all', 'gps.all'))

setnames(ukbb_sans_mhc_gps_dat, c('trait_A', 'trait_B', 'pval', 'gps'), c('trait.A', 'trait.B', 'gps.p.sans_mhc', 'gps.sans_mhc'))

ukbb_all_ldsc_dat[, pair := paste(sort(c(trait.A, trait.B)), collapse = '-'), by = 1:nrow(ukbb_all_ldsc_dat)]
ukbb_sans_mhc_ldsc_dat[, pair := paste(sort(c(trait.A, trait.B)), collapse = '-'), by = 1:nrow(ukbb_sans_mhc_ldsc_dat)]
ukbb_all_sumher_dat[, pair := paste(sort(c(trait.A, trait.B)), collapse = '-'), by = 1:nrow(ukbb_all_sumher_dat)]
ukbb_sans_mhc_sumher_dat[, pair := paste(sort(c(trait.A, trait.B)), collapse = '-'), by = 1:nrow(ukbb_sans_mhc_sumher_dat)]
ukbb_all_hoeffdings_dat[, pair := paste(sort(c(trait.A, trait.B)), collapse = '-'), by = 1:nrow(ukbb_all_hoeffdings_dat)]
ukbb_sans_mhc_hoeffdings_dat[, pair := paste(sort(c(trait.A, trait.B)), collapse = '-'), by = 1:nrow(ukbb_sans_mhc_hoeffdings_dat)]
ukbb_all_gps_dat[, pair := paste(sort(c(trait.A, trait.B)), collapse = '-'), by = 1:nrow(ukbb_all_gps_dat)]
ukbb_sans_mhc_gps_dat[, pair := paste(sort(c(trait.A, trait.B)), collapse = '-'), by = 1:nrow(ukbb_sans_mhc_gps_dat)]

ukbb_all_ldsc_dat[, snp.set := NULL]
ukbb_sans_mhc_ldsc_dat[, snp.set := NULL]
ukbb_all_sumher_dat[, snp.set := NULL]
ukbb_sans_mhc_sumher_dat[, snp.set := NULL]

merge(ukbb_all_ldsc_dat,
      ukbb_sans_mhc_ldsc_dat[, !c('trait.A', 'trait.B')],
      by = 'pair', all = T) %>%
  merge(.,
        ukbb_all_sumher_dat[, !c('trait.A', 'trait.B')],
        by = 'pair', all = T) %>%
  merge(.,
        ukbb_sans_mhc_sumher_dat[, !c('trait.A', 'trait.B')],
        by = 'pair', all = T) %>%
  merge(.,
        ukbb_all_hoeffdings_dat[, !c('trait.A', 'trait.B', 'n', 'Dn', 'scaled')],
        by = 'pair', all = T) %>%
  merge(.,
        ukbb_sans_mhc_hoeffdings_dat[, !c('trait.A', 'trait.B', 'n', 'Dn', 'scaled')],
        by = 'pair', all = T) %>%
  merge(.,
        ukbb_all_gps_dat[, !c('trait.A', 'trait.B', 'n', 'loc', 'loc.sd', 'scale', 'scale.sd', 'shape', 'shape.sd')],
        by = 'pair', all = T) %>%
  merge(.,
        ukbb_sans_mhc_gps_dat[, !c('trait.A', 'trait.B', 'n', 'loc', 'loc.sd', 'scale', 'scale.sd', 'shape', 'shape.sd')],
        by = 'pair', all = T) -> merged_ukbb_dat

merged_ukbb_dat[ ,`:=` (trait.A = unlist(strsplit(pair, '-'))[1], trait.B = unlist(strsplit(pair, '-'))[2]), by = 1:nrow(merged_ukbb_dat)]

merge(merged_ukbb_dat[!is.na(trait.A)],
      trait_metadata_dat[, .(code, short_abbrv, long_abbrv, n_cases, n_controls, group)],
      by.x = 'trait.A', 
      by.y = 'code', allow.cartesian = T)[, 
                 setnames(.SD, c('short_abbrv', 'long_abbrv', 'n_cases', 'n_controls', 'group'), c('abbrv.A', 'long_abbrv.A', 'ncases.A', 'ncontrols.A', 'group.A'))
                 ][] -> int_res
  merge(x = int_res[!is.na(trait.B)], 
        y = trait_metadata_dat[, .(code, short_abbrv, long_abbrv, n_cases, n_controls, group)],
        by.x = 'trait.B',
        by.y = 'code', allow.cartesian = T)[, 
                 setnames(.SD, c('short_abbrv', 'long_abbrv', 'n_cases', 'n_controls', 'group'), c('abbrv.B', 'long_abbrv.B', 'ncases.B', 'ncontrols.B', 'group.B'))
                 ][] -> ukbb_with_meta_dat
  
ukbb_with_meta_dat[, ncases.ukbb.min := min(ncases.A, ncases.B), by = 1:nrow(ukbb_with_meta_dat)]
ukbb_with_meta_dat[, pair.long := paste(sort(c(long_abbrv.A, long_abbrv.B)), collapse = '-'), by = 1:nrow(ukbb_with_meta_dat)]
ukbb_with_meta_dat[, pair.abbrv := paste(sort(c(abbrv.A, abbrv.B)), collapse = '-'), by = 1:nrow(ukbb_with_meta_dat)]

chosen_traits_long_abbrv <- c("lupus",
"T1D",
"CD",
"UC",
"rheumatoid arthritis",
"eczema/derm",
"hypothyroidism",
"hayfever",
"asthma",
"cardiomyopathy",
"endometriosis",
"MD",
"glaucoma",
"leiomyoma",
"IBS",
"cholelithiasis",
"osteoarthritis",
"hypercholesterolaemia")

ukbb_with_meta_dat <- ukbb_with_meta_dat[long_abbrv.A %in% chosen_traits_long_abbrv & long_abbrv.B %in% chosen_traits_long_abbrv]

for(x in c('rg.p.ldsc.all', 'rg.p.ldsc.sans_mhc', 'rg.p.sr.all', 'rg.p.sr.sans_mhc', 'hoeff.p.all', 'hoeff.p.sans_mhc', 'gps.p.all', 'gps.p.sans_mhc')
) {
  ukbb_with_meta_dat[, pstat.trunc := ifelse(pstat < 1e-8, 1e-8, pstat), env = list(pstat = x, pstat.trunc = paste0(x, '.trunc'))]
}
