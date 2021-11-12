library(data.table)
setDTthreads(8)

sum_stats_files <- c("resources/ukbb_sum_stats/20002_1111.gwas.imputed_v3.both_sexes.tsv",
"resources/ukbb_sum_stats/20002_1113.gwas.imputed_v3.both_sexes.tsv",
"resources/ukbb_sum_stats/20002_1154.gwas.imputed_v3.both_sexes.tsv",
"resources/ukbb_sum_stats/20002_1220.gwas.imputed_v3.both_sexes.tsv",
"resources/ukbb_sum_stats/20002_1226.gwas.imputed_v3.both_sexes.tsv",
"resources/ukbb_sum_stats/20002_1286.gwas.imputed_v3.both_sexes.tsv",
"resources/ukbb_sum_stats/20002_1289.gwas.imputed_v3.both_sexes.tsv",
"resources/ukbb_sum_stats/20002_1381.gwas.imputed_v3.both_sexes.tsv",
"resources/ukbb_sum_stats/20002_1452.gwas.imputed_v3.both_sexes.tsv",
"resources/ukbb_sum_stats/20002_1462.gwas.imputed_v3.both_sexes.tsv",
"resources/ukbb_sum_stats/20002_1464.gwas.imputed_v3.both_sexes.tsv",
"resources/ukbb_sum_stats/20002_1465.gwas.imputed_v3.both_sexes.tsv",
"resources/ukbb_sum_stats/20002_1473.gwas.imputed_v3.both_sexes.tsv",
"resources/ukbb_sum_stats/22126.gwas.imputed_v3.both_sexes.tsv",
"resources/ukbb_sum_stats/6148_2.gwas.imputed_v3.both_sexes.tsv",
"resources/ukbb_sum_stats/6148_5.gwas.imputed_v3.both_sexes.tsv",
"resources/ukbb_sum_stats/D25.gwas.imputed_v3.both_sexes.tsv",
"resources/ukbb_sum_stats/I9_IHD.gwas.imputed_v3.both_sexes.tsv",
"resources/ukbb_sum_stats/K51.gwas.imputed_v3.both_sexes.tsv",
"resources/ukbb_sum_stats/K57.gwas.imputed_v3.both_sexes.tsv",
"resources/ukbb_sum_stats/K80.gwas.imputed_v3.both_sexes.tsv")

names(sum_stats_files) <-
  c("20002_1111",
"20002_1113",
"20002_1154",
"20002_1220",
"20002_1226",
"20002_1286",
"20002_1289",
"20002_1381",
"20002_1452",
"20002_1462",
"20002_1464",
"20002_1465",
"20002_1473",
"22126",
"6148_2",
"6148_5",
"D25",
"I9_IHD",
"K51",
"K57",
"K80")

left_dat <- fread(sum_stats_files[1], sep = '\t', header = T, select = c('variant', 'pval'))

setnames(left_dat, c('variant', 'pval'), c('variant', paste0('pval.', names(sum_stats_files[1]))))

for(i in 2:length(sum_stats_files)) {
  right_dat <- fread(sum_stats_files[i], sep = '\t', header = T, select = c('variant', 'pval'))

  setnames(right_dat, c('variant', 'pval'), c('variant', paste0('pval.', names(sum_stats_files[i]))))

  left_dat <- merge(left_dat, right_dat, by = 'variant')
}

fwrite(left_dat, file = 'resources/ukbb_sum_stats/merged_ukbb_sum_stats.tsv.gz', sep = '\t')
