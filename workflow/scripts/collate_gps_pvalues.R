library(magrittr)

results_files <- list.files('results', pattern = '*_gps_pvalue.tsv', full.names = T)

trait_pairs <- gsub("results/(.*)_gps_pvalue.tsv", "\\1", results_files)

strsplit(trait_pairs, "-") %>%
  do.call(rbind, .) %>%
  data.frame %>%
  setNames(., c("Trait_A", "Trait_B")) -> trait_daf

results_files %>%
  lapply(., function(x) read.table(x, sep = '\t', header = T)) %>%
  do.call(rbind, .) %>%
  cbind(., trait_daf) -> daf

trait_name_daf <- read.table('resources/ukbb_sum_stats/traits_codes_abbrv_cases.csv', sep = ',', header = T)

imd_abbrvs <- c('asthma', 'diabetes', 'hypothyroidism', 'eczema/derm', 'CD', 'rheumatoid arthritis', 'hayfever', 'UC')

trait_A_values <- unique(daf$Trait_A)
trait_B_values <- unique(daf$Trait_B)

trait_values <- trait_name_daf$code

gps_mat <- matrix(nrow = length(trait_values), ncol = length(trait_values))
pval_mat <- matrix(nrow = length(trait_values), ncol = length(trait_values))
gps_pval_mat <- matrix(nrow = length(trait_values), ncol = length(trait_values))

rownames(gps_mat) <- trait_values
colnames(gps_mat) <- trait_values
rownames(pval_mat) <- trait_values
colnames(pval_mat) <- trait_values

for(i in 1:length(trait_values)) {
  for(j in i:length(trait_values)) {
    if(i == j) {
      gps <- NA
      pval <- NA
    } else {
      gps <- subset(daf, Trait_A == trait_values[i] & Trait_B == trait_values[j] | Trait_A == trait_values[j] & Trait_B == trait_values[i])$gps
      pval <- subset(daf, Trait_A == trait_values[i] & Trait_B == trait_values[j] | Trait_A == trait_values[j] & Trait_B == trait_values[i])$pval
    }

    gps_mat[i, j] <- gps
    gps_mat[j, i] <- gps
    pval_mat[i, j] <- pval
    pval_mat[j, i] <- pval
  }
}

rownames(gps_mat) <- plyr::mapvalues(rownames(gps_mat), from = trait_name_daf$code, to = trait_name_daf$abbrv)
colnames(gps_mat) <- plyr::mapvalues(colnames(gps_mat), from = trait_name_daf$code, to = trait_name_daf$abbrv)
rownames(pval_mat) <- plyr::mapvalues(rownames(pval_mat), from = trait_name_daf$code, to = trait_name_daf$abbrv)
colnames(pval_mat) <- plyr::mapvalues(colnames(pval_mat), from = trait_name_daf$code, to = trait_name_daf$abbrv)

imd_gps_mat <- gps_mat[rownames(gps_mat) %in% imd_abbrvs, colnames(gps_mat) %in% imd_abbrvs]
imd_pval_mat <- pval_mat[rownames(pval_mat) %in% imd_abbrvs, colnames(pval_mat) %in% imd_abbrvs]


#paste(format(pval_mat, scientific = T, digits = 2), format(gps_mat, digits = 2), sep = '/') %>% stringr::str_replace(., ' ', '')

rg_daf <- read.table('resources/ukbb_sum_stats/traits_rg.tsv', sep = '\t', header = T)

rg_daf$Trait_A <- plyr::mapvalues(rg_daf$p1_desc, trait_name_daf$abbrv, trait_name_daf$code)
rg_daf$Trait_B <- plyr::mapvalues(rg_daf$p2_desc, trait_name_daf$abbrv, trait_name_daf$code)

merged_daf <- merge(daf[c('gps', 'pval', 'Trait_A', 'Trait_B')], rg_daf[c('rg', 'se', 'p', 'Trait_A', 'Trait_B', 'p1_desc', 'p2_desc')], by = c('Trait_A', 'Trait_B'))

merged_daf <- merged_daf[c('p1_desc', 'p2_desc', 'gps', 'pval', 'rg', 'se', 'p')]

names(merged_daf) <- c('Trait_A', 'Trait_B', 'gps', 'pval.gps', 'rg', 'se.rg', 'pval.rg')

merged_daf <- merged_daf[order(merged_daf$pval.gps),]

merged_daf <- merge(merged_daf, trait_name_daf[c('abbrv', 'n_cases')], by.x = 'Trait_A', by.y = 'abbrv')
names(merged_daf)[names(merged_daf) == 'n_cases'] <- 'n_cases.A'

merged_daf <- merge(merged_daf, trait_name_daf[c('abbrv', 'n_cases')], by.x = 'Trait_B', by.y = 'abbrv')
names(merged_daf)[names(merged_daf) == 'n_cases'] <- 'n_cases.B'

merged_daf <- merged_daf[order(merged_daf$gps, decreasing = T),]

merged_daf <- merged_daf[c("Trait_A","Trait_B","gps","pval.gps","rg","se.rg","pval.rg","n_cases.A","n_cases.B")]

write.table(merged_daf, file = 'gps_rg_results.tsv', sep = '\t', row.names = F)

#merged_daf$gps <- format(merged_daf$gps, digits = 1)
#merged_daf$pval.gps <- format(merged_daf$pval.gps, digits = 2)
#merged_daf$rg <- signif(merged_daf$rg, digits = 2)
#
#merged_daf$n_cases.A <- format(merged_daf$n_cases.A, big.mark = ',')
#merged_daf$n_cases.B <- format(merged_daf$n_cases.B, big.mark = ',')

pretty_merged_daf <- merged_daf
pretty_merged_daf$pval.gps <- signif(pretty_merged_daf$pval.gps, 2)
pretty_merged_daf$gps <- round(pretty_merged_daf$gps, 1)
pretty_merged_daf[c('n_cases.A', 'n_cases.B')] <- format(pretty_merged_daf[c('n_cases.A', 'n_cases.B')], big.mark = ',')

write.table(pretty_merged_daf, file = 'pretty_format_gps_rg_results.tsv', sep = '\t', row.names = F)


#daf <- daf[with(daf, order(-gps)),]

#daf$Trait_A_abbrv <- plyr::mapvalues(daf$Trait_A, from = trait_name_daf$code, to = trait_name_daf$abbrv)
#daf$Trait_B_abbrv <- plyr::mapvalues(daf$Trait_B, from = trait_name_daf$code, to = trait_name_daf$abbrv)
