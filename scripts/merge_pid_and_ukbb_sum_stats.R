library(data.table)
setDTthreads(8)

ukbb_dat <- fread("resources/ukbb_sum_stats/merged_ukbb_sum_stats.tsv.gz", sep = '\t', header = T)

pid_dat <- fread("resources/pid.tsv.gz", sep = '\t', header = T, select = c('CHR19', 'BP19', 'REF', 'ALT', 'P'))
pid_dat[, CHR19 := as.character(CHR19)]
pid_dat[, BP19 := as.integer(BP19)]

ukbb_dat[, c("CHR19", "BP19", "ALT", "REF") := tstrsplit(variant, split = ':', keep = 1:4)]
ukbb_dat[, CHR19 := as.character(CHR19)]
ukbb_dat[, BP19 := as.integer(BP19)]

merged_dat <- merge(pid_dat, ukbb_dat, by = c("CHR19", "BP19"), all = T, sort = F)

merged_dat <- merged_dat[(REF.x == REF.y & ALT.x == ALT.y) | (REF.x == ALT.y & REF.y == ALT.x)]
setnames(merged_dat, c("REF.x", "REF.y", "ALT.x", "ALT.y", "P"), c("REF.pid", "REF.ukbb", "ALT.pid", "ALT.ukbb", "pval.pid"))

merged_dat <- merged_dat[!is.na(pval.pid) & !is.na(pval.20002_1111)]

fwrite(merged_dat, file = "resources/pid_ukbb_sum_stats/merged_pid_ukbb_sum_stats.tsv.gz", sep = '\t', col.names = T, row.names = F, quote = F)
