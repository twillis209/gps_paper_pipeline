library(data.table)
library(argparse)

parser <- ArgumentParser(description = 'Merge UKBB and PID summary statistics')
parser$add_argument('-u', '--ukbb_files', type = 'character', nargs = '+', help = 'List of UKBB summary statistics file paths')
parser$add_argument('-p', '--pid_file', type = 'character', help = 'PID GWAS summary statistics file path')
parser$add_argument('-o', '--output_path', type = 'character', help = 'Path to output file', required = T)
parser$add_argument('-nt', '--no_of_threads', type = 'integer', help = 'Number of threads to use', default = 1)

args <- parser$parse_args()

setDTthreads(args$no_of_threads)

args$ukbb_files <- unlist(strsplit(args$ukbb_files, ","))

pid_dat <- fread(args$pid_file, sep = '\t', header = T, select = c('CHR19', 'BP19', 'REF', 'ALT', 'P'))
pid_dat[, CHR19 := as.character(CHR19)]
pid_dat[, BP19 := as.integer(BP19)]
setnames(pid_dat, "P", "pval.pid")

# Get names from ukbb files by pattern-matching
ukbb_ids <- unlist(stringr::str_match_all(args$ukbb_files, "[[:alpha:]]?[[:digit:]]+_?[[:digit:][:alpha:]]+"))

left_dat <- fread(args$ukbb_files[1], sep = '\t', header = T, select = c('variant', 'pval'))

setnames(left_dat, 'pval', paste0('pval.', ukbb_ids[1]))

for(i in 2:length(args$ukbb_files)) {
  ukbb_dat <- fread(args$ukbb_files[i], sep = '\t', header = T, select = c('variant', 'pval'))

  setnames(ukbb_dat, 'pval', paste0('pval.', ukbb_ids[i]))

  left_dat <- merge(left_dat, ukbb_dat, by = 'variant')
}

left_dat[, c("CHR19", "BP19", "ALT", "REF") := tstrsplit(variant, split = ':', keep = 1:4)]
left_dat[, CHR19 := as.character(CHR19)]
left_dat[, BP19 := as.integer(BP19)]

merged_dat <- merge(pid_dat, left_dat, by = c("CHR19", "BP19"), all = T, sort = F)

merged_dat <- merged_dat[(REF.x == REF.y & ALT.x == ALT.y) | (REF.x == ALT.y & REF.y == ALT.x)]
setnames(merged_dat, c("REF.x", "REF.y", "ALT.x", "ALT.y", "P"), c("REF.pid", "REF.ukbb", "ALT.pid", "ALT.ukbb", "pval.pid"), skip_absent = T)

merged_dat <- merged_dat[!is.na(pval.pid) & !is.na(get(paste0('pval.', ukbb_ids[1])))]

fwrite(merged_dat, file = "resources/merged_pid_ukbb_sum_stats.tsv.gz", sep = '\t', col.names = T, row.names = F, quote = F)
