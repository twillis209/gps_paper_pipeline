library(argparse)

parser <- ArgumentParser(description = 'Produce PLINK-compatible ID file from 1000G population metadata')
parser$add_argument('-pa', '--panel_file', type = 'character', help = 'Path to panel file', required = T)
parser$add_argument('-pe', '--ped_file', type = 'character', help = 'Path to ped file', required = T)
parser$add_argument('-o', '--output_file', type = 'character', help = 'Path to fam file', required = T)

args <- parser$parse_args()

panel <- read.table(args$panel_file, header = T)
ped <- read.table(args$ped_file, header = T, sep = '\t')

merged <- merge(panel, ped[c('Individual.ID', 'Paternal.ID', 'Maternal.ID')], by.x = 'sample', by.y = 'Individual.ID', all.x = T)

# Get unrelated European samples
euro <- subset(merged, super_pop == 'EUR' & Paternal.ID == 0 & Maternal.ID == 0)

euro <- euro[c('sample', 'sample', 'Paternal.ID', 'Maternal.ID', 'gender')]

names(euro) <- c('SampleID', 'SampleID', 'FatherID', 'MotherID', 'Sex')

# Fix to work with keep as implemented in plink2
write.table(euro[, c('SampleID')], file = args$output_file, sep = ' ', col.names = F, row.names = F, quote = F)
