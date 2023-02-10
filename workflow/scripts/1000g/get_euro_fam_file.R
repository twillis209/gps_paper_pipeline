panel <- read.table(snakemake@input[['panel_file']], header = T)
ped <- read.table(snakemake@input[['ped_file']], header = T, sep = '\t')

merged <- merge(panel, ped[c('Individual.ID', 'Paternal.ID', 'Maternal.ID')], by.x = 'sample', by.y = 'Individual.ID', all.x = T)

# Get unrelated European samples
euro <- subset(merged, super_pop == 'EUR' & Paternal.ID == 0 & Maternal.ID == 0)

euro <- euro[c('sample', 'sample', 'Paternal.ID', 'Maternal.ID', 'gender')]

names(euro) <- c('SampleID', 'SampleID', 'FatherID', 'MotherID', 'Sex')

# Fix to work with keep as implemented in plink2
write.table(euro[, c('SampleID')], file = snakemake@output[[1]], sep = ' ', col.names = F, row.names = F, quote = F)
