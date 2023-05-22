library(data.table)
setDTthreads(snakemake@threads)

overlap <- fread(snakemake@input[['overlap']], sep = '\t', header = T)
metadata <- fread(snakemake@input[['metadata']], sep = '\t', header = T)

overlap <- overlap[trait_A %in% snakemake@params[['traits_to_keep']] & trait_B %in% snakemake@params[['traits_to_keep']]]

m_dat <- merge(overlap, metadata[, .(code, desc_A = long_abbrv, ncases_A = n_cases, ncontrols_A = n_controls)], by.x = 'trait_A', by.y = 'code')
m_dat <- merge(m_dat, metadata[, .(code, desc_B = long_abbrv, ncases_B = n_cases, ncontrols_B = n_controls)], by.x = 'trait_B', by.y = 'code')

cols_to_convert <- c("A_controls", "B_controls", "A_cases", "B_cases", "AB_controls", "AB_cases")
m_dat[, c(cols_to_convert) := lapply(.SD, as.numeric), .SDcols = cols_to_convert]

m_dat[, f.temp := sqrt(A_cases * B_cases / A_controls / B_controls)]
m_dat[, rho := ( AB_controls * f.temp + AB_cases / f.temp ) / sqrt( (A_controls + A_cases) * (B_controls + B_cases) )]
m_dat[, f.temp := NULL]

fwrite(m_dat, file = snakemake@output[[1]], sep = '\t')
