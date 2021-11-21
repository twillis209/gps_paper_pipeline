library(data.table)
library(argparse)

parser <- ArgumentParser(description = 'Plot rg vs. gps heatmap')
parser$add_argument('-p', '--pvalue_file', type = 'character', help = 'P-value file path')
parser$add_argument('-r', '--rg_file', type = 'character', help = 'rg file path')
parser$add_argument('-t', '--traits', type = 'character', nargs = '+', help = 'Traits in desired order')
parser$add_argument('-o', '--output_path', type = 'character', help = 'Path to output file', required = T)

args <- parser$parse_args()

pvalue_dat <- fread(args$pvalue_file, sep = '\t', header = T)

rg_dat <- fread(args$rg_file, sep = '\t', header = T)

gps_rg_mat <- matrix(nrow = length(args$traits), ncol = length(args$traits))

for(i in 1:length(args$traits)) {
  for(j in i:length(args$traits)) {
    if(i == j) {
      gps_pval <- NA
      rg_pval <- NA
    } else {
      gps_pval <- pvalue_dat[(abbrv_A == args$traits[i] & abbrv_B == args$traits[j]) | (abbrv_A == args$traits[j] & abbrv_B == args$traits[i]), pval]
      rg_pval <- rg_dat[(trait_A == args$traits[i] & trait_B == args$traits[j]) | (trait_A == args$traits[j] & trait_B == args$traits[i]), pval]
    }

    gps_pval <- ifelse(length(gps_pval) > 0, gps_pval, NA)
    rg_pval <- ifelse(length(rg_pval) > 0, rg_pval, NA)

    gps_rg_mat[i, j] <- gps_pval
    gps_rg_mat[j, i] <- rg_pval
  }
}

save(gps_rg_mat, file = args$output_path)
