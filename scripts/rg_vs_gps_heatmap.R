library(data.table)
library(argparse)
library(pheatmap)

parser <- ArgumentParser(description = 'Plot rg vs. gps heatmap')
parser$add_argument('-p', '--pvalue_file', type = 'character', help = 'P-value file path')
parser$add_argument('-r', '--rg_file', type = 'character', help = 'rg file path')
parser$add_argument('-t', '--traits', type = 'character', help = 'Comma-delimited list of traits in desired order')
parser$add_argument('-o', '--output_path', type = 'character', help = 'Path to output file', required = T)

args <- parser$parse_args()

args$traits <- unlist(strsplit(args$traits, ','))

pvalue_dat <- fread(args$pvalue_file, sep = '\t', header = T)

rg_dat <- fread(args$rg_file, sep = '\t', header = T)

gps_rg_mat <- matrix(nrow = length(args$traits), ncol = length(args$traits))
colnames(gps_rg_mat) <- args$traits
rownames(gps_rg_mat) <- args$traits


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

neg_log10_gps_rg_mat <- -log10(gps_rg_mat)

pheatmap(apply(apply(neg_log10_gps_rg_mat, 2, cut, breaks = seq(0,10, by = 2), labels = c(0,2,4,6,8)), 2, as.integer), breaks = seq(0, 10, length.out = 101), cluster_rows = F, cluster_cols = F, fontsize = 13, filename = args$output_path, display_numbers = F, number_format = "%.1f", labels_row = rownames(neg_log10_gps_rg_mat), labels_col = colnames(neg_log10_gps_rg_mat))
