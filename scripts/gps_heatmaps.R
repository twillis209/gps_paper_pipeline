library(data.table)
library(argparse)
library(pheatmap)

parser <- ArgumentParser(description = 'Plot rg vs. gps heatmap')
parser$add_argument('-i', '--gps_file', type = 'character', help = 'gps file path')
parser$add_argument('-t', '--traits', type = 'character', help = 'Comma-delimited list of traits in desired order')
parser$add_argument('-o', '--output', type = 'character', help = 'Path to output file', required = T)

args <- parser$parse_args()

args$traits <- unlist(strsplit(args$traits, ','))

gps_dat <- fread(args$gps_file, sep = '\t', header = T)

gps_mat <- matrix(nrow = length(args$traits), ncol = length(args$traits))
colnames(gps_mat) <- args$traits
rownames(gps_mat) <- args$traits

for(i in 1:length(args$traits)) {
  for(j in i:length(args$traits)) {
    if(i == j) {
      gps <- NA
    } else {
      gps <- gps_dat[(abbrv_A == args$traits[i] & abbrv_B == args$traits[j]) | (abbrv_A == args$traits[j] & abbrv_B == args$traits[i]), pval]
    }

    gps <- ifelse(length(gps) > 0, gps, NA)

    gps_mat[j, i] <- -log10(gps)
    gps_mat[i, j] <- -log10(gps)
  }
}

pheatmap(gps_mat, breaks = seq(0, 10, length.out = 101), cluster_rows = F, cluster_cols = F, fontsize = 13, filename = args$output, display_numbers = F, number_format = "%.1f", labels_row = rownames(gps_mat), labels_col = colnames(gps_mat))
