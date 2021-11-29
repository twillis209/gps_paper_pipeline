library(data.table)
library(argparse)
library(pheatmap)

parser <- ArgumentParser(description = 'Plot rg vs. Hoeffding\'s p-values heatmap')
parser$add_argument('-f', '--hoeffdings_file', type = 'character', help = 'Path to Hoeffding\'s results file')
parser$add_argument('-r', '--rg_file', type = 'character', help = 'rg file path')
parser$add_argument('-t', '--traits', type = 'character', help = 'Comma-delimited list of traits in desired order')
parser$add_argument('-p', '--rg_pvalue_output', type = 'character', help = 'Path to output file', required = T)
parser$add_argument('-e', '--rg_estimate_output', type = 'character', help = 'Path to output file', required = T)

args <- parser$parse_args()

args$traits <- unlist(strsplit(args$traits, ','))

hoeffdings_dat <- fread(args$hoeffdings_file, sep = '\t', header = T)

rg_dat <- fread(args$rg_file, sep = '\t', header = T)

hoeffdings_rg_pval_mat <- matrix(nrow = length(args$traits), ncol = length(args$traits))
colnames(hoeffdings_rg_pval_mat) <- args$traits
rownames(hoeffdings_rg_pval_mat) <- args$traits

hoeffdings_rg_estimate_mat <- matrix(nrow = length(args$traits), ncol = length(args$traits))
colnames(hoeffdings_rg_estimate_mat) <- args$traits
rownames(hoeffdings_rg_estimate_mat) <- args$traits

for(i in 1:length(args$traits)) {
  for(j in i:length(args$traits)) {
    if(i == j) {
      hoeffdings_pval <- NA
      rg_pval <- NA
      rg_estimate <- NA
    } else {
      hoeffdings_pval <- hoeffdings_dat[(abbrv_A == args$traits[i] & abbrv_B == args$traits[j]) | (abbrv_A == args$traits[j] & abbrv_B == args$traits[i]), p.value]
      rg_pval <- rg_dat[(trait_A == args$traits[i] & trait_B == args$traits[j]) | (trait_A == args$traits[j] & trait_B == args$traits[i]), pval]
      rg_estimate <- rg_dat[(trait_A == args$traits[i] & trait_B == args$traits[j]) | (trait_A == args$traits[j] & trait_B == args$traits[i]), rg]
    }

    hoeffdings_pval <- ifelse(length(hoeffdings_pval) > 0, hoeffdings_pval, NA)
    rg_pval <- ifelse(length(rg_pval) > 0, rg_pval, NA)
    rg_estimate <- ifelse(length(rg_estimate) > 0, rg_estimate, NA)

    hoeffdings_rg_pval_mat[i, j] <- -log10(hoeffdings_pval)
    hoeffdings_rg_pval_mat[j, i] <- -log10(rg_pval)
    hoeffdings_rg_estimate_mat[i, j] <- -log10(hoeffdings_pval)
    hoeffdings_rg_estimate_mat[j, i] <- 10*rg_estimate
  }
}

pheatmap(hoeffdings_rg_pval_mat, breaks = seq(0, 10, length.out = 101), cluster_rows = F, cluster_cols = F, fontsize = 13, filename = args$rg_pvalue_output, display_numbers = F, number_format = "%.1f", labels_row = rownames(hoeffdings_rg_pval_mat), labels_col = colnames(hoeffdings_rg_pval_mat))

pheatmap(hoeffdings_rg_estimate_mat, breaks = seq(0, 10, length.out = 101), cluster_rows = F, cluster_cols = F, fontsize = 13, filename = args$rg_estimate_output, display_numbers = F, number_format = "%.1f", labels_row = rownames(hoeffdings_rg_estimate_mat), labels_col = colnames(hoeffdings_rg_estimate_mat))
