library(data.table)
library(argparse)
library(ggplot2)

parser <- ArgumentParser(description = 'Plot heatmap of GPS statistic denominator')
parser$add_argument('-i', '--input_file', type = 'character', help = 'Path to input file')
parser$add_argument('-o', '--output_file', type = 'character', help = 'Path to output file')

args <- c('--input_file', 'results/gps/ukbb/all_pruned_snps/window_1000kb_step_50/20002_1111-20002_1113_ecdf.tsv', '-o',"results/plots/all_pruned_snps/gps_heatmaps/20002_1111-20002_1113.png")

args <- parser$parse_args()

dat <- fread(args$input_file, sep = '\t', header = T)

dat[ , denom := sqrt(F_u*F_v - (F_u^2)*(F_v^2))]
dat[, num := abs(F_uv - F_u*F_v)]
dat[, maximand := sqrt(.N/log(.N))*num/denom]
