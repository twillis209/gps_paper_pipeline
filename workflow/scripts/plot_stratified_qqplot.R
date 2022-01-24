library(argparse)
library(ggplot2)
library(data.table)
setDTthreads(4)

parser <- ArgumentParser(description = 'Plots stratified Q-Q plots for pair of traits')
parser$add_argument('-i', '--input_file', type = 'character', help = 'Path to file containing UKBB summary statistics', required = T)
args <- parser$parse_args()
parser$add_argument('-a', '--trait_a', type = 'character', help = 'Trait A', required = T)
parser$add_argument('-b', '--trait_b', type = 'character', help = 'Trait B', required = T)
parser$add_argument('-o', '--output_file', type = 'character', help = 'Path to output file', required = T)

args <- parser$parse_args()

stratified_qqplot <- function(data_frame, prin_value_label, cond_value_label = NULL, thresholds = c(1, 1e-1, 1e-2, 1e-3, 1e-4)) {

  data_frame$negLogP <- -log10(data_frame[, prin_value_label])

  if(is.null(cond_value_label)) {
    daf <- data_frame[, c(prin_value_label, 'negLogP')]
    daf <- daf[order(daf[,prin_value_label]), ]
    daf$pp <- -log10(ppoints(nrow(daf)))
    daf$threshold <- factor(c(1))
  } else {
    data_frame <- data_frame[, c(prin_value_label, cond_value_label, 'negLogP')]

    dafs <- list()

    for(i in seq_along(thresholds)) {
      daf <- subset(data_frame, get(cond_value_label) < thresholds[i])
      daf <- daf[order(daf[ , prin_value_label]) , ]
      daf$pp <- -log10(ppoints(nrow(daf)))
      daf$threshold <- factor(thresholds)[i]
      dafs[[i]] <- daf
    }

    daf <- do.call(rbind, dafs)
  }

  ggplot(data=daf) + geom_line(aes(x = pp, y = negLogP, group = threshold, colour = threshold)) + geom_abline(intercept=0,slope=1, linetype="dashed")
}

sum_stats_dat <- fread(args$input_file, sep = '\t', header = T, select = c('variant', args$trait_a, args$trait_b))

a_b_plot <- stratified_qqplot(data.frame(sum_stats_dat), prin_value_label = args$trait_a, cond_value_label = args$trait_b)
b_a_plot <- stratified_qqplot(data.frame(sum_stats_dat), prin_value_label = args$trait_b, cond_value_label = args$trait_a)
