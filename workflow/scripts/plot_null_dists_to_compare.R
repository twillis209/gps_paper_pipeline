library(argparse)
library(data.table)
library(ggplot2)

theme_set(theme_bw()+
          theme(
            axis.title = element_text(size=12),
            plot.title=element_text(hjust=0.5, size=12),
            strip.text=element_text(size=10),
            axis.text.x=element_text(size=10, angle=30, color="black"),
            axis.text.y=element_text(size=10, color="black"),
            legend.title=element_text(size=10),
            legend.text=element_text(size=10)
          )
          )

library(ggpubr)

parser <- ArgumentParser(description = 'Plot GPS null distribution using specified null permutations file')
parser$add_argument('-f', '--fitdist_file', type = 'character', help = 'Path to fitted parameters file')
parser$add_argument('-p', '--perm_file', type = 'character', help = 'Path to file contained permuted null GPS statistics')
parser$add_argument('-a', '--exp1_null', type = 'character', help = 'Path to exp1 output plot file', required = T)
parser$add_argument('-b', '--gev_null', type = 'character', help = 'Path to gev output plot file', required = T)
parser$add_argument('-c', '--exp1_gev_combined', type = 'character', help = 'Path to combined output plot file', required = T)
parser$add_argument('-nt', '--no_of_threads', type = 'integer', help = 'Number of threads to use', default = 1)

args <- parser$parse_args()

perm_sample <- scan(args$perm_file, skip = 1)

exp1_pvals <- pexp(perm_sample^-2)

fit_dat <- fread(args$fitdist_file, sep = '\t', header = T)

gev_pvals <- evd::pgev(perm_sample, loc = fit_dat$loc, scale = fit_dat$scale, shape = fit_dat$shape, lower.tail = F)

pl_exp1_pvals_hist <- ggplot(data = data.frame(p = exp1_pvals))+
  geom_histogram(aes(x = p, y = ..count../sum(..count..)), colour = 'black', fill = 'gray', breaks = seq(0, 1, length.out = 21))+
  geom_hline(yintercept = 0.05, linetype = "dashed", col = "blue")+
  xlab('GPS p-value')+
  ylab('Relative frequency')+
  ggtitle('Histogram of GPS p-values\nfrom Exp(1) under null')+
  scale_x_continuous(limits = c(0,1))+
  ylim(0,0.1)

ggsave(plot = pl_exp1_pvals_hist, file = args$exp1_null, units = "in", width = 2.7, height = 3)

pl_gev_pvals_hist <- ggplot(data = data.frame(p = gev_pvals))+
  geom_histogram(aes(x = p, y = ..count../sum(..count..)), colour = 'black', fill = 'gray', breaks = seq(0, 1, length.out = 21))+
  geom_hline(yintercept = 0.05, linetype = "dashed", col = "blue")+
  xlab('GPS p-value')+
  ylab('Relative frequency')+
  ggtitle('Histogram of GPS p-values\nfrom GEV under null')+
  scale_x_continuous(limits = c(0,1))+
  ylim(0,0.1)

ggsave(plot = pl_gev_pvals_hist, file = args$gev_null, units = "in", width = 2.7, height = 3)

ggsave(plot = ggarrange(plotlist = list(pl_exp1_pvals_hist, pl_gev_pvals_hist), ncol = 2, nrow = 1), file = args$exp1_gev_combined, units = "in", width = 5.4, height = 3)
