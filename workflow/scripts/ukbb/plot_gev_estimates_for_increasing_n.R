library(ggplot2)
library(ggpubr)
library(argparse)

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

parser <- ArgumentParser(description = 'Plot GEV parameter estimates as a function of sample size')
parser$add_argument('-f', '--fit_file', type = 'character', help = 'Path to fitted parameters file')
parser$add_argument('-o', '--output_file', type = 'character', help = 'Path to output file')

args <- parser$parse_args()

daf <- read.table(args$fit_file, sep = '\t', header = T)

pl_loc_estimate <- ggplot(daf[c('n', 'loc', 'loc.sd')], aes(x=n, y=loc))+
  geom_line()+
  geom_errorbar(aes(ymin=loc-1.96*loc.sd, ymax=loc+1.96*loc.sd))+
  ggtitle('Location parameter')+
  ylab('Estimate')+
  theme(legend.position = 'none')

pl_scale_estimate <- ggplot(daf[c('n', 'scale', 'scale.sd')], aes(x=n, y=scale))+
  geom_line()+
  geom_errorbar(aes(ymin=scale-1.96*scale.sd, ymax=scale+1.96*scale.sd))+
  ggtitle('Scale parameter')+
  ylab('Estimate')+
  theme(legend.position = 'none')

pl_shape_estimate <- ggplot(daf[c('n', 'shape', 'shape.sd')], aes(x=n, y=shape))+
  geom_line()+
  geom_errorbar(aes(ymin=shape-1.96*shape.sd, ymax=shape+1.96*shape.sd))+
  ggtitle('Shape parameter')+
  ylab('Estimate')+
  theme(legend.position = 'none')

ggsave(plot = ggarrange(plotlist = list(pl_loc_estimate, pl_scale_estimate, pl_shape_estimate), nrow = 1, ncol = 3, common.legend = T), file = args$output_file, width = 8.1, height = 3)
