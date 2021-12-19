library(ggplot2)
library(ggpubr)
library(argparse)
require(scales)

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

parser <- ArgumentParser(description = 'Plot GEV parameter estimates as a function of no. of SNPs')
parser$add_argument('-i', '--fit_files', type = 'character', nargs = '+', help = 'Path to fitted parameters files')
parser$add_argument('-n', '--no_snps', type = 'integer', nargs = '+', help = 'No. of SNPs used to fit each file\' estimates')
parser$add_argument('-o', '--output_file', type = 'character', help = 'Path to output file')

args <- parser$parse_args()

if(length(args$fit_files) != length(args$no_snps)) {
  stop("Lengths of fit_files and no_snps do not match")
}

daf <- data.frame(trait_A = character(), trait_B = character(), no_snps = integer(), loc = numeric(), loc.sd = numeric(), scale = numeric(), scale.sd = numeric(), shape = numeric(), shape.sd = numeric())

for(i in seq_along(args$fit_files)) {
  fit_daf <- read.table(args$fit_files[i], sep = '\t', header = T)

  daf <- rbind(daf, data.frame(trait_A = fit_daf$trait_A,
                               trait_B = fit_daf$trait_B,
                               no_snps = args$no_snps[i],
                               loc = fit_daf$loc,
                               loc.sd = fit_daf$loc.sd,
                               scale = fit_daf$scale,
                               scale.sd = fit_daf$scale.sd,
                               shape = fit_daf$shape,
                               shape.sd = fit_daf$shape.sd))
}

pl_loc_estimate <- ggplot(daf[c('no_snps', 'loc', 'loc.sd')], aes(x=no_snps, y=loc))+
  geom_line()+
  geom_errorbar(aes(ymin=loc-1.96*loc.sd, ymax=loc+1.96*loc.sd))+
  ggtitle('Location parameter')+
  ylab('Estimate')+
  theme(legend.position = 'none')+
  scale_x_continuous(labels = scales::comma)

pl_scale_estimate <- ggplot(daf[c('no_snps', 'scale', 'scale.sd')], aes(x=no_snps, y=scale))+
  geom_line()+
  geom_errorbar(aes(ymin=scale-1.96*scale.sd, ymax=scale+1.96*scale.sd))+
  ggtitle('Scale parameter')+
  ylab('Estimate')+
  theme(legend.position = 'none')+
  scale_x_continuous(labels = scales::comma)

pl_shape_estimate <- ggplot(daf[c('no_snps', 'shape', 'shape.sd')], aes(x=no_snps, y=shape))+
  geom_line()+
  geom_errorbar(aes(ymin=shape-1.96*shape.sd, ymax=shape+1.96*shape.sd))+
  ggtitle('Shape parameter')+
  ylab('Estimate')+
  theme(legend.position = 'none')+
  scale_x_continuous(labels = scales::comma)

ggsave(plot = ggarrange(plotlist = list(pl_loc_estimate, pl_scale_estimate, pl_shape_estimate), nrow = 1, ncol = 3, common.legend = T), file = args$output_file, width = 8.1, height = 3)
