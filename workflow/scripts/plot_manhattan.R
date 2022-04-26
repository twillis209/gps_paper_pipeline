library(data.table)
library(argparse)
library(pidProjCode)

parser <- ArgumentParser(description = 'Prune merged GWAS summary statistics file')
parser$add_argument('-is', '--sum_stats_file', type = 'character', help = 'Path to summary statistics file')
parser$add_argument('-ib', '--bim_file', type = 'character', help = 'Path to bim file')
parser$add_argument('-ir', '--range_file', type = 'character', help = 'Path to range file')
parser$add_argument('-o', '--output_path', type = 'character', help = 'Path to pruned summary statistics file', required = T)
parser$add_argument('-nt', '--no_of_threads', type = 'integer', help = 'Number of threads to use', default = 1)

args <- parser$parse_args()

setDTthreads(args$no_of_threads)

sum_stats_dat <- fread(args$sum_stats_file, sep = '\t', header = T)

asthma_dat[, c('CHR19', 'BP19', 'a0', 'a1') := tstrsplit(variant, split = ':')]
thyroid_dat[, c('CHR19', 'BP19', 'a0', 'a1') := tstrsplit(variant, split = ':')]

asthma_dat[, CHR19 := paste0('chr', CHR19)]
thyroid_dat[, CHR19 := paste0('chr', CHR19)]

asthma_dat <- asthma_dat[pval > 1e-15]
thyroid_dat <- thyroid_dat[pval > 1e-15]

setnames(asthma_dat, 'pval', 'p')
setnames(thyroid_dat, 'pval', 'p')

asthma_granges <- makeGRangesFromDataFrame(data.frame(asthma_dat), start.field = 'BP19', end.field = 'BP19', seqnames.field = 'CHR19', ignore.strand = T, keep.extra.columns = T)
thyroid_granges <- makeGRangesFromDataFrame(data.frame(thyroid_dat), start.field = 'BP19', end.field = 'BP19', seqnames.field = 'CHR19', ignore.strand = T, keep.extra.columns = T)

pp <- getDefaultPlotParams(plot.type = 4)

axis_label_margin <- -0.05

emptyGRanges <- GRanges(P = numeric())

pp <- getDefaultPlotParams(plot.type=4)

pp$leftmargin <- 0.08
pp$rightmargin <- 0.08

r0 <- 0
cex.main <- 1.3
cex.tick <- 0.75
ymax <- 15

width <- 35
height <- 15

png('ukbb_asthma_manhattan.png', width = width, height = height, unit = 'cm', res = 600)
multitrack_manhattan(gRanges = list(granges),
                             axis_labels = 
                             

kp <- plotKaryotype(plot.type = 4, labels.plotter = NULL, plot.params = pp, ideogram.plotter = NULL)
kp <- kpPlotManhattan(kp, data = asthma_granges, points.col = 'brewer.set3', ymax = ymax, r0 = r0)
dev.off()