library(data.table)
library(ggplot2)
library(patchwork)
library(Hmisc)
library(scales)
library(mvtnorm)

theme_set(
  theme_bw()+
  theme(
                                        #text = element_text(family = 'LM Roman 10'),
    axis.title = element_text(size = 12),
    plot.title = element_text(hjust = 0.5, size = 18),
    strip.text = element_text(size = 10),
    axis.text.x = element_text(size = 8, angle = 90, color = "black"),
    axis.text.y = element_text(size = 8, color = "black"),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10),
    plot.tag = element_text(size = 10),
    strip.background = element_rect(fill = 'white')
  )
)

stat_palette <- hue_pal()(5)
gps_col <- stat_palette[1]
li_gps_col <- stat_palette[5]

test_stats <- fread(snakemake@input[['test_statistics']], sep = '\t', header = T)

t1e_dat <- test_stats[, .(true = sum(pvalue <= 0.05), false = sum(pvalue > 0.05)), by = .(rho, zmean, zsd, stat)]

t1e_dat <- t1e_dat[, .(rho, zmean, zsd, stat, true, n = true+false, binconf(x = true, n = true+false, alpha = 0.05, method = 'wilson', return.df = T))]

t1e_dat[, propsig := pnorm(qnorm(5e-8/2), mean = -zmean, sd = sqrt(1+zsd^2))]
t1e_dat[, nsig := pnorm(qnorm(5e-8/2), mean = -zmean, sd = sqrt(1+zsd^2))*400]

#pval_dats <- lapply(snakemake@input[['pvalue_files']], fread, sep = '\t', header = T, col.names = c('p.1', 'p.2'))
#
#zmeans <- c(1, 2, 3, 1, 2, 3)
#zsds <- c(1, 1, 1, 2, 2, 2)
#
#for(i in 1:6) {
#  pval_dats[[i]][, `:=` (zmean = zmeans[i], zsd = zsds[i])]
#  pval_dats[[i]][1:40000, `:=` (null.1 = F, null.2 = F)]
#  pval_dats[[i]][40001:40400, `:=` (null.1 = T, null.2 = F)]
#  pval_dats[[i]][40401:40800, `:=` (null.1 = F, null.2 = T)]
#}
#
#pval_dat <- rbindlist(pval_dats)

simp=function(N00, N01, N10, rho, zmean=ZM, zsd=ZS) {
  S=matrix(c(1, rho, rho, 1), 2, 2)
  Z=rmvnorm(N00+N01+N10, sigma=S)
  if(N01>0) 
    Z[ N00+(1:N01), 2]= Z[ N00+(1:N01), 2] + rnorm(N01,zmean,sd=zsd)
  if(N10>0) 
    Z[ N00+N01+(1:N10), 1]= Z[ N00+N01+(1:N10), 1] + rnorm(N10,zmean,sd=zsd)
  P=2 * pnorm(-abs(Z))
  if(which.max(P[,1])==which.max(P[,2])) { # catch rare complication
    return(simz(N00, N01, N10, rho=r,zmean=5))
  }
  P 
}

p11=simp(40000, 400, 400, rho=0, zmean=1, zsd=1)[,2]

dt=data.table(p=p11, lp=-log10(p11), assoc=rep(c("null","assoc"), times=c(40400, 400)))

lp=function(zmean, zsd, w=400/41200) {
  x=seq(-10,0,by=0.1)
  y0=2*dnorm(x, mean=0, sd=1) #* (1-w)
  y1=2*dnorm(x, mean=-zmean, sd=sqrt(zsd^2 + 1)) #* w
  lp=-log10(pnorm(x)*2)
  data.table(z=x, prob0=y0, prob1=y1, prob=y0+y1, lp=lp, zmean=zmean, zsd=zsd)
}

dt <- melt(rbind(lp(1,1), lp(1,2), lp(2,1), lp(2,2), lp(3,1), lp(3,2)), measure.vars=c("prob0","prob1"), variable.name="associated")

dt[,associated:=c("null","associated")[associated]]

pl1 <- ggplot(dt[zmean > 0]) +
  geom_path(aes(x=lp,y=value,col=associated)) +
  facet_grid(zmean ~ zsd, labeller=label_both) +
  labs(x="-log10 p-value", y="Density") +
  scale_colour_manual(name = '', values=c(null="grey20",associated="slateblue")) +
  geom_vline(xintercept=-log10(5e-8), linetype="dashed", colour="black") +
  geom_text(aes(label=paste0(100*round(propsig,3),"%")),x=15,y=0.4, data = t1e_dat[zmean > 0], col="slateblue") +
  theme(legend.position="bottom")

pl2 <- ggplot(data = t1e_dat[zmean > 0], aes(x = rho, y = PointEst, ymin = Lower, ymax = Upper, col = stat, group = stat)) +
  geom_pointrange(size = 0.05)+
  geom_path()+
  geom_hline(yintercept = 0.05, linetype = 2) +
  scale_colour_discrete("Method")+
  labs(x="between dataset correlation, rho", y="estimated type 1 error rate") +
  facet_grid(zmean ~ zsd, labeller=label_both) +
  ylim(c(0, .375))+
  scale_colour_manual(name = '', values = c('GPS-Exp' = li_gps_col, 'GPS-GEV' = gps_col))+
  theme(legend.position="bottom")+
  ylab('Type 1 error rate')+
  xlab('Between-data set effect estimate correlation')

fig_s10 <- pl1 / pl2+
  plot_annotation(tag_levels = 'A')

ggsave(fig_s10, file = snakemake@output[[1]], width = 8, height = 10)
