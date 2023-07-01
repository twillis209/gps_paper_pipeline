library(data.table)
library(mvtnorm)

##' simulate correlated p-values
##'
##' @param N00 number of snps null for both datasets
##' @param N01 number of snps null for dataset 1, non-null for dataset 2
##' @param N10 number of snps null for dataset 2, non-null for dataset 1
##' @param rho correlation between z scores
##' @param zmean mean z score at non-null snps
##' @return (N00+N01+N10) x 2 matrix of z scores
##' @export 
##' @author Chris Wallace
simp=function(N00, N01, N10, rho, zmean=ZM, zsd=ZS) {
  S=matrix(c(1, rho, rho, 1), 2, 2)

  Z=rmvnorm(N00+N01+N10, sigma=S)

  if(N01>0) {
    Z[ N00+(1:N01), 2]= Z[ N00+(1:N01), 2] + rnorm(N01,zmean,sd=zsd)
  }

  if(N10>0) {
    Z[ N00+N01+(1:N10), 1]= Z[ N00+N01+(1:N10), 1] + rnorm(N10,zmean,sd=zsd)
  }

  P=2 * pnorm(-abs(Z))

  P
}

fwrite(data.table(simp(N00 = snakemake@params[['N00']], N01 = snakemake@params[['N01']], N10 = snakemake@params[['N10']], rho = snakemake@params[['rho']], zmean = snakemake@params[['zmean']], zsd = snakemake@params[['zsd']])), sep = '\t', file = snakemake@output[[1]])
