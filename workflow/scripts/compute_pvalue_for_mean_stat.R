library(sn)
library(data.table)

stat <- fread(snakemake@input[['mean_stat']], sep = '\t')
perms <- scan(snakemake@input[['permutations']], skip = 1, what = double())

sn_fit <- selm(-log10(y) ~ 1, data = data.frame(y = perms))

pval <- psn(-log10(stat$GPS), xi = sn_fit@param$dp[1], omega = sn_fit@param$dp[2], alpha = sn_fit@param$dp[3])

fwrite(data.table(Trait_A = stat$Trait_A, Trait_B = stat$Trait_B, xi = sn_fit@param$dp[1], omega = sn_fit@param$dp[2], alpha = sn_fit@param$dp[3], GPS = stat$GPS, pval = pval), sep = '\t', file = snakemake@output[[1]])
