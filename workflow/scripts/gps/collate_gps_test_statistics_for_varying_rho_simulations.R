library(data.table)
library(stringr)

snakemake@input[['gev']]
snakemake@input[['li']]

# returns atomic vector, we want indices 2:6
str_match(x, "results/rho_simulations/(?<rho>0_\\d+)/(?<zmean>\\d)_(?<zsd>\\d)/(?<rep>\\d+)/3000_draws/gps_pvalue.tsv")

str_match(x, "results/rho_simulations/(?<rho>0_\d+)/(?<zmean>\d)_(?<zsd>\d)/(?<rep>\d+)/li_gps_pvalue.tsv")
