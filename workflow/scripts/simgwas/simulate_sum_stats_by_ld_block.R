library(data.table)
library(simGWAS)
library(argparse)
library(parallel)
library(magrittr)

simulated_z_score_par <- function(exp_z_score, ld_mat, nrep=1, ncores = 1){
    sim_z_score <- mvnfast::rmvn(n = nrep, mu = exp_z_score, sigma = ld_mat, ncores 
= ncores)
    if(nrep==1)
        return(c(sim_z_score))
    sim_z_score
}

parser <- ArgumentParser(description = 'Uses simGWAS to simulate GWAS summary statistics for a single LD block.')
parser$add_argument('--hap_file', type = 'character', help = 'Path to haplotype file')
parser$add_argument('--leg_file', type = 'character', help = 'Path to legend file')
parser$add_argument('--bim_file', type = 'character', help = 'Path to bim file')
parser$add_argument('--ld_mat_file', type = 'character', help = 'Path to LD matrix file')
parser$add_argument('--chr_no', type = 'integer', help = 'Number of chromosome')
parser$add_argument('--causal_variant_ind', type = 'integer', nargs = '*', help = 'Indices of causal variants within LD block.')
parser$add_argument('--effect_size', type = 'character', help = 'Determines odds ratios for causal variants', required = T)
parser$add_argument('--no_controls', type = 'integer', help = 'No. of controls', default = 10000)
parser$add_argument('--no_cases', type = 'integer', help = 'No. of cases', default = 10000)
parser$add_argument('--no_reps', type = 'integer', help = 'No. of replicates', default = 1)
parser$add_argument('--seed', type = 'integer', help = 'Seed for RNG', default = 1)
parser$add_argument('-o', '--output_file', type = 'character', help = 'Path to output file', required = T)
parser$add_argument('-nt', '--no_of_threads', type = 'integer', help = 'Number of threads to use', default = 1)

args <- parser$parse_args()

setDTthreads(args$no_of_threads)

set.seed(args$seed)

leg_dat <- fread(file = args$leg_file, sep = ' ', header = T)
hap_dat <- fread(file = args$hap_file, sep = ' ', header = F)
bim_dat <- fread(file = args$bim_file, sep = '\t', header = F, col.names = c('chr', 'bim.id', 'Cm', 'bp', 'A1', 'A2'))
load(file = args$ld_mat_file)

if(is.null(args$causal_variant_ind)) {
  args$causal_variant_ind <- sample(1:ncol(ld_mat), size = 1)
}

odds_ratios <- list('null' = 1,
                    'tiny' = 1.02,
                    'small' = 1.05,
                    'infinitesimal' = 1.1,
                    'medium' = 1.2,
                    'large' = 1.4,
                    'vlarge' = 2)

# TODO does not work, do not reimplement without fixing
#args$odds_ratios <- sample(odds_ratios[[args$effect_size]], size = length(args$causal_variant_ind), replace = T)

if(args$causal_variant_ind > ncol(ld_mat)) {
  args$causal_variant_ind <- ncol(ld_mat)/2
}

if(args$effect_size %in% names(odds_ratios)) {
  args$odds_ratios <- log(odds_ratios[[args$effect_size]])
} else if(startsWith(args$effect_size, 'random_')) {
  stringr::str_match(args$effect_size, '(\\d+)-(\\d+)')[1,2:3] %>%
    paste(., collapse = '.') %>%
    as.numeric -> h2_M
  cv_maf <- leg_dat[args$causal_variant_ind, EUR]
  args$odds_ratios <- rnorm(length(args$causal_variant_ind), mean = 0, sd = sqrt(h2_M/(2*(cv_maf*(1-cv_maf)))))
} else {
  stop(sprintf("Effect size %s must be one of the following: '%s.", args$effect_size, paste(names(odds_ratios), collapse = ',')))
}

for(j in 1:(ncol(hap_dat)-2)) {
  set(hap_dat, j = j, value = as.numeric(hap_dat[[j]]))
}

hap_mat <- as.matrix(hap_dat[,1:(ncol(hap_dat)-2)])

hap_meta_dat <- hap_dat[, (ncol(hap_dat)-1):ncol(hap_dat)]

names(hap_meta_dat) <- c('rs', 'block')

rm(hap_dat)

freq_dat <- data.table(t(hap_mat)+1)

rm(hap_mat)

colnames(freq_dat) <- hap_meta_dat$rs

# Required by the make_GenoProbList function below (at least)
freq_dat[, Probability := 1/.N]

chosen_snps <- colnames(ld_mat)

cv_ind <- args$causal_variant_ind

cv_snp <- chosen_snps[cv_ind]

sub_freq_dat <- freq_dat[, c(colnames(ld_mat), 'Probability'), with = F]

geno_probs <- make_GenoProbList(snps = colnames(ld_mat),
                                W = colnames(ld_mat)[cv_ind],
                                freq = sub_freq_dat)

zexp <- expected_z_score(N0 = args$no_controls,
                          N1 = args$no_cases,
                          snps = chosen_snps,
                          W = cv_snp,
                          gamma.W = args$odds_ratios,
                          freq = sub_freq_dat,
                          GenoProbList = geno_probs)

zsim <- simulated_z_score_par(exp_z_score = zexp, ld_mat = ld_mat, nrep = args$no_reps, ncores = args$no_of_threads)

# Both vbetasim and betasim overflow if passed as integers
vbetasim <- simulated_vbeta(N0 = as.numeric(args$no_controls),
                            N1 = as.numeric(args$no_cases),
                            snps = chosen_snps,
                            W = cv_snp,
                            gamma.W = args$odds_ratios,
                            freq = sub_freq_dat,
                            nrep = args$no_reps)

betasim <- zsim * sqrt(vbetasim)

if(args$no_reps == 1 ) {
  result_dat <- data.table(leg_dat[rs %in% chosen_snps, .(id, position, block, a0, a1, TYPE, EUR)], zexp, zsim, vbetasim = t(vbetasim), betasim = t(betasim))
} else {
  # Add p-values for multiple reps
  result_dat <- data.table(leg_dat[rs %in% chosen_snps, .(id, position, block, a0, a1, TYPE, EUR)], zexp, t(zsim), t(vbetasim), t(betasim))
}

names(result_dat) <- c('id', 'position', 'block', 'a0', 'a1', 'TYPE', 'EUR', 'zexp',
                       paste0('zsim.', 1:args$no_reps),
                       paste0('vbetasim.', 1:args$no_reps),
                       paste0('betasim.', 1:args$no_reps))

for(j in 1:args$no_reps) {
  result_dat[, c(paste0('p.', j)) := 2*pnorm(abs(get(paste0('zsim.', j))), lower.tail = F)]
}

result_dat[, or := 1]
result_dat[cv_ind, or := exp(args$odds_ratios)]
setnames(result_dat, 'or', 'chosen_or')

result_dat[, `:=` (ncases = args$no_cases, ncontrols = args$no_controls)]

# This merging step confines our SNPs to those present in the bim file; we keep the same order of a0 and a1 alleles, even if this differs in the bim file
result_dat <- merge(result_dat, bim_dat[, .(bim.id, bp, A1, A2)], by.x = 'position', by.y = 'bp', all.x = T)

result_dat <- result_dat[(a0 == A2 & a1 == A1) | (a0 == A1 & a1 == A2)]

result_dat[, c("A1", "A2", "bim.id") := NULL]

result_dat[, chr := args$chr_no]

result_dat[, id := tstrsplit(id, split = ':')[[1]]]

# NB: We add 'block_effect_size' to the end of this
cols <- c("position", "a0", "a1", "id", "block", "TYPE", "EUR", "zexp", paste0("zsim.", 1:args$no_reps), paste0("vbetasim.", 1:args$no_reps), paste0("betasim.", 1:args$no_reps), paste0("p.", 1:args$no_reps), "chosen_or", "ncases", "ncontrols", "chr")

result_dat <- result_dat[, ..cols]

result_dat[, block_effect_size := args$effect_size]

fwrite(result_dat, file = args$output_file, sep = '\t', col.names = F)
