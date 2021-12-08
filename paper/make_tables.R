library(data.table)
library(xtable)

ukbb_traits_cases_controls <- function() {

  traits_dat <- fread('resources/ukbb_sum_stats/traits_codes_abbrv_cases.tsv', sep = '\t', header = T)

  capitalise <- function(x) {
    paste0(toupper(substr(x, 1,1)), substr(x, 2, length(x)))
  }

  traits_dat[, Trait := capitalise(abbrv)]

  traits_dat[, Trait := plyr::mapvalues(from = c('Eczema/derm', 'Emphysema/chronic bronc', 'IHD', 'CD', 'UC', 'IBS', 'MD'), to = c('Eczema/dermatitis', 'Emphysema/chronic bronchitis', 'Ischaemic heart disease', 'Crohn\'s disease', 'Ulcerative colitis', 'Irritable bowel syndrome', 'Macular degeneration'), Trait)]

  traits_dat <- traits_dat[!(Trait %in% c('Diabetes', 'Pid'))]

  traits_dat[, n_cases := format(n_cases, big.mark = ',')]
  traits_dat[, n_controls := format(n_controls, big.mark = ',')]
  traits_dat <- traits_dat[order(n_cases)]

  setnames(traits_dat, c('n_cases', 'n_controls'), c('No. of cases', 'No. of controls'))

  x <- xtable(traits_dat[,c('Trait', 'No. of cases', 'No. of controls')])

  align(x) <- "ll|r|r"

  caption(x) <- "Selected UK Biobank traits used to study the properties of the GPS test in this work."

  label(x) <- "tab:ukbb_cases_controls"

  print(x, include.rownames = F, hline.after = c(0))
}
