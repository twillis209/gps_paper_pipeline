import pandas as pd
import os

exec(open('workflow/rules/simgwas/simgwas_one_chrom_export_functions.py', 'r').read())
exec(open('workflow/rules/simgwas/simgwas_export_functions.py', 'r').read())

ldsc_files = get_one_chrom_test_files('results/simgwas/one_chrom_simulation_parameters.tsv', reps = 400, filetype = 'ldsc', subset = 's60')
existing_ldsc_files = [x for x in ldsc_files if os.path.exists(x)]
ldsc_daf = compile_ldsc_results_into_daf(existing_ldsc_files)
ldsc_daf.to_csv('s60_chr1_ldsc.tsv', sep = '\t', index = False)

sumher_files = get_one_chrom_test_files('results/simgwas/one_chrom_simulation_parameters.tsv', reps = 400, filetype = 'sumher', subset = 's60')
existing_sumher_files = [x for x in sumher_files if os.path.exists(x)]
sumher_daf = compile_sumher_results_into_daf(existing_sumher_files)
sumher_daf.to_csv('s60_chr1_sumher.tsv', sep = '\t', index = False)

gps_files = get_one_chrom_test_files('results/simgwas/one_chrom_simulation_parameters.tsv', reps = 400, filetype = 'gps', draws = 3000, subset = 's60')
existing_gps_files = [x for x in gps_files if os.path.exists(x)]
gps_daf = compile_gps_results_into_daf(existing_gps_files)
gps_daf.to_csv('s60_chr1_gps.tsv', sep = '\t', index = False)

li_gps_files = get_one_chrom_test_files('results/simgwas/one_chrom_simulation_parameters.tsv', reps = 400, filetype = 'li_gps', subset = 's60')
existing_li_gps_files = [x for x in li_gps_files if os.path.exists(x)]
li_gps_daf = compile_li_gps_results_into_daf(existing_li_gps_files)
li_gps_daf.to_csv('s60_chr1_li_gps.tsv', sep = '\t', index = False)

hoeffdings_files = get_one_chrom_test_files('results/simgwas/one_chrom_simulation_parameters.tsv', reps = 400, filetype = 'hoeffdings', subset = 's60')
existing_hoeffdings_files = [x for x in hoeffdings_files if os.path.exists(x)]
hoeffdings_daf = compile_hoeffdings_results_into_daf(existing_hoeffdings_files)
hoeffdings_daf.to_csv('s60_chr1_hoeffdings.tsv', sep = '\t', index = False)

theo_files = get_one_chrom_test_files('results/simgwas/one_chrom_simulation_parameters.tsv', reps = 400, filetype = 'theo', subset = 's60')
existing_theo_files = [x for x in theo_files if os.path.exists(x)]
theo_daf = compile_theo_results_into_daf(existing_theo_files)
theo_daf.to_csv('s60_chr1_theo.tsv', sep = '\t', index = False)
