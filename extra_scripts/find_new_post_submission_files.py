import pandas as pd
import os

exec(open('workflow/rules/simgwas/simgwas_export_functions.py').read())

ncases = [500, 1000, 5000, 10000, 100000]
ncontrols = [10000, 10000, 10000, 10000, 100000]
s400_shared_blocks = ['s0', 's100', 's200', 's300', 's400']
m25_shared_blocks = ['m0', 'm5', 'm10', 'm15', 'm20', 'm25']
m50_shared_blocks = ['m0', 'm10', 'm20', 'm30', 'm40', 'm50']
s200_m25_shared_blocks = ['s0-m0', 's100-m0', 's100-m15', 's100-m25', 's200-m0', 's200-m15', 's200-m25']

li_gps_s400_files = get_test_files("results/simgwas/simulation_parameters.tsv", reps = 400, filetype = 'li_gps', subset = f"a_blocks == \'s400\' & shared_blocks in {s400_shared_blocks} & ncases_A in {ncases} & ncontrols_A in {ncontrols}")
existing_li_gps_s400_files = [x for x in li_gps_s400_files if os.path.exists(x)]

li_gps_m25_files = get_test_files("results/simgwas/simulation_parameters.tsv", reps = 400, filetype = 'li_gps', subset = f"a_blocks == \'m25\' & shared_blocks in {m25_shared_blocks} & ncases_A in {ncases} & ncontrols_A in {ncontrols}")
existing_li_gps_m25_files = [x for x in li_gps_m25_files if os.path.exists(x)]

li_gps_m50_files = get_test_files("results/simgwas/simulation_parameters.tsv", reps = 400, filetype = 'li_gps', subset = f"a_blocks == \'m50\' & shared_blocks in {m50_shared_blocks} & ncases_A in {ncases} & ncontrols_A in {ncontrols}")
existing_li_gps_m50_files = [x for x in li_gps_m50_files if os.path.exists(x)]

li_gps_s200_m25_files = get_test_files("results/simgwas/simulation_parameters.tsv", reps = 400, filetype = 'li_gps', subset = f"a_blocks == \'s200-m25\' & shared_blocks in {s200_m25_shared_blocks} & ncases_A in {ncases} & ncontrols_A in {ncontrols}")

existing_li_gps_s200_m25_files = [x for x in li_gps_s200_m25_files if os.path.exists(x)]

sub_ncases = [1000, 5000, 10000]
sub_ncontrols = [10000, 10000, 10000]
sub_s400_shared_blocks = ['s0', 's200', 's400']

mean_gps_stat_s400_files = get_test_files("results/simgwas/simulation_parameters.tsv", reps = 400, filetype = 'mean_stat', subset = f"a_blocks == \'s400\' & shared_blocks in {sub_s400_shared_blocks} & ncases_A in {sub_ncases} & ncontrols_A in {sub_ncontrols}", draws = '10000')
