import sys

exec(open('workflow/rules/simgwas/simgwas_export_functions.py', 'r').read())

effect_blocks = ['m25', 'm50', 's400', 's200-m25']

if not sys.argv[1]:
    raise Exception("Must specify effect block as first argument")
elif sys.argv[1] not in effect_blocks:
    raise Exception(f"Must specify effect block in {', '.join(effect_blocks)}")
else:
    effect_block = sys.argv[1]

subset = f"a_blocks == \'{effect_block}\'"

if len(sys.argv) > 2:
    tests = sys.argv[2]
else:
    tests = 'ldsc,sumher,gps,hoeffdings'

if 'ldsc' in tests:
    ldsc_files = get_existing_test_files("results/simgwas/simulation_parameters.tsv", reps = 400, filetype = 'ldsc', subset = subset)

    ldsc_out = f"results/ldsc/simgwas/400_reps/randomised/existing_compiled_{effect_block}_results.tsv"

    ldsc_daf = compile_ldsc_results_into_daf(ldsc_files)
    ldsc_daf.to_csv(ldsc_out, sep = '\t', index = False)

    print('ldsc compiled')

if 'sumher' in tests:
    sumher_files = get_existing_test_files("results/simgwas/simulation_parameters.tsv", reps = 400, filetype = 'sumher', subset = subset)

    sumher_out = f"results/ldak/ldak-thin/simgwas/400_reps/randomised/rg/existing_compiled_{effect_block}_results.tsv"

    sumher_daf = compile_sumher_results_into_daf(sumher_files)
    sumher_daf.to_csv(sumher_out, sep = '\t', index = False)

    print('sumher compiled')

if 'hoeffdings' in tests:
    hoeffdings_files = get_existing_test_files("results/simgwas/simulation_parameters.tsv", reps = 400, filetype = 'hoeffdings', subset = subset)

    hoeffdings_out = f"results/hoeffdings/simgwas/400_reps/randomised/existing_compiled_{effect_block}_results.tsv"

    hoeffdings_daf = compile_hoeffdings_results_into_daf(hoeffdings_files)
    hoeffdings_daf.to_csv(hoeffdings_out, sep = '\t', index = False)

    print('hoeffdings compiled')

if 'gps' in tests:
    gps_files = get_existing_test_files("results/simgwas/simulation_parameters.tsv", reps = 400, filetype = 'gps', subset = subset)

    gps_out = f"results/gps/simgwas/400_reps/randomised/existing_compiled_{effect_block}_results.tsv"

    gps_daf = compile_gps_results_into_daf(gps_files)
    gps_daf.to_csv(gps_out, sep = '\t', index = False)

    print('gps compiled')

#    theo_files = [y for y in get_theo_rg_files('results/simgwas/simulation_parameters.tsv', reps = 400, subset = subset) if os.path.exists(y)]
#
#    theo_out = f"results/ldsc/simgwas/400_reps/randomised/existing_compiled_{x}_theo_rg.tsv"
#
#    theo_daf = compile_theo_rg_results_into_daf(theo_files)
#    theo_daf.to_csv(theo_out, sep = '\t', index = False)
