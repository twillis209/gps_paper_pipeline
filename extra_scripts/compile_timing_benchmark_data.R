library(data.table)


naive_dat <- fread('benchmarks/results/gps/simgwas/400_reps/randomised/10000_10000_10000_10000/s400_s400_s0/window_1000kb_step_50/pert_0/naive/compiled.txt')
lw_dat <- fread('benchmarks/results/gps/simgwas/400_reps/randomised/10000_10000_10000_10000/s400_s400_s0/window_1000kb_step_50/pert_100/lw/compiled.txt')
pp_dat <- fread('benchmarks/results/gps/simgwas/400_reps/randomised/10000_10000_10000_10000/s400_s400_s0/window_1000kb_step_50/pert_1/pp/compiled.txt')

summary(naive_dat)
summary(pp_dat)
summary(lw_dat)
