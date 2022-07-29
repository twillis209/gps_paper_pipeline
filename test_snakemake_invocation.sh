#!/bin/bash

snakemake -j 30 --local-cores 2 --cores 10 --default-resources runtime=5 mem_mb=3420 tmpdir='tmp' --group-components simulate=100 --use-conda --scheduler greedy --profile "$HOME/.config/snakemake/slurm" --rerun-incomplete --retries 3 results/simgwas/simulated_sum_stats/whole_genome_sum_stats/400_reps/10000_10000/null_sum_stats.tsv.gz
