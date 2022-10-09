# Pipeline for the manuscript 'Accurate detection of shared genetic architecture from GWAS summary statistics in the small-sample context'

This repository contains a `snakemake` pipeline for the reproduction of the results in the manuscript named above.

## Running the pipeline to obtain the manuscript results

### Running real data analyses

The file `workflow/rules/export.smk` contains rules which specify the main real data outputs of the pipeline used in the manuscript.

* `ukbb_with_mhc`: outputs the GPS, Hoeffding's, LDSC, and SumHer test statistics for the MHC-inclusive analyses of real UKBB GWAS data sets
* `ukbb_sans_mhc`: outputs the GPS, Hoeffding's, LDSC, and SumHer test statistics for the sans-MHC analyses of real UKBB GWAS data sets
* `asthma_emphysema_gof`: outputs the data sets used to generate Figure 1 and Supplementary Figures 1, 2, and 3

### Running `simGWAS` simulations and analyses

Running the simulations and their downstream analyses is extremely compute-intensive. At the time of writing, construction of the task graph by `snakemake` is single-threaded in a single instance and generating the task graph for only one of the five simulation regimes (actually four in the pipeline, see note below) would itself take >10 hours. We therefore recommend using `snakemake`'s `batch` feature to run each regime in tenths to cut down on the task graph generation time.

For example, to invoke `snakemake` to run the first tenth of the small-effect regime analyses:

```
snakemake run_s400_simulations_and_analyses --batch run_s400_simulations_and_analyses=1/10
```

Batches one to nine should be non-overlapping, so multiple instances of `snakemake` can be run in parallel by specifying the `--nolock` option in the invocation. This effectively parallelises task graph generation.

#### Notes on the five simulation regimes

We specify the simulation regimes in the pipeline by labels which give their effect size and the number of causal variants. 

* Small-effect regime: `s400` (`400` `s`mall-effect variants)
* Large-effect I regime: `m25` (`25` `m`edium-effect variants)
* Large-effect II regime: `m50` (`50` `m`edium-effect variants)
* Mixed-effect regimes I and II: `s200-m25` (`200` `s`mall-effect and `25` `m`edium-effect variants)

'Large' is aliased to 'medium' here as our large-effect regimes were formerly medium-effect regimes.

We split the two mixed-effect regimes from a single set, hence `s200-m25` accounting for both here.

Result for each regime can be run with one of the following rules in `workflow/rules/simgwas/simgwas_export.smk`:

* `run_m25_simulations_and_analyses`
* `run_m50_simulations_and_analyses`
* `run_s400_simulations_and_analyses`
* `run_s200-m25_simulations_and_analyses`

#### Disk usage

We generate simulated GWAS data sets by combining LD block-sized simulation results. Whilst we keep the latter, the former are marked as temporary with `temp` and `snakemake` should delete them once they're no longer required for execution. Each replicate of two GWAS data sets occupies about ~2G disk space at most during execution; how long this is occupied with depend on how many threads you allocate to each job (see below).

If you are making use of `snakemake`'s `group` feature to bundle together related jobs for cluster submission, you should take care that the files marked `temp` are actually deleted during execution; `temp` and `group`, together, seem not to be an entirely mature feature combination and `snakemake` will not always delete temporary files when it should. The upshot of this is that you should keep an eye on your disk usage on shared drives whilst the pipeline is running. It may be necessary to delete the files from `results/simgwas/simulated_sum_stats/whole_genome_sum_stats` and `results/simgwas/simulated_sum_stats/munged_sum_stats` if `snakemake` is not able.

#### Resource allocation

If running `snakemake` in cluster execution mode (it's hard to imagine you would be doing otherwise with this pipeline), you can set the resources for each rule in your profile file. The resource specifications given in the `smk` are to be taken as recommendations rather than hard and fast specifications; what is suitable for execution on your own cluster will depend on the compute resources available to you. 

We ran jobs in bundles and used the relatively new `group-resources` feature when running the pipeline. An excerpt from our profile `config.yaml` file is below:

```
local-cores: 8
cores: 12
keep-going: True
use-conda: True
scheduler: greedy
default-resources:
  - runtime=5
  - mem_mb=3420
  - tmpdir='tmp'
resources:
  - concurrent_sans_permute_jobs=200
set-threads:
  - merge_randomised_simulated_sum_stats=12
  - prune_merged_randomised_simulated_sum_stats=12
  - process_combined_simgwas_sum_stats=12
  - calculate_theoretical_rg=1
  - permute_sim_pair=12
set-resources:
  - merge_randomised_simulated_sum_stats:
      - runtime=25
  - prune_merged_randomised_simulated_sum_stats:
      - runtime=25
  - process_combined_simgwas_sum_stats:
      - runtime=25
group-components:
  - simulate=24
  - permutation=1
  - ldsc_hoeffding_sumher_gps_sans_permutation=1
  - calculate_theoretical_rg=1000
```

Each `group` submission was run with 12 cores, hence our use of 12 threads for most of the jobs belonging to the `ldsc_hoeffding_sumher_gps_sans_permutation` group, which carries out almost all steps of a single replicate's analysis but for the GPS permutation testing and GEV fitting steps. 

Note the user-defined `concurrent_sans_permute_jobs` resource. We used this to limit the number of concurrent `ldsc_hoeffding_sumher_gps_sans_permutation` group submissions to 200. Note that this is *per `snakemake` instance*; if you are running the pipeline in batches with concurrent `snakemake` instances, each of these will run a maximum of 200 such jobs concurrently, so take care not to fill up your hard drive!

## Dependencies

### `conda` environments

`environment.yaml` specifies the `conda` environment used to run the python scripts in the pipeline with the exception of the LDSC script, which has its own `conda` environment at `workflow/rules/ldsc/envs/ldsc.ymal`. The first environment should be activated prior to launching the pipeline whilst `snakemake` itself should manage the one-timeinstallation and activation of the LDSC environment through the `conda` directive in the relevant rules.

Note that `snakemake` is not contained within the main `conda` environment as the version of `snakemake` with which the pipeline was developed, `7.14.0`, was not available through `conda` at the time of development. The pipeline should be compatible with versions of `snakemake` `>=7.14.0`.

At the moment, we are not using `conda` (or any other environment/package manager) to manage the R package dependencies of the scripts in the pipeline.

### LDSC

LDSC should be installed as directed on the [software's repository](https://github.com/bulik/ldsc). To use the pipeline as written, create an environment variable `ldsc` which points to the directory in which the `ldsc.py` script is located. The `conda` environment required to run LDSC is located at `workflow/rules/ldsc/envs/ldsc.ymal` and should be managed by `snakemake` during execution.

### GPS

The software to compute the GPS statistic and p-values for the test is a mixture of R and C++. The C++ component is contained in the submodule `workflow/scripts/gps_cpp`. Its programs are built with `cmake` and `make`. To build the software `cd` into `gps_cpp` and run:

```
mkdir build
cd build
cmake ..
make
```

This should build the C++ programs the pipeline depends on.

### SumHer

Install LDAK (SumHer's parent software suite) as per the directions on the [software's website](http://dougspeed.com/downloads2/). To use the pipeline as written, create an environment variable `ldakRoot` which points to the directory in which the LDAK binary is located.
