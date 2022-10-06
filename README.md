# Pipeline for the manuscript 'Accurate detection of shared genetic architecture from GWAS summary statistics in the small-sample context'

## Running the pipeline to obtain the manuscript results

The file `workflow/rules/export.smk` contains rules which specify the main outputs of the pipeline used in the manuscript.

* `ukbb_with_mhc`: outputs the GPS, Hoeffding's, LDSC, and SumHer test statistics for the MHC-inclusive analyses of real UKBB GWAS data sets
* `ukbb_sans_mhc`: outputs the GPS, Hoeffding's, LDSC, and SumHer test statistics for the sans-MHC analyses of real UKBB GWAS data sets
* `asthma_emphysema_gof`: outputs the data sets used to analyse the dependence of parameter estimates of the generalised extreme value distribution (GEVD) on the number of permutations and SNPs, and to examine the goodness-of-fit of the GEVD to the 
* ``

## Dependencies

### `conda` environments

`environment.yaml` specifies the `conda` environment used to run the python scripts in the pipeline with the exception of the LDSC script, which has its own `conda` environment at `workflow/rules/ldsc/envs/ldsc.ymal`. The first environment should be activated prior to launching the pipeline whilst `snakemake` itself should manage the use of the LDSC environment through the `conda` directive in the relevant rules.

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
