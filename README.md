rrbgen
======
**__Current Version: 0.0.6__**
Release date: Dec 20, 2018

A lightweight limited functionality R bgen read/write library

[![Build Status](https://travis-ci.org/rwdavies/rrbgen.svg)](https://travis-ci.org/rwdavies/rrbgen)

This library supports "v1.3" of bgen defined [here](http://www.well.ox.ac.uk/~gav/bgen_format/). It supports reading and writing using 8, 16, 24 or 32 bits per probability, using Layout = 2 and CompressedSNPBlocks = 1, for bi-allelic SNPs with samples of ploidy 2. Any other format specifications may crash unexpectedly without properly defined error.

### Quick start on Linux and Mac

Conda instructions below. To install from source, go to the releases page, download the latest release, and install

```
git clone --recursive https://github.com/rwdavies/rrbgen.git
cd rrbgen
./scripts/install-r-dependencies.R
cd releases
wget https://github.com/rwdavies/rrbgen/releases/download/0.0.6/rrbgen_0.0.6.tar.gz ## or curl -O
R CMD INSTALL rrbgen_0.0.6.tar.gz
```

To install the latest development code in the repository, use `./scripts/build-and-install.sh`

To install an older release, either download an older release from the Github releases page, or use the older releases included with the repository in the `releases` folder

If you see errors like "error while loading shared libraries: libmpc.so.2: cannot open shared object file: No such file or directory", then please run the following before the R CMD INSTALL
```
./scripts/install-package-dependencies.sh
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:`pwd`/install/lib/
```
If you have any other installation problems or suggestions please report them as a github issue.

### Install using conda

rrbgen can be installed using [conda](https://conda.io/miniconda.html). Full tutorials can be found elsewhere, but briefly, something like this should work
```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
conda install rrbgen -c defaults -c bioconda -c conda-forge
source activate
R -e 'library("rrbgen")'
```

### Example commands in R
```
library("rrbgen")

## compare against ./external/bgen/example.gen
## see also examples in rrbgen/tests/testthat/test-read.R and test-write.R
bgen_file <- "./external/bgen/example.16bits.bgen"

## sample names only
sample_names <- rrbgen_load_samples(bgen_file)

## variant information only
var_info <- rrbgen_load_variant_info(bgen_file)

## load sample names, variant information, and genotype probabilities
out <- rrbgen_load(bgen_file)

## load a subset based on sample names and variant information
## note as implemented this requires a full pass over the data (although not decompression)
out <- rrbgen_load(bgen_file,
    vars_to_get = c("SNPID_2", "SNPID_10", "SNPID_3"),
    samples_to_get = c("sample_100", "sample_010", "sample_035")
)

```
