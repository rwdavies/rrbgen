rrbgen
======
**__Current Version: 0.0.4__**
Release date: July 13, 2018

A lightweight limited functionality R bgen read/write library

[![Build Status](https://travis-ci.org/rwdavies/rrbgen.svg)](https://travis-ci.org/rwdavies/rrbgen)

This library supports "v1.3" of bgen defined [here](http://www.well.ox.ac.uk/~gav/bgen_format/). It supports reading and writing using 8, 16, 24 or 32 bits per probability, using Layout = 2 and CompressedSNPBlocks = 1, for bi-allelic SNPs with samples of ploidy 2. Any other format specifications may crash unexpectedly without properly defined error.

### Quick start on Linux and Mac

```
git clone --recursive https://github.com/rwdavies/rrbgen.git
cd rrbgen
./scripts/install-r-dependencies.R
R CMD INSTALL ./releases/rrbgen_0.0.4.tar.gz
```
If you see errors like "error while loading shared libraries: libmpc.so.2: cannot open shared object file: No such file or directory", then please run the following before the R CMD INSTALL
```
./scripts/install-package-dependencies.sh
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:`pwd`/install/lib/
```
If you have any other installation problems or suggestions please report them as a github issue.

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
