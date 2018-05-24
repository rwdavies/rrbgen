## tests using external data from bgen library
## note - these are run in rrbgen/rrbgen/tests/testthat
## need to go back 3 directories
## library('testthat'); setwd("~/Dropbox/rrbgen/rrbgen/tests/testthat")
external_bgen_dir <- "../../../external/bgen/"
bgen_file <- file.path(external_bgen_dir, "example.16bits.bgen")
gen_file <- file.path(external_bgen_dir, "example.gen")    
sample_file <- file.path(external_bgen_dir, "example.sample")
gen <- read.table(gen_file)
sample <- read.table(sample_file, sep = " ")
sample_ids_from_sample_file <- as.character(sample[-c(1:2), ])

## source("~/Dropbox/rrbgen/rrbgen/R/functions.R")

test_that("can test things", {

    ## source("~/Dropbox/rrbgen/rrbgen/R/functions.R");     close(to.read)
    to.read <- file(bgen_file, "rb")

    ## header block
    bgen_header <- load_bgen_header(to.read)
    Layout <- bgen_header$Layout
    CompressedSNPBlocks <- bgen_header$CompressedSNPBlocks
    L_H <- bgen_header$L_H
    M <- bgen_header$M
    N <- bgen_header$N
    SampleIdentifiers <- bgen_header$SampleIdentifiers
    expect_equal(bgen_header$N, 500)

    ## sample identifier block
    out <- load_bgen_sample_identifier_block(to.read, L_H + 4, SampleIdentifiers, N)
    sample_names <- out$sample_names
    expect_equal(
        (sample_ids_from_sample_file),
        as.character(sample_names)
    )

    ## variant
    out <- load_variant_identifying_data_for_one_snp(to.read, Layout)
    variant_info <- out$variant_info
    num_K_alleles <- out$num_K_alleles
    expect_equal(as.integer(variant_info["chr"]), gen[1, 1])
    expect_equal(as.character(variant_info["snpid"]), as.character(gen[1, 2]))
    expect_equal(as.character(variant_info["rsid"]), as.character(gen[1, 3]))
    expect_equal(as.character(variant_info["position"]), as.character(gen[1, 4]))
    expect_equal(as.character(variant_info["ref"]), as.character(gen[1, 5]))
    expect_equal(as.character(variant_info["alt"]), as.character(gen[1, 6]))

    ## genotypes
    out <- load_genotypes_for_one_snp(to.read, num_K_alleles, N, CompressedSNPBlocks) 

    tolerance <- 1e-4
    ## check a few 
    expect_equal(
        is.na(out$gen_probs[, 1]),
        as.logical(gen[1, 7:9] == 0)
    )
    ## 
    expect_equal(
        as.logical(abs(out$gen_probs[, 2] - gen[1, 10:12]) > tolerance),
        rep(FALSE, 3)
    )
    ## 
    expect_equal(
        as.logical(abs(out$gen_probs[, 3] - gen[1, 13:15]) > tolerance),
        rep(FALSE, 3)
    )

    close(to.read)

})


test_that("can load sample names", {
    
    sample_names <- rrbgen_load_samples(bgen_file)
    expect_equal(
        (sample_ids_from_sample_file),
        as.character(sample_names)
    )

})
