test_that("can write a simple bgen header", {

    L_SI <- 100 ## ignored in this context since no sample names yet
    free <- NULL ## ignored
    M <- 10
    N <- 20
    SampleIdentifiers <- 1
    Layout <- 2
    CompressedSNPBlocks <- 1

    bgen_file <- tempfile()
    to.write <- file(bgen_file, "wb")    
    L_H <- write_bgen_header(
        to.write,
        L_SI,
        M,
        N,
        SampleIdentifiers,
        Layout,
        CompressedSNPBlocks,
        free
    )
    close(to.write)

    to.read <- file(bgen_file, "rb")
    bgen_header <- load_bgen_header(to.read)
    close(to.read)

    ## check
    expect_equal(bgen_header$offset, L_H + L_SI) ## by definition
    expect_equal(bgen_header$L_H, 20)
    expect_equal(bgen_header$M, M)
    expect_equal(bgen_header$CompressedSNPBlocks, CompressedSNPBlocks)
    expect_equal(bgen_header$Layout, Layout)
    expect_equal(bgen_header$SampleIdentifiers, SampleIdentifiers)

})


test_that("can write a bgen file with header and sample names", {

    sample_names <- c("samp1", "jimmy445", "samp3")

    bgen_file <- tempfile()
    rrbgen_write(
        bgen_file,
        sample_names
    )
    sample_names_from_bgen <- rrbgen_load_samples(bgen_file)
    
    expect_equal(
        as.character(sample_names_from_bgen),
        sample_names
    )
    
})


test_that("can write a bgen file with header, sample names and SNP information", {

    sample_names <- c("samp1", "jimmy445", "samp3")
    var_info <- make_fake_var_info(12)
    CompressedSNPBlocks <- 1 ## compressed
    
    bgen_file <- tempfile()

    rrbgen_write(
        bgen_file,
        sample_names = sample_names,
        var_info = var_info
    )
    
    sample_names_from_bgen <- rrbgen_load_samples(bgen_file)
    expect_equal(
        as.character(sample_names_from_bgen),
        as.character(sample_names)
    )
    
    loaded_var_info <- rrbgen_load_variant_info(bgen_file)

    expect_equal(
        loaded_var_info,
        var_info
    )


})

test_that("can write a full bgen file", {

    skip("wer")
    
   library("testthat"); setwd("~/Dropbox/rrbgen/rrbgen/R/"); source("read-functions.R") ;source("write-functions.R"); source("test-drivers.R")
    
    sample_names <- c("edgar", "gsp", "silva", "lesnar")
    var_info <- make_fake_var_info(8)
    var_ids <- var_info[, "snpid"]

    set.seed(234)
    gp <- make_fake_gp(sample_names, var_ids, random_fraction = 0.05)

    bgen_file <- tempfile()

    ## for (CompressedSNPBlocks in c(0, 1)) {

        CompressedSNPBlocks <- 0
        rrbgen_write(
            bgen_file,
            sample_names = sample_names,
            var_info = var_info,
            gp = gp
        )

        out <- rrbgen_load(bgen_file)
        gp <- out$gp

        ## OK THIS RUNS
        ## NOW CHECK
        ## I AM NOT 100% SURE IF THERE IS AN OFF-BY-1 error in the offsets
        ## for the genotype in THIS BIT OF CODE
        ## Genotype block (Layout 2) (Minus C and D)
        

        ## CHECK!
    ## }

    
})
