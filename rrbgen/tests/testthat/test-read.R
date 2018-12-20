## tests using external data from bgen library
## note - these are run in rrbgen/rrbgen/tests/testthat
## need to go back 3 directories
## library("testthat"); setwd("~/Dropbox/rrbgen/rrbgen/R/"); source("from-stitch.R"); source("read-functions.R") ;source("write-functions.R"); source("test-drivers.R");   library('testthat'); setwd("~/Dropbox/rrbgen/rrbgen/tests/testthat"); library("rrbgen")
external_bgen_dir <- "../../../external/bgen/"
bits <- c(8, 16, 24, 32)
external_bgen_files <- sapply(bits, function(bit) {
    file.path(external_bgen_dir, paste0("example.", bit, "bits.bgen"))
})
names(external_bgen_files) <- as.character(bits)
gen_file <- file.path(external_bgen_dir, "example.gen")    
sample_file <- file.path(external_bgen_dir, "example.sample")
gen <- read.table(gen_file)
sample <- read.table(sample_file, sep = " ")
sample_ids_from_sample_file <- as.character(sample[-c(1:2), ])


## check for files, if they do not exist, do not run these tests
## these tests are meant to be run with external data
skip_if_no_files <- function() {
    for(file in external_bgen_files) {
        if (file.exists(file) == FALSE) {
            skip("External files not available")
        }
    }
}


test_that("can test things", {

    skip_if_no_files()
    
    ##
    ## library("testthat"); setwd("~/Dropbox/rrbgen/rrbgen/R/"); source("read-functions.R") ;source("write-functions.R"); source("test-drivers.R");  setwd("~/Dropbox/rrbgen/rrbgen/tests/testthat")    ;      close(to.read); library("rrbgen")
    to.read <- file(external_bgen_files["16"], "rb")

    ## header block
    bgen_header <- load_bgen_header(to.read)
    Layout <- bgen_header$Layout
    CompressedSNPBlocks <- bgen_header$CompressedSNPBlocks
    L_H <- bgen_header$L_H
    M <- bgen_header$M
    N <- bgen_header$N
    SampleIdentifiers <- bgen_header$SampleIdentifiers
    offset <- bgen_header$offset
    expect_equal(bgen_header$N, 500)

    ## sample identifier block
    out <- load_bgen_sample_identifier_block(to.read, L_H + 4, SampleIdentifiers, N)
    sample_names <- out$sample_names
    L_SI <- out$L_SI
    expect_equal(
        (sample_ids_from_sample_file),
        as.character(sample_names)
    )
    
    ## variant
    out <- load_variant_identifying_data_for_one_snp(to.read, L_SI + L_H + 4, Layout, CompressedSNPBlocks)
    variant_info <- out$variant_info
    num_K_alleles <- out$num_K_alleles
    L_vid <- out$L_vid
    C <- out$C
    D <- out$D
    expect_equal(as.integer(variant_info["chr"]), gen[1, 1])
    expect_equal(as.character(variant_info["varid"]), as.character(gen[1, 2]))
    expect_equal(as.character(variant_info["rsid"]), as.character(gen[1, 3]))
    expect_equal(as.character(variant_info["position"]), as.character(gen[1, 4]))
    expect_equal(as.character(variant_info["ref"]), as.character(gen[1, 5]))
    expect_equal(as.character(variant_info["alt"]), as.character(gen[1, 6]))

    ## genotypes - last two 4's are for C and D
    out <- load_genotypes_for_one_snp(
        to.read,
        binary_start = (offset + 4) + L_vid + 4 + 4,
        num_K_alleles = num_K_alleles,
        N = N,
        CompressedSNPBlocks = CompressedSNPBlocks,
        C = C
    ) 

    tolerance <- 1e-4
    ## check a few 
    expect_equal(
        is.na(out$gen_probs[1, ]),
        as.logical(gen[1, 7:9] == 0)
    )
    ## 
    expect_equal(
        as.logical(abs(out$gen_probs[2, ] - gen[1, 10:12]) > tolerance),
        rep(FALSE, 3)
    )
    ## 
    expect_equal(
        as.logical(abs(out$gen_probs[3, ] - gen[1, 13:15]) > tolerance),
        rep(FALSE, 3)
    )

    close(to.read)

})


test_that("can load sample names", {
    
    skip_if_no_files()
    
    sample_names <- rrbgen_load_samples(external_bgen_files["16"])
    expect_equal(
        (sample_ids_from_sample_file),
        as.character(sample_names)
    )

})

test_that("can load SNPs", {
    
    skip_if_no_files()
    
    var_info <- rrbgen_load_variant_info(external_bgen_files["16"])

    expect_equal(as.integer(var_info[, 1]), gen[, 1])
    for(i in 2:6) {
        expect_equal(var_info[, i], as.character(gen[, i]))
    }

})

test_that("can load everything", {
    
    skip_if_no_files()
    
    for(B_bit_prob in c(8, 16, 24, 32)) {

        out <- rrbgen_load(external_bgen_files[as.character(B_bit_prob)])
        gp <- out$gp

        tolerance <- acceptable_tolerance(B_bit_prob)         

        gen1 <- gen[, 6 + seq(1, ncol(gen) - 6, 3)]
        gen2 <- gen[, 6 + seq(2, ncol(gen) - 6, 3)]
        gen3 <- gen[, 6 + seq(3, ncol(gen) - 6, 3)]

        expect_equal(max(abs(gen1 - gp[, , 1]), na.rm = TRUE) < tolerance, TRUE)
        expect_equal(max(abs(gen2 - gp[, , 2]), na.rm = TRUE) < tolerance, TRUE)
        expect_equal(max(abs(gen3 - gp[, , 3]), na.rm = TRUE) < tolerance, TRUE)
        
        ## check missing
        expect_equal(is.na(gp[, , 1]), (gen1 == 0) & (gen2 == 0) & (gen3 == 0), check.attributes = FALSE)

    }

})


test_that("can read, write then re-read a file", {

    skip_if_no_files()
    
    B_bit_prob <- 16
    out <- rrbgen_load(external_bgen_files[as.character(B_bit_prob)])
    gp <- out$gp
    sample_names <- out$sample_names
    var_info <- out$var_info
    CompressedSNPBlocks <- 1

    bgen_file <- tempfile()
    rrbgen_write(
        bgen_file,
        sample_names = sample_names,
        var_info = var_info,
        gp = gp,
        CompressedSNPBlocks = CompressedSNPBlocks,
        B_bit_prob = B_bit_prob
    )

    out2 <- rrbgen_load(bgen_file)
    expect_equal(out2$sample_names, out$sample_names)
    expect_equal(out2$var_info, out$var_info)
    expect_equal(out2$gp, out$gp)    ## hmm, exactly the same? good, yes?

})


test_that("requesting a SNP that does not exist throws an error", {

    skip_if_no_files()
    
    var_info <- rrbgen_load_variant_info(external_bgen_files["16"])
    which_vars <- 5:8
    B_bit_prob <- 8
    expect_error(
        rrbgen_load(
            bgen_file = external_bgen_files[as.character(B_bit_prob)],
            vars_to_get = c("betsy", var_info[which_vars, "varid"], "jimmy")
        )
    )

})


test_that("can load subset of SNPs by varid", {

    skip_if_no_files()
    
    var_info <- rrbgen_load_variant_info(external_bgen_files["16"])
    which_vars <- c(5, 8, 4, 6)
    B_bit_prob <- 8    

    for(B_bit_prob in c(8, 16, 24, 32)) {
    
        out <- rrbgen_load(
            bgen_file = external_bgen_files[as.character(B_bit_prob)],
            vars_to_get = var_info[which_vars, "varid"]
        )
        tolerance <- acceptable_tolerance(B_bit_prob)         
        gp <- out$gp

        expect_equal(as.integer(out$var_info[, 1]), gen[which_vars, 1])

        gen1 <- gen[which_vars, 6 + seq(1, ncol(gen) - 6, 3)]
        gen2 <- gen[which_vars, 6 + seq(2, ncol(gen) - 6, 3)]
        gen3 <- gen[which_vars, 6 + seq(3, ncol(gen) - 6, 3)]
        
        expect_equal(max(abs(gen1 - gp[, , 1]), na.rm = TRUE) < tolerance, TRUE)
        expect_equal(max(abs(gen2 - gp[, , 2]), na.rm = TRUE) < tolerance, TRUE)
        expect_equal(max(abs(gen3 - gp[, , 3]), na.rm = TRUE) < tolerance, TRUE)
        
        ## check missing
        expect_equal(is.na(gp[, , 1]), (gen1 == 0) & (gen2 == 0) & (gen3 == 0), check.attributes = FALSE)

    }

})

test_that("requesting a sample that does not exist throws an error", {

    skip_if_no_files()
    
    expect_error(
        rrbgen_load(
            bgen_file = external_bgen_files[as.character(8)],
            samples_to_get = c("betsy", "jimmy")
        ),
        "Cannot find some of the samples. The first such missing samples are:betsy, jimmy"
    )

})


test_that("requesting a sample more than once throws an error", {

    skip_if_no_files()
    
    expect_error(
        rrbgen_load(
            bgen_file = external_bgen_files[as.character(8)],
            samples_to_get = c("samp1", "samp2", "samp1")
        ),
        "Please do not request the same sample more than once. The following sample was requested 2 times:samp1"
    )

})

test_that("can load subset of samples", {

    skip_if_no_files()
    
    sample_names <- rrbgen_load_samples(external_bgen_files["16"])
    which_samples <- c(15, 92, 101, 45)
    which_samples <- 1:4
    samples_to_get <- sample_names[which_samples]
    B_bit_prob <- 8    

    for(B_bit_prob in c(8, 16, 24, 32)) {
    
        out <- rrbgen_load(
            bgen_file = external_bgen_files[as.character(B_bit_prob)],
            samples_to_get = samples_to_get
        )
        tolerance <- acceptable_tolerance(B_bit_prob)         
        gp <- out$gp

        expect_equal(as.integer(out$var_info[, 1]), gen[, 1])

        gen1 <- gen[, 6 + 3 * which_samples - 2]
        gen2 <- gen[, 6 + 3 * which_samples - 1]
        gen3 <- gen[, 6 + 3 * which_samples - 0]
        
        expect_equal(max(abs(gen1 - gp[, , 1]), na.rm = TRUE) < tolerance, TRUE)
        expect_equal(max(abs(gen2 - gp[, , 2]), na.rm = TRUE) < tolerance, TRUE)
        expect_equal(max(abs(gen3 - gp[, , 3]), na.rm = TRUE) < tolerance, TRUE)
        
        ## check missing
        expect_equal(is.na(gp[, , 1]), (gen1 == 0) & (gen2 == 0) & (gen3 == 0), check.attributes = FALSE)

    }

})

test_that("can load subset of samples and SNPs", {

    skip_if_no_files()
    
    sample_names <- rrbgen_load_samples(external_bgen_files["16"])
    which_samples <- c(15, 92, 101, 45)
    samples_to_get <- sample_names[which_samples]    
    var_info <- rrbgen_load_variant_info(external_bgen_files["16"])
    which_vars <- c(5, 8, 4, 6)
    vars_to_get <- var_info[which_vars, "varid"]

    B_bit_prob <- 8

    for(B_bit_prob in c(8, 16, 24, 32)) {
    
        out <- rrbgen_load(
            bgen_file = external_bgen_files[as.character(B_bit_prob)],
            vars_to_get = vars_to_get,
            samples_to_get = samples_to_get
        )
        tolerance <- acceptable_tolerance(B_bit_prob)         
        gp <- out$gp

        expect_equal(as.integer(out$var_info[, 1]), gen[which_vars, 1])

        gen1 <- gen[which_vars, 6 + 3 * which_samples - 2]
        gen2 <- gen[which_vars, 6 + 3 * which_samples - 1]
        gen3 <- gen[which_vars, 6 + 3 * which_samples - 0]
        
        expect_equal(max(abs(gen1 - gp[, , 1]), na.rm = TRUE) < tolerance, TRUE)
        expect_equal(max(abs(gen2 - gp[, , 2]), na.rm = TRUE) < tolerance, TRUE)
        expect_equal(max(abs(gen3 - gp[, , 3]), na.rm = TRUE) < tolerance, TRUE)
        
        ## check missing
        expect_equal(is.na(gp[, , 1]), (gen1 == 0) & (gen2 == 0) & (gen3 == 0), check.attributes = FALSE)

    }

})
