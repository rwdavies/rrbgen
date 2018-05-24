## tests using external data from bgen library
## note - these are run in rrbgen/rrbgen/tests/testthat
## need to go back 3 directories
external_bgen_dir <- "../../../external/bgen/"
bgen_file <- file.path(external_bgen_dir, "example.16bits.bgen")
gen_file <- file.path(external_bgen_dir, "example.gen")    
gen <- read.table(gen_file)

test_that("can test something", {

    ## close(to.read)
    to.read <- file(bgen_file, "rb")
    bgen_header <- load_bgen_header(to.read)
    close(to.read)

    ## expect_equal(1, 2)    
    expect_equal(bgen_header$N, 500)

})
