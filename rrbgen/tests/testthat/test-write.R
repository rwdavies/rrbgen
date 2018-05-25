test_that("can write a simple bgen header", {

    offset <- 0 ## meaningless in this context
    free <- NULL ## ignored
    L_H <- 20 ## free being NULL makes this 20 by default
    M <- 10
    N <- 20
    SampleIdentifiers <- 1
    Layout <- 2
    CompressedSNPBlocks <- 1

    bgen_file <- tempfile()
    to.write <- file(bgen_file, "wb")    
    write_bgen_header(
        to.write,
        offset,
        M,
        N,
        SampleIdentifiers,
        Layout,
        CompressedSNPBlocks,
        free = NULL
    )
    close(to.write)

    to.read <- file(bgen_file, "rb")
    bgen_header <- load_bgen_header(to.read)
    close(to.read)

    ## check
    expect_equal(bgen_header$offset, offset)
    expect_equal(bgen_header$L_H, 20)
    expect_equal(bgen_header$M, M)
    expect_equal(bgen_header$CompressedSNPBlocks, CompressedSNPBlocks)
    expect_equal(bgen_header$Layout, Layout)
    expect_equal(bgen_header$SampleIdentifiers, SampleIdentifiers)

})


test_that("can write a bgen file with header and sample names", {

    sample_names <- c("samp1", "jimmy445", "samp3")
    offset <- 0 ## still meaningless in this context
    M <- 10 ## still meaningless in this context
    N <- length(sample_names)
    L_H <- 20 ## since free is NULL
    
    bgen_file <- tempfile()
    to.write <- file(bgen_file, "wb")
    ## header
    write_bgen_header(
        to.write,
        offset,
        M,
        N
    )
    write_bgen_sample_identifier_block(
        to.write,
        binary_start = L_H + 4,
        sample_names
    )
    close(to.write)

    sample_names_from_bgen <- rrbgen_load_samples(bgen_file)
    
    expect_equal(
        as.character(sample_names_from_bgen),
        sample_names
    )
    

})
