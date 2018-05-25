test_that("can write a simple bgen header", {

    offset <- 10
    free <- NULL
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
