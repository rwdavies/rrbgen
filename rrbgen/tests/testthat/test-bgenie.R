bgenie_command <- "bgenie_v1.3_static1"
skip_bgenie_tests <- Sys.which("bgenie_v1.3_static1") == ""

test_that("can use bgenie to generate correct p-values", {

    if (skip_bgenie_tests)
        skip("bgenie not in path")

    expect_equal(1, 0)
    
})
