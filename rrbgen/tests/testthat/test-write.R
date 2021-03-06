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

    set.seed(234)    
    sample_names <- c("edgar", "gsp", "silva", "lesnar")
    CompressedSNPBlocks <- 1

    for(N in c(4, 1)) {
        for(M in c(8, 1)) {
            for(nCores in c(1, 4)) {

                var_info <- make_fake_var_info(M)
                var_ids <- var_info[, "varid"]
                
                set.seed(100)
                gp <- make_fake_gp(sample_names[1:N], var_ids, random_fraction = 0.05)
                
                bgen_file <- tempfile()
                
                for(B_bit_prob in c(8, 16, 24, 32)) {    
                    
                    rrbgen_write(
                        bgen_file,
                        sample_names = sample_names[1:N],
                        var_info = var_info,
                        gp = gp,
                        CompressedSNPBlocks = CompressedSNPBlocks,
                        B_bit_prob = B_bit_prob,
                        nCores = nCores
                    )
                
                    out <- rrbgen_load(bgen_file)
                    loaded_gp <- out$gp
                    
                    expect_equal(dimnames(gp)[[1]], dimnames(loaded_gp)[[1]])
                    expect_equal(dimnames(gp)[[2]], as.character(dimnames(loaded_gp)[[2]])) ## argh
                    expect_equal(dimnames(gp)[[3]], dimnames(loaded_gp)[[3]])
                    tolerance <- acceptable_tolerance(B_bit_prob) 
                    
                    expect_equal(sum(abs(gp - loaded_gp) > tolerance, na.rm = TRUE), 0)
                    
                    expect_equal(as.logical(is.na(gp) ), as.logical(is.na(loaded_gp)))
                    
                }

            }
        }

    }
        
})


## intended write case for STITCH
## lists of gp_raw come from different cores
## at the end, write intelligently
## note that list_of_gp_raw_t
## is a list, where each entry is a matrix
## where SNPs are the columns and subjects are the individuals
## where for B_bit_prob and B_bit_prob_divide_8
## where list_of_gp_raw_t[[1]][1:B_bit_prob_divide_8, 1]
##       are the hom ref bits, and 
## where list_of_gp_raw_t[[1]][B_bit_prob_divide_8 + 1:B_bit_prob_divide_8, 1]
##       are the het bits
## and so forth
test_that("writing a full bgen file is the same using either gp or list of gp_raw ", {

    ##   library("testthat"); setwd("~/Dropbox/rrbgen/rrbgen/R/"); source("read-functions.R") ;source("write-functions.R"); source("from-stitch.R"); source("test-drivers.R"); library("rrbgen")

    set.seed(234)    
    sample_names <- paste0("samp", 1:40)
    var_info <- make_fake_var_info(100)
    var_ids <- var_info[, "varid"]
    CompressedSNPBlocks <- 1
    
    gp <- make_fake_gp(sample_names, var_ids, random_fraction = 0)
    B_bit_prob <- 24
    gp[1, 1, ] <- c(0.985363562311, 0.014634971374, 0.000001466315)

    for(nCores in c(1, 4)) {
        for(B_bit_prob in c(8, 16, 24, 32)) {
        
            list_of_gp_raw_t <- convert_gp_to_list_of_raw(
                gp,
                nCores = nCores,
                B_bit_prob = B_bit_prob
            )
            
            bgen_file_list_of_gp_raw_t <- tempfile()
            bgen_file_gp <- tempfile()    
            
            rrbgen_write(
                bgen_file_gp,
                gp = gp,
                sample_names = sample_names,
                var_info = var_info,
                CompressedSNPBlocks = CompressedSNPBlocks,
                B_bit_prob = B_bit_prob
            )
            
            rrbgen_write(
                bgen_file_list_of_gp_raw_t,
                sample_names = sample_names,
                var_info = var_info,
                list_of_gp_raw_t = list_of_gp_raw_t,        
                CompressedSNPBlocks = CompressedSNPBlocks,
                B_bit_prob = B_bit_prob
            )
            
            out_gp <- rrbgen_load(bgen_file_gp)
            out_list_of_gp_raw_t <- rrbgen_load(bgen_file_list_of_gp_raw_t)
            
            expect_equal(out_list_of_gp_raw_t, out_gp)
            ## floating point numbers don't "have" to be within the bounds
            ## but rounded in R should
            expect_equal(0 <= min(round(out_gp$gp)), TRUE)
            expect_equal(max(round(out_gp$gp)) <= 1, TRUE)        
            
        }
    }
        
})


test_that("can write to bgen in stages", {

    ##   library("testthat"); setwd("~/Dropbox/rrbgen/rrbgen/R/"); source("read-functions.R") ;source("write-functions.R"); source("from-stitch.R"); source("test-drivers.R"); library("rrbgen")    
    set.seed(824)    
    sample_names <- paste0("samp", 1:40)
    n_snps <- 100
    var_info <- make_fake_var_info(n_snps)
    var_ids <- var_info[, "varid"]
    
    gp <- make_fake_gp(sample_names, var_ids, random_fraction = 0)
    B_bit_prob <- 24
    gp[1, 1, ] <- c(0.985363562311, 0.014634971374, 0.000001466315)

    for(nCores in c(1, 4)) {
        for(B_bit_prob in c(8, 16, 24, 32)) {

            list_of_gp_raw_t <- convert_gp_to_list_of_raw(
                gp,
                nCores = nCores,
                B_bit_prob = B_bit_prob
            )
            
            bgen_file_full <- tempfile()
            bgen_file_parts <- tempfile()    
            
            ## subset on the fly
            f <- function(w) {
                lapply(list_of_gp_raw_t, function(x) x[, w])
            }

            ## write normally
            rrbgen_write(
                bgen_file_full,
                list_of_gp_raw_t = f(1:100),
                sample_names = sample_names,
                var_info = var_info[1:100, ],
                B_bit_prob = B_bit_prob
            )

            ## write header
            out <- rrbgen_write(
                bgen_file_parts,
                sample_names = sample_names,
                B_bit_prob = B_bit_prob,
                close_bgen_file = FALSE,
                header_M = 100
            )
            bgen_file_connection <- out$bgen_file_connection
            previous_offset <- out$final_binary_length
            
            ## write first third
            out <- rrbgen_write(
                bgen_file_connection = bgen_file_connection,
                previous_offset = previous_offset,
                add_to_bgen_connection = TRUE,
                close_bgen_file = FALSE,
                list_of_gp_raw_t = f(1:10),
                sample_names = sample_names,
                var_info = var_info[1:10, ],
                B_bit_prob = B_bit_prob
            )
            bgen_file_connection <- out$bgen_file_connection
            previous_offset <- out$final_binary_length

            ## write second half
            out2 <- rrbgen_write(
                bgen_file_connection = bgen_file_connection,
                previous_offset = previous_offset,
                add_to_bgen_connection = TRUE,
                close_bgen_file = FALSE,
                list_of_gp_raw_t = f(11:20),
                sample_names = sample_names,
                var_info = var_info[11:20, ],
                B_bit_prob = B_bit_prob
            )
            previous_offset <- out2$final_binary_length
            
            ## write third part
            rrbgen_write(
                bgen_file_connection = bgen_file_connection,
                previous_offset = previous_offset,
                add_to_bgen_connection = TRUE,
                close_bgen_file = TRUE,
                list_of_gp_raw_t = f(21:100),
                sample_names = sample_names,
                var_info = var_info[21:100, ],
                B_bit_prob = B_bit_prob
            )

            out_full <- rrbgen_load(bgen_file_full)
            out_parts <- rrbgen_load(bgen_file = bgen_file_parts)
            
            expect_equal(out_full, out_parts)
            ## floating point numbers don't "have" to be within the bounds
            ## but rounded in R should
            expect_equal(0 <= min(round(out_parts$gp)), TRUE)
            expect_equal(max(round(out_parts$gp)) <= 1, TRUE)        
            
        }
    }
        
})

