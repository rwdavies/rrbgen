test_that("writing, reading, then writing is equivalent in R", {

    gp_sub <- array(NA, c(9, 3))
    dimnames(gp_sub)[[1]] <- list("sam", "jim", "jon", "bar", "bill", "fred", "missy", "weirdy", "weirdier")
    dimnames(gp_sub)[[2]] <- list("hom_ref", "het", "hom_alt")
    gp_sub[1, ] <- c(0, 0, 1)
    gp_sub[2, ] <- c(1, 0, 0)
    gp_sub[3, ] <- c(0.5, 0.5, 0)
    gp_sub[4, ] <- c(0.5, 0, 0.5)
    gp_sub[5, ] <- c(0.1, 0.2, 0.7)
    gp_sub[6, ] <- c(0.7, 0.1, 0.2)
    gp_sub[7, ] <- c(NA, NA, NA)
    gp_sub[8, ] <- c(0.4721318,  0.2002350,  0.3276332 )
    gp_sub[9, ] <- c(0.985363562311, 0.014634971374, 0.000001466315)
    N <- 9
    is_missing <- array(FALSE, N)
    is_missing[7] <- TRUE    

    for(B_bit_prob in c(8, 16, 24, 32)) {
        data_raw_for_probs <- make_raw_data_vector_for_probabilities(
            gp_sub = gp_sub,
            B_bit_prob  = B_bit_prob
        )
        
        gp_written <- convert_raw_probabilities_to_double_probabilities(
            data_raw_for_probs,
            N,
            B_bit_prob,
            is_missing
        )
        
        data_raw_for_probs2 <- make_raw_data_vector_for_probabilities(
            gp_sub = gp_written,
            B_bit_prob  = B_bit_prob
        )
        
        expect_equal(
            data_raw_for_probs,
            data_raw_for_probs2
        )

    }
    
})


test_that("can write matrix of triplet of probabilities to raw vector", {

    gp_sub <- array(NA, c(9, 3))
    dimnames(gp_sub)[[1]] <- list("sam", "jim", "jon", "bar", "bill", "fred", "missy", "weirdy", "weirdier")
    dimnames(gp_sub)[[2]] <- list("hom_ref", "het", "hom_alt")
    gp_sub[1, ] <- c(0, 0, 1)
    gp_sub[2, ] <- c(1, 0, 0)
    gp_sub[3, ] <- c(0.5, 0.5, 0)
    gp_sub[4, ] <- c(0.5, 0, 0.5)
    gp_sub[5, ] <- c(0.1, 0.2, 0.7)
    gp_sub[6, ] <- c(0.7, 0.1, 0.2)
    gp_sub[7, ] <- c(NA, NA, NA)
    gp_sub[8, ] <- c(0.4721318, 0.2002350, 0.3276332)
    gp_sub[9, ] <- c(0.985363562311, 0.014634971374, 0.000001466315)    
    N <- 9
    is_missing <- array(FALSE, N)
    is_missing[7] <- TRUE        
    
    for(B_bit_prob in c(8, 16, 24, 32)) {
            
        v <- make_raw_data_vector_for_probabilities(
            gp_sub = gp_sub,
            B_bit_prob  = B_bit_prob
        )

        ## convert to work on 
        v2 <- rcpp_make_raw_data_vector_for_probabilities(
            gp_sub = gp_sub,
            B_bit_prob  = B_bit_prob
        )
        
        expect_equal(v, v2)
    }
    
})


test_that("can read raw to form probabilities", {

    gp_sub <- array(NA, c(9, 3))
    dimnames(gp_sub)[[1]] <- list("sam", "jim", "jon", "bar", "bill", "fred", "missy", "weirdy", "weirder")
    dimnames(gp_sub)[[2]] <- list("hom_ref", "het", "hom_alt")
    gp_sub[1, ] <- c(0, 0, 1)
    gp_sub[2, ] <- c(1, 0, 0)
    gp_sub[3, ] <- c(0.5, 0.5, 0)
    gp_sub[4, ] <- c(0.5, 0, 0.5)
    gp_sub[5, ] <- c(0.1, 0.2, 0.7)
    gp_sub[6, ] <- c(0.7, 0.1, 0.2)
    gp_sub[7, ] <- c(NA, NA, NA)
    gp_sub[8, ] <- c(0.4721318, 0.2002350, 0.3276332)
    gp_sub[9, ] <- c(0.985363562311, 0.014634971374, 0.000001466315)        
    N <- 9
    is_missing <- array(FALSE, N)
    is_missing[7] <- TRUE
    
    for(B_bit_prob in c(8, 16, 24, 32)) {
            
        data_raw_for_probs <- make_raw_data_vector_for_probabilities(
            gp_sub = gp_sub,
            B_bit_prob  = B_bit_prob
        )
        
        gp1 <- convert_raw_probabilities_to_double_probabilities(
            data_raw_for_probs,
            N,
            B_bit_prob,
            is_missing
        )

        gp2 <- rcpp_convert_raw_probabilities_to_double_probabilities(
            data_raw_for_probs,
            N,
            B_bit_prob,
            is_missing
        )

        expect_equal(gp1, gp2)
        expect_equal((min(gp1, na.rm = TRUE) + 1e-8) > 0, TRUE)
        expect_equal((max(gp2, na.rm = TRUE) - 1e-8) < 1, TRUE)        
        
    }
    
})
