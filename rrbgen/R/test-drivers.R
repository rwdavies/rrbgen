make_fake_var_info <- function(M) {
    var_info <- array(NA, c(M, 6))
    colnames(var_info) <- c("chr", "snpid", "rsid", "position", "ref", "alt")
    var_info[, "chr"] <- as.character(12)
    var_info[, "snpid"] <- paste0("affx-", round(100 * runif(M)))    
    var_info[, "rsid"] <- paste0("rs", round(100 * runif(M)))
    var_info[, "position"] <- 1:M
    var_info[, "ref"] <- sample(c("A", "C", "G", "T"), M, replace = TRUE)
    var_info[, "alt"] <- c("A", "C", "G", "T")[5 - match(var_info[, "ref"], c("A", "C", "G", "T"))]
    return(var_info)
}

make_fake_gp <- function(sample_names, var_ids, random_fraction = 0.05) {
    N <- length(sample_names)
    M <- length(var_ids)
    gp <- array(runif(M * N * 3), c(M, N, 3))
    dimnames(gp)[[1]] <- var_ids
    dimnames(gp)[[2]] <- sample_names
    dimnames(gp)[[3]] <- c("hom_ref", "het", "hom_alt")
    ## now fill with random! include some NAs
    for(i_snp in 1:M) {
        for(i_sample in 1:N) {
            if (runif(1) < random_fraction) {
                gp[i_snp, i_sample, ] <- NA
            } else {
                gp[i_snp, i_sample, ] <- gp[i_snp, i_sample, ] / sum(gp[i_snp, i_sample, ])
            }
        }
    }
    ##
    return(gp)
}


acceptable_tolerance <- function(B_bit_prob) {
    ## tolerance not quite working here for higher values
    ## my assumption of gen being "written" to bgen might not be right
    tolerance <- 1 / (2 ** B_bit_prob - 1) / 2 ## I think this is right
    tolerance <- as.numeric(c("8" = tolerance * 2, "16" = tolerance * 8, "24" = 5e-5, "32" = 5e-5)[as.character(B_bit_prob)])
    return(tolerance)
}
