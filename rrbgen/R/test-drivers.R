make_fake_var_info <- function(M) {
    var_info <- array(NA, c(M, 6))
    colnames(var_info) <- c("chr", "varid", "rsid", "position", "ref", "alt")
    var_info[, "chr"] <- as.character(12)
    var_info[, "varid"] <- paste0("affx-", round(100 * runif(M)))    
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

## gp is nice and all
## but STITCH will look like this
## list_of_gp_raw_t_hom_ref where rows = individuals & cols = snps
## list_of_gp_raw_t_het     where rows = individuals & cols = snps
convert_gp_to_list_of_raw <- function(gp, nCores = 1, B_bit_prob = 8) {
    M <- dim(gp)[[1]]
    N <- dim(gp)[[2]]
    ## nCores here doesn't mean cores
    ## but this is how this will be generated in STITCH
    sampleRange <- getSampleRange(N, nCores)
    ## this part emulates STITCH
    list_of_gp_raw_t <- parallel::mclapply(1:length(sampleRange), mc.cores = nCores, function(i_range) {
        whoR <- sampleRange[[i_range]]
        who <- whoR[1]:whoR[2]
        N2 <- length(who)
        to_out <- array(as.raw(0), c(2 * N2 * (B_bit_prob / 8), M))        
        ## now fill this in
        for(i_sample in 1:length(who)) {
            ## every SNP for this individual
            gp_t <- t(gp[, who[i_sample], ])
            rcpp_place_gp_t_into_output(
                gp_t,
                to_out,
                i_sample,
                nSNPs = ncol(gp_t),
                B_bit_prob
            )
        }
        return(to_out)        
    })
    return(list_of_gp_raw_t)
}


acceptable_tolerance <- function(B_bit_prob) {
    ## tolerance not quite working here for higher values
    ## my assumption of gen being "written" to bgen might not be right
    tolerance <- 1 / (2 ** B_bit_prob - 1) / 2 ## I think this is right
    tolerance <- as.numeric(c("8" = tolerance * 2, "16" = tolerance * 8, "24" = 5e-5, "32" = 5e-5)[as.character(B_bit_prob)])
    return(tolerance)
}


convert_gp_to_dosage <- function(gp, i_snp, normalize = TRUE) {
    dosage <- gp[i_snp, , 2] * 1 + gp[i_snp, , 3] * 2
    if (normalize) {
        dosage[is.na(dosage)] <- mean(dosage, na.rm = TRUE)
        dosage <- dosage - mean(dosage)
    }
    return(dosage)
}

make_fake_pheno_file <- function(gp, i_snp = 4, h2_g = 0.5) {
    dosage <- convert_gp_to_dosage(gp, i_snp)     
    y_snp <- dosage / sd(dosage) * sqrt(h2_g)
    y_e <- rnorm(dim(gp)[2], mean = 0, sd = sqrt(1 - h2_g))
    y <- y_snp + y_e
    ##
    pheno_file <- tempfile()
    m <- matrix(y, ncol = 1)
    colnames(m) <- "pheno1"
    write.table(
        m,
        file = pheno_file,
        row.names = FALSE,
        col.names = TRUE,
        quote = FALSE
    )
    return(
        list(
            pheno_file = pheno_file,
            y = y
        )
    )
}
