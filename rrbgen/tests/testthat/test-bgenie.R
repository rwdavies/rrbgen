bgenie_command <- "bgenie_v1.3_static1"
skip_bgenie_tests <- Sys.which("bgenie_v1.3_static1") == ""

if (1 == 0) {
    
    dir <- "/data/smew1/rdavies/rrbgen/rrbgen/R/"
    ##dir <- "~/personal/proj/rrbgen/rrbgen/R/"
    setwd(dir)
    library("testthat");    source("read-functions.R") ;source("write-functions.R"); source("test-drivers.R"); source("from-stitch.R"); library("rrbgen")

}


test_that("can use bgenie to generate correct p-values", {

    if (skip_bgenie_tests) {
        skip("bgenie not in path")
    }

    N <- 100
    M <- 10
    i_snp <- 4 ## this is the causal SNP
    h2_g <- 0.5 ## approx heritability of SNP
    CompressedSNPBlocks <- 1
    random_fraction <- 0 ## bgenie is NOT supported for missing data at this time
    
    set.seed(234)    
    sample_names <- paste0("samp", 1:N)
    var_info <- make_fake_var_info(M)
    var_ids <- var_info[, "varid"]
    
    gp <- make_fake_gp(sample_names, var_ids, random_fraction = random_fraction)
    out <- make_fake_pheno_file(gp, i_snp = i_snp, h2_g = h2_g)
    pheno_file <- out$pheno_file
    y <- out$y
        
    for(B_bit_prob  in c(8, 16, 24, 32)) {

        bgen_file <- tempfile()
        
        rrbgen_write(
            bgen_file,
            sample_names = sample_names,
            var_info = var_info,
            gp = gp,
            CompressedSNPBlocks = CompressedSNPBlocks,
            B_bit_prob = B_bit_prob
        )
        
        out <- rrbgen_load(bgen_file)
        gp_loaded <- out$gp
        
        dosage_file <- tempfile()
        out <- system(paste0(
            bgenie_command,
            " --bgen ", bgen_file,
            " --dosage ", dosage_file        
        ), intern = TRUE)
        
        out_file <- tempfile()
        out <- system(paste0(
            bgenie_command,
            " --bgen ", bgen_file,
            " --pheno ", pheno_file,    
            " --out ", out_file
        ), intern = TRUE)
        
        stats <- read.table(paste0(out_file, ".gz"), header = TRUE)
        bgenie_dosage <- read.table(paste0(dosage_file, ".gz"))        
        
        dosage <- convert_gp_to_dosage(gp_loaded, i_snp, normalize = FALSE)
        check <- dosage - bgenie_dosage[i_snp, -c(1:7)]
        expect_equal(max(abs(check), na.rm = TRUE) > 1e-4, FALSE) ## close enough
        
        ## now check coefficients of other SNPs
        for(i_snp in 1:M) {
            dosage <- convert_gp_to_dosage(gp_loaded, i_snp, normalize = FALSE)    
            c <- coefficients(summary(lm(y ~ dosage)))
            expect_equal(abs(stats[i_snp, "pheno1_beta"] - c["dosage", "Estimate"]) > 1e-2, FALSE)
            expect_equal(abs(stats[i_snp, "pheno1_se"] - c["dosage", "Std. Error"]) > 1e-2, FALSE)
            expect_equal(abs(stats[i_snp, "pheno1_t"] - c["dosage", "t value"]) > 1e-1, FALSE)
        }

    }
    
})
