#' @title Load entire bgen file
#' @param bgen_file Path to bgen file to load
#' @param gp_names_col Which column to use to label the variants in the first dimension of the genotype probabilities. Note that bgen has two variant ID columns, a variant identifier (here labelled varid) and rsid (here labelled rsid)
#' @param vars_to_get Only load these variants from file (from varid column). Note that this requires a complete pass over the data, although this should be efficient, as nothing is decompressed
#' @param samples_to_get Only load these samples from file. Note that this requires a complete pass over the data, although this should be efficient, as nothing is decompressed
#' @return gp, sample_names and var_info
#' @author Robert Davies
#' @export
rrbgen_load <- function(
    bgen_file,
    gp_names_col = "varid",
    vars_to_get = NULL,
    samples_to_get = NULL
) {
    if (! (gp_names_col %in% c("varid", "rsid"))) {
        stop("Please select gp_in_names as one of varid or rsid, to select how to name dimnames(gp)[[1]]")
    }
    ## here load all samples and variants, unless instructed otherwise
    to.read <- file(bgen_file, "rb")
    ## header
    bgen_header <- load_bgen_header(to.read)
    offset <- bgen_header$offset
    L_H <- bgen_header$L_H
    Layout <- bgen_header$Layout    
    SampleIdentifiers <- bgen_header$SampleIdentifiers
    M <- bgen_header$M
    N <- bgen_header$N
    CompressedSNPBlocks <- bgen_header$CompressedSNPBlocks    
    ## samples
    out <- load_bgen_sample_identifier_block(to.read, L_H + 4, SampleIdentifiers, N)
    sample_names <- out$sample_names
    ## snp info
    out <- load_variant_identifying_data_for_all_snps(to.read, offset, Layout, M, CompressedSNPBlocks)
    var_info <- out$var_info
    per_var_offset <- out$per_var_offset
    per_var_C <- out$per_var_C
    per_var_L_vid <- out$per_var_L_vid
    per_var_L_gdb <- out$per_var_L_gdb
    per_var_num_K_alleles <- out$per_var_num_K_alleles
    ## determine what we need - SNPs
    if (is.null(vars_to_get) == FALSE) {
        which_vars_to_load <- match(vars_to_get, var_info[, "varid"])
        if (sum(is.na(which_vars_to_load)) > 0) {
            stop(paste0(
                "Cannot find some of the variants. The first such missing variants are:",
                paste0(head(vars_to_get[is.na(which_vars_to_load)]), collapse = ", ")
            ))
        }
        ## re-cast some variables
        M <- length(which_vars_to_load)
        var_info <- var_info[which_vars_to_load, , drop = FALSE]
    } else {
        which_vars_to_load <- 1:M
    }
    ## determine what samples we need
    if (is.null(samples_to_get) == FALSE) {
        t <- table(samples_to_get)
        if (sum(t > 1) > 0) {
            stop(paste0("Please do not request the same sample more than once. The following sample was requested ", t[which.max(t > 1)], " times:", names(t)[which.max(t > 1)]))
        }
        which_samples_to_load <- match(samples_to_get, sample_names)
        if (sum(is.na(which_samples_to_load)) > 0) {
            stop(paste0(
                "Cannot find some of the samples. The first such missing samples are:",
                paste0(head(samples_to_get[is.na(which_samples_to_load)]), collapse = ", ")
            ))
        }
        ## now - make extraction vector
        sample_names <- sample_names[which_samples_to_load] ## can re-cast here
        N_to_load <- length(which_samples_to_load)        
    } else {
        required_data_raw_for_samples <- NULL
        which_samples_to_load <- NULL
        N_to_load <- N
    }
    required_data_raw_for_samples <- NULL ## initialize as raw as requires B_bit_prob which we cannot see until we load a SNP    
    ## load genotypes
    gp <- array(NA, c(M, N_to_load, 3))
    for(i_var in 1:length(which_vars_to_load)) {
        var_to_load <- which_vars_to_load[i_var]
        binary_start <- per_var_offset[var_to_load] + per_var_L_vid[var_to_load]
        ## so need to skip C and/or D
        if (CompressedSNPBlocks == 1) {
            binary_start <- binary_start + 8
        } else {
            binary_start <- binary_start + 4
        }
        ## now this is reading genotypes not full genotype data block        
        C <- per_var_C[var_to_load]
        num_K_alleles <- per_var_num_K_alleles[var_to_load]
        out <- load_genotypes_for_one_snp(
            to.read = to.read,
            binary_start = binary_start,
            C = C,
            num_K_alleles = num_K_alleles,
            N = N,
            N_to_load = N_to_load,
            CompressedSNPBlocks = CompressedSNPBlocks,
            required_data_raw_for_samples = required_data_raw_for_samples,
            which_samples_to_load = which_samples_to_load
        )
        gp[i_var, , ] <- out$gen_probs
        if (i_var == 1) {
            ## will be NULL anyway if not requesting all samples
            required_data_raw_for_samples <- out$required_data_raw_for_samples
        }
    }
    ## add names
    dimnames(gp)[[3]] <- c("hom_ref", "het", "hom_alt")    
    dimnames(gp)[[2]] <- sample_names
    dimnames(gp)[[1]] <- var_info[, gp_names_col]    
    close(to.read)
    return(
        list(
            gp = gp,
            sample_names = sample_names,
            var_info = var_info
        )
    )
}


#' @title Load samples names only from a bgen file
#' @param bgen_file Path to bgen file to load
#' @return sample_names
#' @author Robert Davies
#' @export
rrbgen_load_samples <- function(
    bgen_file
) {
    to.read <- file(bgen_file, "rb")
    ## header
    bgen_header <- load_bgen_header(to.read)
    offset <- bgen_header$offset
    L_H <- bgen_header$L_H
    SampleIdentifiers <- bgen_header$SampleIdentifiers    
    N <- bgen_header$N
    ## samples
    out <- load_bgen_sample_identifier_block(to.read, L_H + 4, SampleIdentifiers, N)
    sample_names <- out$sample_names
    close(to.read)
    return(sample_names)
}


#' @title Load variant information only from a bgen file
#' @param bgen_file Path to bgen file to load
#' @return var_info
#' @author Robert Davies
#' @export
rrbgen_load_variant_info <- function(
    bgen_file
) {
    to.read <- file(bgen_file, "rb")
    ## header
    bgen_header <- load_bgen_header(to.read)
    offset <- bgen_header$offset
    Layout <- bgen_header$Layout
    CompressedSNPBlocks <- bgen_header$CompressedSNPBlocks
    M <- bgen_header$M
    ## load all variants
    out <- load_variant_identifying_data_for_all_snps(to.read, offset, Layout, M, CompressedSNPBlocks)
    close(to.read)
    return(out$var_info)
}



## to.read is a file connection
load_bgen_header <- function(to.read) {
    ## header
    seek(to.read, where = 0)
    offset <- readBin(to.read, integer(), endian = "little")
    L_H <- readBin(to.read, integer(), endian = "little")
    M <- readBin(to.read, integer(), endian = "little") ## number of variant data blocks
    N <- readBin(to.read, integer(), endian = "little") ## number of samples
    ## check magic
    magic <- readBin(to.read, size = 1, "raw", n = 4, endian = "little")
    if ((rawToChar(magic) != "bgen") & (rawToChar(magic) != "")) {
        stop(paste0("magic bytes does not equal bgen:", rawToChar(magic)))
    }
    ## I think free is freer than this
    free <- readBin(to.read, integer(),  endian = "little", n = (L_H - 20) / 4)
    ## 
    flag <- readBin(to.read, integer(), endian = "little")
    flag_bits <- intToBits(flag)
    CompressedSNPBlocks <- sum(as.integer(flag_bits[1:2]) * (2 ** (0:1)))
    Layout <- sum(as.integer(flag_bits[3:6]) * (2 ** (0:3)))
    SampleIdentifiers <- as.integer(flag_bits[32])
    return(
        list(
            offset = offset,
            L_H = L_H,
            M = M,
            N = N,
            CompressedSNPBlocks = CompressedSNPBlocks,
            Layout = Layout,
            SampleIdentifiers = SampleIdentifiers
        )
    )
}


load_bgen_sample_identifier_block <- function(to.read, binary_start, SampleIdentifiers, N) {
    ##
    seek(to.read, where = binary_start)
    sample_names <- NULL
    L_SI <- 0
    ## sample identifier block
    if (SampleIdentifiers == 1) {
        L_SI <- readBin(to.read, integer(), endian = "little")
        new_N <- readBin(to.read, integer(), endian = "little")
        if (new_N != N) {
            stop("bad interpretation")
        }
        sample_names <- array(NA, N)
        L_si <- array(NA, N)
        for(i_sample in 1:N) {
            L_si[i_sample] <- readBin(to.read, size = 2, "integer", n = 1, endian = "little")
            name_si <- readBin(to.read, size = 1, "raw", n = L_si[i_sample], endian = "little")
            sample_namei <- rawToChar(name_si)
            sample_names[i_sample] <- sample_namei
        }
        L_SI <- 8 + 2 * N + sum(L_si) ## length of sample identifier block        
    }
    return(
        list(
            sample_names = sample_names,
            L_SI = L_SI
        )
    )
}


load_variant_identifying_data_for_all_snps <- function(to.read, offset, Layout, M, CompressedSNPBlocks) {
    ## first one starts at offset + 4
    ## subsequent ones depend!
    per_var_offset <- array(NA, M)
    per_var_offset[1] <- offset + 4
    per_var_C <- array(NA, M)
    per_var_L_vid <- array(NA, M)
    per_var_L_gdb <- array(NA, M)
    per_var_num_K_alleles <- array(NA, M)
    ## 
    var_info <- array(NA, c(M, 6))
    colnames(var_info) <- c("chr", "varid", "rsid", "position", "ref", "alt")
    ##
    for(i_var in 1:M) {
        out <- load_variant_identifying_data_for_one_snp(to.read, per_var_offset[i_var], Layout, CompressedSNPBlocks)
        per_var_C[i_var] <- out$C
        per_var_L_vid[i_var] <- out$L_vid
        per_var_num_K_alleles[i_var] <- out$num_K_alleles
        per_var_L_gdb[i_var] <- out$L_gdb
        if (i_var < M) {
            per_var_offset[i_var + 1] <- per_var_offset[i_var] + out$L_vid + out$L_gdb
        }
        var_info[i_var, ] <- out$variant_info
    }
    return(
        list(
            var_info = var_info,
            per_var_offset = per_var_offset,
            per_var_C = per_var_C,
            per_var_L_vid = per_var_L_vid,
            per_var_L_gdb = per_var_L_gdb,
            per_var_num_K_alleles = per_var_num_K_alleles
        )
    )
}



## Note I am including C and/or D from the genotype data block in this
load_variant_identifying_data_for_one_snp <- function(to.read, binary_start, Layout, CompressedSNPBlocks) {
    seek(to.read, where = binary_start)
    ## Variant data blocks
    if (Layout == 1) {
        N <- readBin(to.read, integer(), endian = "little")
    }
    ## variant ID
    L_id <- readBin(to.read, size = 2, "integer", n = 1, endian = "little")
    var_id <- readBin(to.read, size = 1, "raw", n = L_id, endian = "little")
    var_id_char <- rawToChar(var_id)
    ## rsid
    L_rsid <- readBin(to.read, size = 2, "integer", n = 1, endian = "little")
    var_rsid <- readBin(to.read, size = 1, "raw", n = L_rsid, endian = "little")
    var_rsid_char <- rawToChar(var_rsid)
    ## chromosome (char)
    L_chr <- readBin(to.read, size = 2, "integer", n = 1, endian = "little")
    var_chr <- readBin(to.read, size = 1, "raw", n = L_chr, endian = "little")
    var_chr_char <- rawToChar(var_chr)
    ## variant position
    var_position <- readBin(to.read, size = 4, "integer", n = 1, endian = "little")
    ## re-produce 
    ## number of K-alleles
    if (Layout == 2) {
        num_K_alleles <- readBin(to.read, size = 2, "integer", n = 1, endian = "little")
    }
    alleles <- array(NA, num_K_alleles)
    L_ai <- array(NA, num_K_alleles)
    for(i_allele in 1:num_K_alleles) {
        L_ai[i_allele] <- readBin(to.read, size = 4, "integer", n = 1, endian = "little")
        ai_char <- rawToChar(readBin(to.read, size = 1, "raw", n = L_ai, endian = "little"))
        alleles[i_allele] <- ai_char
    }
    if (num_K_alleles > 2) {
        stop("num_K_alleles not supported for >2")
    }
    ## make a summary here, similar to .gen file
    variant_info <- c(
        chr = var_chr_char,
        varid = var_id_char,
        rsid = var_rsid_char,
        position = var_position,
        ref = alleles[1],
        alt = alleles[2]
    )
    ## now, also, read next 4 bytes, to get C and D, from the Genotype Data block
    C <- readBin(to.read, size = 4, "integer", n = 1, endian = "little") ## length of rest of the data for this variant
    if (Layout == 2) {    
        if (CompressedSNPBlocks > 0) {
            D <- readBin(to.read, size = 4, "integer", n = 1, endian = "little") ## length after de-compression
        } else {
            D <- C
        }
    }
    ## so L_vid is length of variant identification block
    if (Layout == 1) {
        L_vid <- 16 + 4 * num_K_alleles + L_id + L_rsid + L_chr + sum(L_ai)
    } else if (Layout == 2) {
        L_vid <- 12 + 4 * num_K_alleles + L_id + L_rsid + L_chr + sum(L_ai)
    }
    ## and L_gdb is the length of the genotype data block
    ## together, L_vidand L_gdb are the entire length of the information for that variant 
    L_gdb <- 4 ## for C
    if (CompressedSNPBlocks == 1) {
        ## store D then genotype probability data of length C
        L_gdb <- L_gdb + (C - 4) + 4
    } else if (CompressedSNPBlocks == 0) {
        ## store no D 
        L_gdb <- L_gdb + C
    }
    return(
        list(
            variant_info = variant_info,
            num_K_alleles = num_K_alleles,
            C = C,
            D = D,
            L_vid = L_vid,
            L_gdb = L_gdb
        )
    )
}


## does not include C or D
## this comes from 
load_genotypes_for_one_snp <- function(
    to.read,
    binary_start,
    C,    
    num_K_alleles,
    N,
    N_to_load = NULL,
    CompressedSNPBlocks,
    required_data_raw_for_samples = NULL,
    which_samples_to_load = NULL
) {
    if (is.null(N_to_load)) {
        N_to_load <- N
    }
    seek(to.read, where = binary_start)    
    ## Genotype data block
    data_compressed <- readBin(to.read, size = 1, "raw", n = C - 4, endian = "little")
    if (CompressedSNPBlocks == 1) {
        ## Indicates SNP block probability data is compressed using zlib's compress() function.
        data <- memDecompress(data_compressed, type = "gzip")
        dataC <- rawConnection(data)
    } else if (CompressedSNPBlocks == 0) {
        dataC <- rawConnection(data_compressed)
    }
    new_N <- readBin(dataC, size = 4, "integer", n = 1, endian = "little")
    if (new_N != N) {
        stop("bad interpretation")
    }
    new_num_K_alleles <- readBin(dataC, size = 2, "integer", n = 1, endian = "little")
    if (new_num_K_alleles != num_K_alleles) {
        stop("bad interpretation")
    }
    p_min <- readBin(dataC, size = 1, "integer", n = 1, endian = "little")
    p_max <- readBin(dataC, size = 1, "integer", n = 1, endian = "little")
    N_bytes <- readBin(dataC, size = 1, "raw", n = N, endian = "little")
    ## last bit is missingness
    out <- convert_ploidy_byte(N_bytes) 
    is_missing <- out$is_missing
    ploidy <- out$ploidy
    ##
    phased_flag <- as.integer(readBin(dataC, size = 1, "raw", n = 1, endian = "little"))
    if ((phased_flag != 1) & (phased_flag != 0)) {
        stop("problem with file, phased flag wrong")
    }
    B_bit_prob <- as.integer(readBin(dataC, size = 1, "raw", n = 1, endian = "little"))
    if ((B_bit_prob < 1) | (32 < B_bit_prob)) {
        stop("bit prob is outside range")
    }
    if (phased_flag == 1) {
        ## haplotype probabilities
        stop("haplotype parsing not written")
    } else if (phased_flag == 0) {
        ##
        ## remaining_bytes <- D - 4 - 2 - 1 - 1 - N - 1 - 1
        ## dosage <- array(NA, N)
        if ((B_bit_prob %in% c(8, 16, 24, 32)) == FALSE) {
            stop("non multiple of 8 B_bit_prob not supported")
        }
        if (sum(ploidy[1:N] != 2)) {
            stop("this code does not support non-2 ploidy")
        }
        ## 
        data_raw_for_probs <- readBin(dataC, size = 1, "raw", n = 2 * N * (B_bit_prob / 8), endian = "little", signed = FALSE)
        if (is.null(which_samples_to_load) == FALSE) {
            ## only do this once
            if (is.null(required_data_raw_for_samples)) {
                required_data_raw_for_samples <- unlist(lapply(which_samples_to_load, function(samp) {
                    o <- 2 * B_bit_prob / 8 * (samp - 1)
                    o + 1:(2* B_bit_prob / 8)
                }))
            }
            ##
            data_raw_for_probs <- data_raw_for_probs[required_data_raw_for_samples]
            is_missing <- is_missing[which_samples_to_load]
        }
        ## only convert some of them!
        gen_probs <- rcpp_convert_raw_probabilities_to_double_probabilities(
            data_raw_for_probs,
            N_to_load,
            B_bit_prob,
            is_missing
        )
    }
    close(dataC)
    rm(data)
    return(
        list(
            gen_probs = gen_probs,
            required_data_raw_for_samples = required_data_raw_for_samples
        )
    )
}



convert_ploidy_byte <- function(N_bytes) {
    N <- length(N_bytes)
    is_missing <- array(NA, N)
    ploidy <- array(NA, N)
    ## raw 02 (as.raw(2))and raw 82 (as.raw(130)) are ploidy 2 not missing and yes missing, respectively
    not_missing_ploidy_2 <- N_bytes == as.raw(2)
    is_missing[not_missing_ploidy_2] <- FALSE
    ploidy[not_missing_ploidy_2] <- 2
    ## 
    yes_missing_ploidy_2 <- N_bytes == as.raw(130)
    is_missing[yes_missing_ploidy_2] <- TRUE
    ploidy[yes_missing_ploidy_2] <- 2
    if (sum(is.na(is_missing)) > 0) {
        m <- 2 ** (0:6)
        w <- 1:7
        for(i_sample in which(is_missing)) {
            bits <- rawToBits(N_bytes[i_sample])
            is_missing[i_sample] <- as.logical(bits[8])
            ploidy[i_sample] <- sum(m * as.integer(bits[w]))
        }
    }
    return(
        list(
            is_missing = is_missing,
            ploidy = ploidy
        )
    )
}


convert_raw_probabilities_to_double_probabilities <- function(
    data_raw_for_probs,
    N,
    B_bit_prob,
    is_missing
) {
    gen_probs <- array(NA, c(N, 3))
    B_bit_prob_divide_8 <- B_bit_prob / 8
    w_last <- 1:B_bit_prob_divide_8
    last_byte_used <- 0
    w <- 1:B_bit_prob
    m <- 2 ** (0:(B_bit_prob - 1))
    denom <- (2 ** B_bit_prob - 1)
    for(iSample in 1:N) {
        b_hom_ref <- data_raw_for_probs[last_byte_used + w_last]
        b_het <- data_raw_for_probs[last_byte_used + w_last + B_bit_prob_divide_8]
        last_byte_used <- last_byte_used + B_bit_prob_divide_8 * 2        
        ##x_hom_ref <- data_raw_for_probs[2 * (iSample - 1) + 1]
        ##x_het <- data_raw_for_probs[2 * (iSample - 1) + 2]
        if (is_missing[iSample]) {
            ## dosage[iSample] <- NA
            gen_probs[iSample, 1:3] <- NA
        } else {
            x_hom_ref <- convert_raw_to_int(b_hom_ref, B_bit_prob, w, m)
            ## x_hom_ref <- sum(as.integer(rawToBits(b_hom_ref)[1:B_bit_prob]) * (2 ** (0:(B_bit_prob - 1)))    )
            p_hom_ref <- x_hom_ref / denom
            ## x_het <- sum(as.integer(rawToBits(b_het)[1:B_bit_prob]) * (2 ** (0:(B_bit_prob - 1)))    )
            x_het <- convert_raw_to_int(b_het, B_bit_prob, w, m)
            p_het <- x_het / denom
            p_hom_alt <- 1 - p_hom_ref - p_het
            ## dosage[iSample] <- p_het * 2 * p_hom_alt
            gen_probs[iSample, 1:3] <- c(p_hom_ref, p_het, p_hom_alt)
        }
    }
    return(gen_probs)
}

## R does not seem to do arbitrary byte lengths
## so this *should* work
## just ugly
convert_raw_to_int <- function(x, B_bit_prob, w, m) {
    z <- rawToBits(x)[w]
    return(sum(as.integer(z) * m))
}
