#' @title Write variant information to a bgen file
#' @param bgen_file Path to bgen file to be written to
#' @param sample_names character vector of sample names
#' @param var_info var_info (6 first columns of gen file, or chr, varid, rsid, position (1-based), ref, alt)
#' @param gp genotype probabilities as 3-dimensional array, where dim1 = variants, dim2 = samples, dim3 = hom ref, het, hom alt
#' @param list_of_gp_raw_t list of raw data. bespoke format. to understand, investigate the tests
#' @param free free text for header. Should probably remain NULL. Not read if written in rrbgen_read
#' @param Layout Only supported 2 (see bgen spec)
#' @param CompressedSNPBlocks Only supported 1. Whether to (1) or not to (0) compress genotype probabilities
#' @param B_bit_prob How many bits to use to encode genotype probabilities (only supported 8, 16, 24, 32)
#' @author Robert Davies
#' @export
rrbgen_write <- function(
    bgen_file,
    sample_names = NULL,
    var_info = NULL,
    gp = NULL,
    list_of_gp_raw_t = NULL,
    free = NULL,
    Layout = 2,
    CompressedSNPBlocks = 1,
    B_bit_prob = 8,
    close_bgen_file = TRUE,
    add_to_bgen_connection = FALSE,
    bgen_file_connection = NULL,
    previous_offset = NULL,
    header_M = NULL
) {
    if (! (B_bit_prob %in% c(8, 16, 24, 32))) {
        stop("For simplicity, rrbgen can only write probabilities using B bits per entry")
    }
    if (is.null(gp) == FALSE) {
        if (is.null(sample_names) == FALSE) {
            if (length(sample_names) != dim(gp)[[2]]) {
                stop("Incompatible sample names length and gp size")
            }
        }
    }
    ## start with prepratory work
    if (is.null(sample_names) == FALSE) {
        N <- length(sample_names)
        out <- prepare_bgen_sample_identifier_block(sample_names)
        sample_names_as_raw <- out$sample_names_as_raw
        L_Si <- out$L_Si
        L_SI <- out$L_SI
        SampleIdentifiers <- 1        
    } else {
        stop("this has not been written")
    }
    ##
    if (is.null(var_info)) {
        M <- 0
    } else {
        M <- nrow(var_info)
    }
    ## for now, this is fixed
    per_var_num_K_alleles <- rep(2, M)    
    ## 
    if ((is.null(gp) == FALSE) | (is.null(list_of_gp_raw_t) == FALSE)) {
        ## take gps or gp raw, make per_var_C, per_var_D
        out <- prepare_genotypes_for_all_snps(
            gp = gp,
            list_of_gp_raw_t = list_of_gp_raw_t,
            M = M,
            N = N,
            CompressedSNPBlocks = CompressedSNPBlocks,
            B_bit_prob = B_bit_prob,
            per_var_num_K_alleles = per_var_num_K_alleles
        )
        per_var_C <- out$per_var_C
        per_var_D <- out$per_var_D
        data_list <- out$data_list
        dataC_list <- out$dataC_list
    } else {
        ## no genotype probabilities, so nothing to store
        if (CompressedSNPBlocks > 0) {
            per_var_C <- array(4, M) ## to store D only
            per_var_D <- array(0, M) ## uncompressed data
        } else {
            per_var_C <- array(0, M) ## then C and D are same thing
            per_var_D <- array(0, M) ## and are zero
        }
        dataC_list <- NULL
    }
    ## open connection
    if (add_to_bgen_connection == FALSE) {
        to.write <- file(bgen_file, "wb")
    } else {
        to.write <- bgen_file_connection
    }
    ## header - note, could put this lower if prepared header, but need L_H here
    if (add_to_bgen_connection == FALSE) {
        if (is.null(header_M)) {
            ## if are adding to file, then need total number of SNPs up front
            header_M <- M
        }
        L_H <- write_bgen_header(
            to.write,
            L_SI = L_SI,
            header_M = header_M,
            N = N,
            SampleIdentifiers = SampleIdentifiers,
            Layout = Layout,
            CompressedSNPBlocks = CompressedSNPBlocks,
            free = free
        )
        ## 
        offset <- L_H + L_SI
    } else {
        offset <- previous_offset - 4
        ## since +4 is added later
        ## not sure this is the cleanest way to do this, this was surprisingly tricky previously
    }
    ## 
    if (is.null(var_info) == FALSE) {
        var_info_raw_list <- prepare_variant_identifying_data_for_all_snps(
            var_info = var_info,
            offset = offset,
            Layout = Layout,
            CompressedSNPBlocks = CompressedSNPBlocks,
            per_var_num_K_alleles = per_var_num_K_alleles,
            per_var_C = per_var_C
        )
        final_binary_length <- var_info_raw_list$final_binary_length
    } else {
        final_binary_length <- offset + 4
    }
    if (! add_to_bgen_connection) {
        ## write samples
        write_bgen_sample_identifier_block(
            to.write = to.write,
            L_H = L_H,
            sample_names_as_raw = sample_names_as_raw,
            L_Si = L_Si,
            L_SI = L_SI
        )
    }
    if (is.null(var_info) == FALSE) {
        ## note - gp irrelevant at this point, effectively stored in var_info_raw_list
        write_variant_identifying_data_and_genotypes_for_all_snps(
            to.write = to.write,
            offset = offset,
            M = M,
            dataC_list = dataC_list,
            var_info = var_info,
            per_var_C = per_var_C,
            per_var_D = per_var_D,
            per_var_num_K_alleles = per_var_num_K_alleles,
            var_info_raw_list = var_info_raw_list,
            Layout = Layout,
            CompressedSNPBlocks = CompressedSNPBlocks
        )
    }
    ## done writing
    if (close_bgen_file) {
        close(to.write)
        return(NULL)        
    } else {
        ## what do I need?
        return(
            list(
                bgen_file_connection = to.write,
                final_binary_length = final_binary_length
            )
        )
    }
}

write_bgen_header <- function(
    to.write,
    L_SI,
    header_M,
    N,
    SampleIdentifiers = 1,
    Layout = 2,
    CompressedSNPBlocks = 1,
    free = NULL
) {
    ## check flags OK
    if ((CompressedSNPBlocks %in% c(0, 1, 2)) == FALSE) {
        stop(paste0("CompressedSNPBlocks must be either 0, 1 or 2 but you have selected:", CompressedSNPBlocks))
    }
    if ((Layout %in% c(1, 2)) == FALSE) {
        stop(paste0("Layout must be either 1 or 2 but you have selected:", Layout))
    }
    if ((SampleIdentifiers %in% c(0, 1)) == FALSE) {
        stop(paste0("SampleIdentifiers must be either 0 or 1 but you have selected:", SampleIdentifiers))
    }
    ## probably not necessary as this is under my control, but just in case
    if ((header_M < 0) | (round(header_M) != header_M)) {
        stop(paste0("The number of SNPs M must be an integer but you have selected:", header_M))
    }
    if ((N < 0) | (round(N) != N)) {
        stop(paste0("The number of samples N must be an integer but you have selected:", N))
    }
    ## 
    if (is.null(free)) {
        L_H <- 20
    } else {
        ## not 100% sure this work
        if (class(free) != "raw") {
            stop("free must be of class raw")
        }
        L_H <- length(raw) + 20
    }
    ## writing part now
    seek(to.write, where = 0)
    offset <- L_H + L_SI    
    writeBin(as.integer(offset), to.write, endian = "little")
    writeBin(as.integer(L_H), to.write, endian = "little")    
    writeBin(as.integer(header_M), to.write, endian = "little")
    writeBin(as.integer(N), to.write, endian = "little")
    magic_raw <- as.raw(sapply(c("b", "g", "e", "n"), function(x) charToRaw(x)))    
    writeBin(magic_raw, to.write, size = 1, endian = "little")
    if (L_H > 20) {
        ## Not sure this will work?
        writeBin(free, to.write, endian = "little")
    }
    ## flag bits
    flag_bits <- as.raw(rep(0, 32))
    ## 
    flag_bits[1:2] <- intToBits(CompressedSNPBlocks)[1:2]
    flag_bits[3:6] <- intToBits(Layout)[1:4]
    flag_bits[32] <- as.raw(SampleIdentifiers)
    writeBin(packBits(flag_bits), to.write, endian = "little")
    return(L_H)
}


prepare_bgen_sample_identifier_block <- function(sample_names, SampleIdentifiers = 1) {
    if (SampleIdentifiers == 1) {
        sample_names_as_raw <- lapply(sample_names, charToRaw)
        L_Si <- lapply(sample_names_as_raw, length) ## per-sample
        L_SI <- as.integer(8 + 2 * length(sample_names) + sum(unlist(L_Si))) ## total
    } else {
        sample_names_as_raw <- NULL
        L_Si <- NULL
        L_SI <- 0
    }
    return(
        list(
            sample_names_as_raw = sample_names_as_raw,
            L_Si = L_Si,
            L_SI = L_SI
        )
    )
}

write_bgen_sample_identifier_block <- function(
    to.write,
    L_H,
    sample_names_as_raw,
    L_Si,
    L_SI,
    SampleIdentifiers = 1
) {
    ##
    seek(to.write, where = L_H + 4)
    if (SampleIdentifiers == 1) {
        ## first, convert and build rest
        N <- length(sample_names_as_raw)        
        ## can now write
        writeBin(L_SI, to.write, endian = "little")
        writeBin(as.integer(N), to.write, endian = "little")
        for(iSample in 1:N) {
            writeBin(L_Si[[iSample]], to.write, size = 2, endian = "little")
            writeBin(sample_names_as_raw[[iSample]], to.write, size = 1, endian = "little")
        }
    }
    return(NULL)
}



prepare_variant_identifying_data_for_all_snps <- function(
    var_info,
    offset,
    Layout,
    CompressedSNPBlocks,
    per_var_num_K_alleles,
    per_var_C
) {
    M <- nrow(var_info)
    per_var_varid_raw <- lapply(var_info[, "varid"], charToRaw)
    per_var_rsid_raw <- lapply(var_info[, "rsid"], charToRaw)
    per_var_chr_raw <- lapply(var_info[, "chr"], charToRaw)
    per_var_ref_allele_raw <- lapply(var_info[, "ref"], charToRaw)
    per_var_alt_allele_raw <- lapply(var_info[, "alt"], charToRaw)
    ## now calculate length of each variant identifying block including C (and D)
    per_var_L_vid <- array(0, M)
    ## Variant data blocks
    if (Layout == 1) {
        per_var_L_vid <- per_var_L_vid + 4 ## first optional N
    }
    per_var_L_vid <- per_var_L_vid + 16 +
        4 * sapply(per_var_num_K_alleles, length) +
        sapply(per_var_varid_raw, length) +
        sapply(per_var_rsid_raw, length) +
        sapply(per_var_chr_raw, length) + 
        sapply(per_var_ref_allele_raw, length) + 
        sapply(per_var_alt_allele_raw, length)
    ## separately, store length of genotype data block
    per_var_L_gdb <- rep(4, M) ## for C
    if (CompressedSNPBlocks == 1) {
        ## store D then genotype probability data of length C
        per_var_L_gdb <- per_var_L_gdb + (per_var_C - 4) + 4
    } else if (CompressedSNPBlocks == 0) {
        ## store no D 
        per_var_L_gdb <- per_var_L_gdb + per_var_C
    }
    ## 
    per_var_offset <- array(NA, M)
    per_var_offset[1] <- offset + 4
    if (M > 1) {
        for(i_var in 2:M) {
            per_var_offset[i_var] <- per_var_offset[i_var - 1] + per_var_L_vid[i_var - 1] + per_var_L_gdb[i_var - 1]
        }
    }
    final_binary_length <- per_var_offset[M] + per_var_L_vid[M] + per_var_L_gdb[M]
    return(
        list(
            per_var_varid_raw = per_var_varid_raw,
            per_var_rsid_raw = per_var_rsid_raw,
            per_var_chr_raw = per_var_chr_raw,
            per_var_num_K_alleles = per_var_num_K_alleles,
            per_var_ref_allele_raw = per_var_ref_allele_raw,
            per_var_alt_allele_raw = per_var_alt_allele_raw,
            per_var_offset = per_var_offset,
            per_var_L_vid = per_var_L_vid,
            per_var_L_gdb = per_var_L_gdb,
            final_binary_length = final_binary_length
        )
    )
}


write_variant_identifying_data_and_genotypes_for_all_snps <- function(
    to.write,
    offset,
    M,
    dataC_list,
    var_info_raw_list,
    var_info,
    per_var_C,
    per_var_D,
    per_var_num_K_alleles,
    Layout = 2,
    CompressedSNPBlocks = 1
) {
    per_var_varid_raw <- var_info_raw_list$per_var_varid_raw
    per_var_rsid_raw <- var_info_raw_list$per_var_rsid_raw
    per_var_chr_raw <- var_info_raw_list$per_var_chr_raw
    per_var_ref_allele_raw <- var_info_raw_list$per_var_ref_allele_raw
    per_var_alt_allele_raw <- var_info_raw_list$per_var_alt_allele_raw
    per_var_offset <- var_info_raw_list$per_var_offset
    per_var_L_vid <- var_info_raw_list$per_var_L_vid
    per_var_L_gdb <- var_info_raw_list$per_var_L_gdb
    for(i_var in 1:M) {
        ## Variant identifying data (Plus C and D)
        write_variant_identifying_data_for_one_snp(
            to.write,
            binary_start = per_var_offset[i_var],
            C = per_var_C[i_var],
            D = per_var_D[i_var],
            varid_raw = per_var_varid_raw[[i_var]],
            rsid_raw = per_var_rsid_raw[[i_var]],
            chr_raw = per_var_chr_raw[[i_var]],
            position = var_info[i_var, "position"],
            num_K_alleles = per_var_num_K_alleles[i_var],
            ref_allele_raw = per_var_ref_allele_raw[[i_var]],
            alt_allele_raw = per_var_alt_allele_raw[[i_var]],
            Layout = Layout,
            CompressedSNPBlocks = CompressedSNPBlocks
        )
        ## so assume compressed bit is length C not C-4 
        ## Genotype block (Layout 2) (Minus C and D)
        if (is.null(dataC_list) == FALSE) {
            where <- per_var_offset[i_var] + per_var_L_vid[i_var] + 4
            if (CompressedSNPBlocks == 1) {
                where <- where + 4
            }
            seek(to.write, where = where)
            writeBin(dataC_list[[i_var]], to.write, size = 1, endian = "little")
        }
    }
    return(
        list(
            per_var_offset = per_var_offset
        )
    )
}




## Note I am including C and/or D from the genotype data block in this
write_variant_identifying_data_for_one_snp <- function(
    to.write,
    binary_start,
    C,
    D,
    varid_raw,
    rsid_raw,
    chr_raw,
    position,    
    num_K_alleles,
    ref_allele_raw,
    alt_allele_raw,
    Layout,
    CompressedSNPBlocks
) {
    seek(to.write, where = binary_start)
    if (Layout == 1) {
        writeBin(as.integer(N), to.write, endian = "little")
    }
    for(to_write_raw in list(varid_raw, rsid_raw, chr_raw)) {
        writeBin(as.integer(length(to_write_raw)), to.write, size = 2, endian = "little")
        writeBin(to_write_raw, to.write, size = 1, endian = "little")
    }
    ## variant position
    writeBin(as.integer(position), to.write, size = 4, endian = "little")
    ## number K of alleles - for now only support bi-allelic
    if (num_K_alleles != 2) {
        stop("functionality not written for >2 alleles")
    }
    if (Layout == 2) {
        writeBin(as.integer(num_K_alleles), to.write, size = 2, endian = "little")
    }
    ## write alleles and their length
    for(allele_raw in c(ref_allele_raw, alt_allele_raw)) {
        writeBin(as.integer(length(allele_raw)), to.write, size = 4, endian = "little")
        writeBin(allele_raw, to.write, size = 1, endian = "little")
    }
    ## now, also, write C and D, from the Genotype Data block
    writeBin(as.integer(C), to.write, size = 4, endian = "little")
    ##
    if (Layout == 2) {
        if (CompressedSNPBlocks > 0) {
            writeBin(as.integer(D), to.write, size = 4, endian = "little")
        }
    }
    return(NULL)
}



prepare_genotypes_for_all_snps <- function(
    gp = NULL,
    list_of_gp_raw_t = NULL,
    M,
    N,
    CompressedSNPBlocks,
    B_bit_prob,
    per_var_num_K_alleles
) {
    if ((as.integer(is.null(gp)) + as.integer(is.null(list_of_gp_raw_t))) != 1) {
        stop("Please write using either one of gp or list_of_gp_raw_t")
    }
    ## build container?
    per_var_C <- array(0, M)
    per_var_D <- array(0, M)
    data_list <- as.list(1:M)
    dataC_list <- as.list(1:M)    
    for(i_snp in 1:M) {
        out <- prepare_genotypes_for_one_snp(
            gp = gp,
            list_of_gp_raw_t = list_of_gp_raw_t,
            i_snp = i_snp,
            CompressedSNPBlocks = CompressedSNPBlocks,
            B_bit_prob = B_bit_prob,
            N = N,
            num_K_alleles = per_var_num_K_alleles[i_snp]
        )
        per_var_C[i_snp] <- out$C
        per_var_D[i_snp] <- out$D
        data_list[[i_snp]] <- out$data
        dataC_list[[i_snp]] <- out$dataC
    }
    return(
        list(
            per_var_C = per_var_C,
            per_var_D = per_var_D,            
            data_list = data_list,
            dataC_list = dataC_list
        )
    )
}

## if compressed
## do compression (and possibly binary conversion) here
## note this does not include variables C or D
## those are considered part of variant ID block
prepare_genotypes_for_one_snp <- function(
    gp = NULL,
    list_of_gp_raw_t = NULL,
    i_snp,
    CompressedSNPBlocks,
    B_bit_prob,
    N,
    num_K_alleles
) {
    ## make whole store here
    v <- vector(mode = "raw")            
    non_gp_stuff_length <- 4 + 2 + 1 + 1 + N + 1 + 1
    data <- vector(mode = "raw", length = non_gp_stuff_length + 2 * N * (B_bit_prob / 8))
    ## do pre-genotype probability things
    data[1:4] <- writeBin(as.integer(N), v, size = 4, endian = "little")
    ## number_of_alleles
    data[5:6] <- writeBin(as.integer(num_K_alleles), v, size = 2, endian = "little")
    ## min then max ploidy
    data[7] <- writeBin(as.integer(2), v, size = 1, endian = "little") ## P_min
    data[8] <- writeBin(as.integer(2), v, size = 1, endian = "little") ## P_max
    ## missing and ploidy
    local_offset <- 4 + 2 + 1 + 1
    ##
    ## two options here
    ## one - using genotype probabilities
    ## other, using weird STITCH format
    ## 
    if (is.null(gp) == FALSE) {
        missing <- is.na(gp[i_snp, , 1])
        data[local_offset + 1:N] <- make_ploidy_raw(missing)
        ## finally, last 2 things
        phased_flag <- 0
        data[4 + 2 + 1 + 1 + N + 1] <- writeBin(as.integer(phased_flag), v, size = 1, endian = "little")
        data[4 + 2 + 1 + 1 + N + 1 + 1] <- writeBin(as.integer(B_bit_prob), v, size = 1, endian = "little")
        ## genotype probabilities
        data[non_gp_stuff_length + 1:(2 * N * (B_bit_prob / 8))] <- rcpp_make_raw_data_vector_for_probabilities(
            gp = gp,
            B_bit_prob = B_bit_prob,
            i_snp_1_based = i_snp
        )
    } else {
        data[local_offset + 1:N] <- as.raw(2) ## 2 ploidy, no missing
        phased_flag <- 0
        data[4 + 2 + 1 + 1 + N + 1] <- writeBin(as.integer(phased_flag), v, size = 1, endian = "little")
        data[4 + 2 + 1 + 1 + N + 1 + 1] <- writeBin(as.integer(B_bit_prob), v, size = 1, endian = "little")
        ## genotpe probabilities - already in list form!
        c <- non_gp_stuff_length
        for(i_sheet in 1:length(list_of_gp_raw_t)) {
            X <- nrow(list_of_gp_raw_t[[i_sheet]])
            data[c + 1:X] <- list_of_gp_raw_t[[i_sheet]][, i_snp]
            c <- c + X
        }
    }
    ## 
    if (CompressedSNPBlocks == 1) {
        dataC <- memCompress(data, type = "gzip")        
    } else {
        dataC <- data
    }
    C <- length(dataC) + 4
    D <- length(data)
    return(
        list(
            C = C,
            D = D,
            dataC = dataC,
            data = data
        )
    )
}


make_ploidy_raw <- function(missing) {
    N <- length(missing)
    ## there are only two values under current assumptions
    ## make this naive but fast for now
    not_missing_ploidy <- rawToBits(as.raw(as.integer(2))) ## assume everyone 2 for now
    missing_ploidy <- not_missing_ploidy
    missing_ploidy[8] <- as.raw(1)
    not_missing_ploidy <- packBits(not_missing_ploidy)    
    missing_ploidy <- packBits(missing_ploidy)
    data_local <- vector("raw", N)
    data_local[missing == TRUE] <- missing_ploidy
    data_local[missing == FALSE] <- not_missing_ploidy
    return(data_local)
    ##w <- 1:8
    ##for(i_sample in 1:N) {
    ##    ## do not apply check
    ##    ploidy <- default_ploidy
    ##    if (missing[i_sample]) {
    ##        ploidy[8] <- as.raw(1)
    ##    }
    ##    ## last byte is missing
    ##    data_local[8 * (i_sample - 1) + w] <- ploidy
    ## }
    ## v <- vector(mode = "raw", N)        
    ##return(writeBin(packBits(data_local), v, size = 1, endian = "little"))
}



## prob_mat has rows = samples
## columns are hom_ref, het, hom_alt
make_raw_data_vector_for_probabilities <- function(
   gp,
   B_bit_prob,
   i_snp_1_based
) {
    N <- dim(gp)[[2]]
    data_local <- vector(mode = "raw", length = N * (B_bit_prob / 8))
    ## constants
    B_bit_prob_divide_8 <- (B_bit_prob / 8)
    const_2_bit_minus_1 <- (2 ** B_bit_prob - 1)
    const_where <- list(
        1:8,
        9:16,
        17:24,
        25:32
    )
    last_byte_used <- 0
    ## be careful - no seek used below 
    ## do not change for loops without checking format specifications
    for(i_sample in 1:N) {
        ## get int representatiosn using rule
        if (is.na(gp[i_snp_1_based, i_sample, 1]))  {
            v <- rep(0, 3)
        } else {
            v <- gp[i_snp_1_based, i_sample, ] * (const_2_bit_minus_1)
            v2 <- v - floor(v)
            F <- as.integer(round(sum(v2)))
            if (F > 0) {
                o <- order(v2, decreasing = TRUE)            
                for(i in 1:F) {
                    v[o[i]] <- ceiling(v[o[i]])
                }
                for(i in F:3) {
                    v[o[i]] <- floor(v[o[i]])
                }
            }
        }
        for(i_gen in 1:2) {
            x_hom_ref <- v[i_gen]
            x_bits <- robbie_intTobits(x_hom_ref, B_bit_prob) 
            for(i_B in 1:(B_bit_prob_divide_8)) {
                data_local[last_byte_used + i_B] <-
                    packBits(x_bits[const_where[[i_B]]], type = "raw")
            }
            last_byte_used <- last_byte_used + B_bit_prob_divide_8            
        }
    }
    return(data_local)
}



## want unsigned integer
## so B_bit_prob = 8, can do 0 through 
robbie_intTobits <- function(x, B_bit_prob) {
    if (is.na(x)) {
        return(vector("raw", B_bit_prob))
    }
    if (B_bit_prob == 32) {
        ## 2 ** 32 
        if (x >= (2147483648)) {
            to_out <- intToBits(x - 2147483648)
            to_out[32] <- as.raw(1)
            return(to_out)
        } else {
            return(intToBits(x))
        }
    } else {
        return(intToBits(x))
    }
}

## if I ever want to do this cleaner...
## 
##    if (class(x) != "integer") {
##        stop("x must be an integer")
##    }
##    if (x < 0 | x > 2 ** B_bit_prob) {
##        stop("x out of range")
##    }
##    out <- vector(mode = "raw", B_bit_prob)
##    if (x == 0) {
##        return(out)
##    }
##    for(iB in (B_bit_prob - 1):1) {
##        2 ** iB## was here
##    }
## }
