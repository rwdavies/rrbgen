rrbgen_write <- function(
    bgen_file,
    sample_names = NULL,
    var_info = NULL,
    gp = NULL,
    gp_raw = NULL,
    free = NULL,
    Layout = 2,
    CompressedSNPBlocks = 1,
    B_bit_prob = 8
) {
    if (! (B_bit_prob %in% c(8, 16, 24, 32))) {
        stop("For simplicity, rrbgen can only write probabilities using B bits per entry")
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
    if (is.null(gp) == FALSE) {
        ## take gps or gp raw, make per_var_C, per_var_D
        out <- prepare_genotypes_for_all_snps(
            gp = gp,
            gp_raw = NULL,
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
        per_var_C <- array(0, M)
        per_var_D <- array(0, M)
    }
    ## open connection
    to.write <- file(bgen_file, "wb")
    ## header - note, could put this lower if prepared header, but need L_H here
    L_H <- write_bgen_header(
        to.write,
        L_SI = L_SI,
        M = M,
        N = N,
        SampleIdentifiers = SampleIdentifiers,
        Layout = Layout,
        CompressedSNPBlocks = CompressedSNPBlocks,
        free = free
    )
    offset <- L_H + L_SI    
    ##     
    if (is.null(var_info) == FALSE) {
        M <- as.integer(nrow(var_info))
        var_info_raw_list <- prepare_variant_identifying_data_for_all_snps(
            var_info = var_info,
            offset = offset,
            Layout = Layout,
            CompressedSNPBlocks = CompressedSNPBlocks,
            per_var_num_K_alleles = per_var_num_K_alleles
        )
    }
    ## samples
    write_bgen_sample_identifier_block(
        to.write = to.write,
        L_H = L_H,
        sample_names_as_raw = sample_names_as_raw,
        L_Si = L_Si,
        L_SI = L_SI
    )
    if (is.null(var_info) == FALSE) {
        ## note - gp irrelevant at this point, effectively stored in var_info_raw_list
        write_variant_identifying_data_and_genotypes_for_all_snps(
            to.write = to.write,
            offset = offset,
            M = M,
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
    close(to.write)
    return(NULL)
}

write_bgen_header <- function(
    to.write,
    L_SI,
    M,
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
    if ((M < 0) | (round(M) != M)) {
        stop(paste0("The number of SNPs M must be an integer but you have selected:", M))
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
    writeBin(as.integer(M), to.write, endian = "little")
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
    per_var_num_K_alleles
) {
    M <- nrow(var_info)
    per_var_snpid_raw <- lapply(var_info[, "snpid"], charToRaw)
    per_var_rsid_raw <- lapply(var_info[, "rsid"], charToRaw)
    per_var_chr_raw <- lapply(var_info[, "chr"], charToRaw)
    per_var_ref_allele_raw <- lapply(var_info[, "ref"], charToRaw)
    per_var_alt_allele_raw <- lapply(var_info[, "alt"], charToRaw)
    ## now calculate length of each variant identifying block including C (and D)
    L_vid <- array(0, M)
    ## Variant data blocks
    if (Layout == 1) {
        L_vid <- L_vid + 4 ## first optional N
    }
    if (Layout == 2) {
        if (CompressedSNPBlocks > 0) {
            L_vid <- L_vid + 4 ## to store D
        }
    }
    L_vid <- L_vid + 16 +
    4 * sapply(per_var_num_K_alleles, length) +
        sapply(per_var_snpid_raw, length) +
        sapply(per_var_rsid_raw, length) +
        sapply(per_var_chr_raw, length) + 
        sapply(per_var_ref_allele_raw, length) + 
        sapply(per_var_alt_allele_raw, length)
    ##
    per_var_offset <- array(NA, M)
    per_var_offset[1] <- offset + 4
    for(i_var in 2:M) {
        per_var_offset[i_var] <- per_var_offset[i_var - 1] + L_vid[i_var - 1]
    }
    return(
        list(
            per_var_snpid_raw = per_var_snpid_raw,
            per_var_rsid_raw = per_var_rsid_raw,
            per_var_chr_raw = per_var_chr_raw,
            per_var_num_K_alleles = per_var_num_K_alleles,
            per_var_ref_allele_raw = per_var_ref_allele_raw,
            per_var_alt_allele_raw = per_var_alt_allele_raw,
            per_var_offset = per_var_offset
        )
    )
}


write_variant_identifying_data_and_genotypes_for_all_snps <- function(
    to.write,
    offset,
    M,
    var_info_raw_list,
    var_info,
    per_var_C,
    per_var_D,
    per_var_num_K_alleles,
    Layout = 2,
    CompressedSNPBlocks = 1
) {
    per_var_snpid_raw <- var_info_raw_list$per_var_snpid_raw
    per_var_rsid_raw <- var_info_raw_list$per_var_rsid_raw
    per_var_chr_raw <- var_info_raw_list$per_var_chr_raw
    per_var_ref_allele_raw <- var_info_raw_list$per_var_ref_allele_raw
    per_var_alt_allele_raw <- var_info_raw_list$per_var_alt_allele_raw
    per_var_offset <- var_info_raw_list$per_var_offset
    for(i_var in 1:M) {
        ## Variant identifying data (Plus C and D)
        write_variant_identifying_data_for_one_snp(
            to.write,
            binary_start = per_var_offset[i_var],
            C = per_var_C[i_var],
            D = per_var_D[i_var],
            snpid_raw = per_var_snpid_raw[[i_var]],
            rsid_raw = per_var_rsid_raw[[i_var]],
            chr_raw = per_var_chr_raw[[i_var]],
            position = var_info[i_var, "position"],
            num_K_alleles = per_var_num_K_alleles[i_var],
            ref_allele_raw = per_var_ref_allele_raw[[i_var]],
            alt_allele_raw = per_var_alt_allele_raw[[i_var]],
            Layout = Layout,
            CompressedSNPBlocks = CompressedSNPBlocks
        )
        ## Genotype block (Layout 2) (Minus C and D)
        ## TODO 
        ##write_genotype_probabilities_for_one_snp(
        ## 
        ## )
    }
    return(NULL)
}




## Note I am including C and/or D from the genotype data block in this
write_variant_identifying_data_for_one_snp <- function(
    to.write,
    binary_start,
    C,
    D,
    snpid_raw,
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
    for(to_write_raw in list(snpid_raw, rsid_raw, chr_raw)) {
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
    if (Layout == 2) {
        writeBin(as.integer(D), to.write, size = 4, endian = "little")        
    }
    return(NULL)
}



prepare_genotypes_for_all_snps <- function(
    gp = NULL,
    gp_raw = NULL,
    M,
    N,
    CompressedSNPBlocks,
    B_bit_prob,
    per_var_num_K_alleles
) {
    if ((as.integer(is.null(gp)) + as.integer(is.null(gp_raw))) != 1) {
        stop("Please write using either one of gp or gp_raw")
    }
    ## build container?
    per_var_C <- array(0, M)
    per_var_D <- array(0, M)
    data_list <- as.list(1:M)
    dataC_list <- as.list(1:M)    
    for(i_snp in 1:M) {
        out <- prepare_genotypes_for_one_snp(
            gp = gp,
            gp_raw = gp_raw,
            i_snp = i_snp,
            CompressedSNPBlocks = CompressedSNPBlocks,
            B_bit_prob = B_bit_prob,
            N = N,
            num_K_alleles = per_var_num_K_alleles[i_snp]
        )
        per_var_C[i_snp] <- out$C
        per_var_D[i_snp] <- out$D
        data_list <- out$data
        dataC_list <- out$dataC
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
    gp_raw = NULL,
    i_snp,
    CompressedSNPBlocks,
    B_bit_prob,
    N,
    num_K_alleles
) {
    ## make whole store here
    non_gp_stuff_length <- 4 + 2 + 1 + 1 + N + 1 + 1
    data <- vector(mode = "raw", length = non_gp_stuff_length + N * (B_bit_prob / 8))
    v <- vector(mode = "raw")
    ## do pre-genotype probability things
    data[1:4] <- writeBin(as.integer(N), v, size = 4, endian = "little")
    ## number_of_alleles
    data[5:6] <- writeBin(as.integer(num_K_alleles), v, size = 2, endian = "little")
    ## min then max ploidy
    data[7] <- writeBin(as.integer(2), v, size = 1, endian = "little") ## P_min
    data[8] <- writeBin(as.integer(2), v, size = 1, endian = "little") ## P_max
    ## missing and ploidy
    local_offset <- 4 + 2 + 1 + 1
    default_ploidy <- rawToBits(as.raw(as.integer(2))) ## assume everyone 2 for now
    ## 
    for(i_sample in 1:N) {
        ## do not apply check
        x <- gp[i_snp, i_sample, ]
        s <- sum(is.na(x))
        ploidy <- default_ploidy
        ploidy[8] <- as.raw(s)
        ## last byte is missing
        data[local_offset + i_sample] <- writeBin(packBits(ploidy), v, size = 1, endian = "little")
    }
    ## finally, last 2 things
    phased_flag <- 0
    data[4 + 2 + 1 + 1 + N + 1] <- writeBin(as.integer(phased_flag), v, size = 1, endian = "little")
    data[4 + 2 + 1 + 1 + N + 1 + 1] <- writeBin(as.integer(B_bit_prob), v, size = 1, endian = "little")
    ## now, work on genotype probabilities
    B_bit_prob_divide_8 <- (B_bit_prob / 8)
    const_2_bit <- (2 ** B_bit_prob)
    const_where <- list(
        1:8,
        9:16,
        17:24,
        25:32
    )
    last_byte_used <- non_gp_stuff_length
    ## be careful - no seek used below 
    ## do not change for loops without checking format specifications
    for(i_sample in 1:N) {
        ## hom_ref then alt
        for(i_gen in 1:2) {
            x_hom_ref <- round(gp[i_snp, i_sample, c("hom_ref", "het")[i_gen]] * const_2_bit)
            x_bits <- intToBits(x_hom_ref)
            for(i_B in (B_bit_prob / 8)) {
                ## (2 * B_bit_prob_divide_8) * (i_sample - 1) + 2 * (i_B - 1) + i_gen
                data[last_byte_used + 1:B_bit_prob_divide_8] <-
                    packBits(x_bits[const_where[[i_B]]], type = "raw")
                last_byte_used <- last_byte_used + B_bit_prob_divide_8
            }
        }
    }
    ## 
    if (CompressedSNPBlocks == 1) {
        dataC <- memCompress(data, type = "gzip")        
    } else {
        dataC <- data
    }
    C <- length(dataC)
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


write_genotype_probabilities_for_one_snp <- function(
    to.write,
    binary_start,
    CompressedSNPBlocks
) {
    seek(to.write, where = binary_start)
    ## write stuff here
    return(NULL)
}

