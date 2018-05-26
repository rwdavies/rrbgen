rrbgen_write <- function(
    bgen_file,
    sample_names = NULL,
    var_info = NULL,
    gp = NULL,
    free = NULL,
    Layout = 2,
    CompressedSNPBlocks = 1
) {
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
    ##
    if (is.null(gp) & is.null(var_info) == FALSE) {
        per_var_C <- array(0, M) ## dummy for now
        per_var_D <- array(0, M) ## dummy for now
    }
    ## open connection
    to.write <- file(bgen_file, "wb")
    ## header    
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
    ## samples
    write_bgen_sample_identifier_block(
        to.write,
        L_H,
        sample_names_as_raw,
        L_Si,
        L_SI
    )
    if (is.null(var_info) == FALSE) {
        offset <- L_H + L_SI
        write_variant_identifying_data_for_all_snps(
            to.write = to.write,
            offset = offset,
            var_info = var_info,
            per_var_C = per_var_C,
            per_var_D = per_var_D
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
    CompressedSNPBlocks
) {
    M <- nrow(var_info)
    per_var_snpid_raw <- lapply(var_info[, "snpid"], charToRaw)
    per_var_rsid_raw <- lapply(var_info[, "rsid"], charToRaw)
    per_var_chr_raw <- lapply(var_info[, "chr"], charToRaw)
    per_var_num_K_alleles <- rep(2, M)
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


write_variant_identifying_data_for_all_snps <- function(
    to.write,
    offset,
    var_info,
    per_var_C,
    per_var_D,   
    Layout = 2,
    CompressedSNPBlocks = 1
) {
    M <- as.integer(nrow(var_info))
    out <- prepare_variant_identifying_data_for_all_snps(
        var_info,
        offset,
        Layout,
        CompressedSNPBlocks
    )
    per_var_snpid_raw <- out$per_var_snpid_raw
    per_var_rsid_raw <- out$per_var_rsid_raw
    per_var_chr_raw <- out$per_var_chr_raw
    per_var_num_K_alleles <- out$per_var_num_K_alleles
    per_var_ref_allele_raw <- out$per_var_ref_allele_raw
    per_var_alt_allele_raw <- out$per_var_alt_allele_raw
    per_var_offset <- out$per_var_offset
    ##     
    for(i_var in 1:M) {
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


## if compressed
## do compression (and possibly binary conversion) here
## note this does not include variables C or D
## those are considered part of variant ID block
prepare_genotypes_for_one_snp <- function(gp = NULL, gp_raw = NULL, CompressedSNPBlocks) {
    if (CompressedSNPBlocks == 1) {
        ## Indicates SNP block probability data is compressed using zlib's compress() function.
        data <- memDecompress(data_compressed, type = "gzip")
        dataC <- rawConnection(data)
    } else if (CompressedSNPBlocks == 0) {
        dataC <- rawConnection(data_compressed)
    }
    ##
    ## TODO - DO THIS PROPERLY
    ## WILL INCLUDE per-SNP C and per-SNP D
    ##
    return(
        list(
            gen_probs = gen_probs
        )
    )
}
