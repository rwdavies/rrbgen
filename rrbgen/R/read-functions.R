## here is the general purpose API
## later, instruct only loading certain samples or snps
## have separate samples and SNPs below
rrbgen_load <- function(
    bgen_file,
    gp_names_col = "snpid"
) {
    if (! (gp_names_col %in% c("snpid", "rsid"))) {
        stop("Please select gp_in_names as one of snpid or rsid, to select how to name dimnames(gp)[[1]]")
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
    snp_info <- out$snp_info
    per_var_offset <- out$per_var_offset
    per_var_C <- out$per_var_C
    per_var_L_vid <- out$per_var_L_vid
    per_var_num_K_alleles <- out$per_var_num_K_alleles
    ## load genotypes
    gp <- array(NA, c(M, N, 3))
    for(i_var in 1:M) {
        binary_start <- per_var_offset[i_var] + per_var_L_vid[i_var]
        C <- per_var_C[i_var]
        num_K_alleles <- per_var_num_K_alleles[i_var]
        out <- load_genotypes_for_one_snp(to.read, binary_start, num_K_alleles, N, CompressedSNPBlocks, C)
        gp[i_var, , ] <- out$gen_probs
    }
    ## add names
    dimnames(gp)[[1]] <- snp_info[, gp_names_col]
    dimnames(gp)[[2]] <- sample_names
    dimnames(gp)[[3]] <- c("hom_ref", "het", "hom_alt")
    return(
        list(
            gp = gp,
            sample_names = sample_names
        )
    )
}


## just load the sample names
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


## just load the variant information
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
    return(out$snp_info)
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
    per_var_num_K_alleles <- array(NA, M)
    ## 
    snp_info <- array(NA, c(M, 6))
    colnames(snp_info) <- c("chr", "snpid", "rsid", "position", "ref", "alt")
    ##     
    for(i_var in 1:M) {
        out <- load_variant_identifying_data_for_one_snp(to.read, per_var_offset[i_var], Layout, CompressedSNPBlocks)
        per_var_C[i_var] <- out$C
        per_var_L_vid[i_var] <- out$L_vid
        per_var_num_K_alleles[i_var] <- out$num_K_alleles
        if (i_var < M) {
            per_var_offset[i_var + 1] <- per_var_offset[i_var] + out$L_vid
            if (CompressedSNPBlocks == 1) {
                per_var_offset[i_var + 1] <- per_var_offset[i_var + 1] + out$C - 4
            } else if (CompressedSNPBlocks == 0) {
                per_var_offset[i_var + 1] <- per_var_offset[i_var + 1] + out$C
            }
        }
        snp_info[i_var, ] <- out$variant_info
    }
    return(
        list(
            snp_info = snp_info,
            per_var_offset = per_var_offset,
            per_var_C = per_var_C,
            per_var_L_vid = per_var_L_vid,
            per_var_num_K_alleles = per_var_num_K_alleles
        )
    )
}



## Note I am including C and/or D from the genotype data block in this
load_variant_identifying_data_for_one_snp <- function(to.read, binary_start, Layout, CompressedSNPBlocks) {
    seek(to.read, where = binary_start)
    L_vid <- 0
    ## Variant data blocks
    if (Layout == 1) {
        N <- readBin(to.read, integer(), endian = "little")
        L_vid <- L_vid + 4
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
        snpid = var_id_char,
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
            L_vid <- L_vid + 4
        } else {
            D <- C
        }
    }
    ## get the total number of bytes this used. note this now includes C and D
    L_vid <- L_vid + 16 + 4 * num_K_alleles + L_id + L_rsid + L_chr + sum(L_ai)
    ## I think this is an error - this is captured above
    ## if (Layout == 1) {
    ##     L_vid <- L_vid + 4
    ## }
    return(
        list(
            variant_info = variant_info,
            num_K_alleles = num_K_alleles,
            C = C,
            D = D,
            L_vid = L_vid
        )
    )
}


## does not include C or D
## this comes from 
load_genotypes_for_one_snp <- function(to.read, binary_start, num_K_alleles, N, CompressedSNPBlocks, C) {
    seek(to.read, where = binary_start)    
    ## Genotype data block
    data_compressed <- readBin(to.read, size = 1, "raw", n = C, endian = "little")
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
    is_missing <- as.logical(sapply(N_bytes, function(byte) rawToBits(byte)[8]))
    ploidy <- sapply(N_bytes, function(byte) sum(2 ** (0:6) * as.integer(rawToBits(byte)[1:7])))
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
        gen_probs <- array(NA, c(N, 3))
        if ((B_bit_prob %in% c(8, 16, 24, 32)) == FALSE) {
        stop("non multiple of 8 B_bit_prob not supported")
        }
        for(iSample in 1:N) {
            if (ploidy[iSample] != 2) {
                stop("this code does not support non-2 ploidy")
            }
            b_hom_ref <- readBin(dataC, size = 1, "raw", n = B_bit_prob / 8, endian = "little")
            b_het <- readBin(dataC, size = 1, "raw", n = B_bit_prob / 8, endian = "little")        
            if (is_missing[iSample]) {
                ## dosage[iSample] <- NA
                gen_probs[iSample, 1:3] <- NA
            } else {
                x <- sum(as.integer(rawToBits(b_hom_ref)) * (2 ** (0:(B_bit_prob - 1)))    )
                p_hom_ref <- x / (2 ** B_bit_prob)
                x <- sum(as.integer(rawToBits(b_het)) * (2 ** (0:(B_bit_prob - 1)))    )
                p_het <- x / (2 ** B_bit_prob)
                p_hom_alt <- 1 - p_hom_ref - p_het
                ## dosage[iSample] <- p_het * 2 * p_hom_alt
                gen_probs[iSample, 1:3] <- c(p_hom_ref, p_het, p_hom_alt)
            }
        }
    }
    close(dataC)
    rm(data)
    return(
        list(
            gen_probs = gen_probs
        )
    )
}
