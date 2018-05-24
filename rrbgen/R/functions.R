## here is the general purpose API
## have separate samples and SNPs below
rrbgen_load <- function(
    file,
    samples = NULL,
    snps = NULL
) {
    print("TODO")
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


## just load the SNP information
rrbgen_load_snps <- function(
    file
) {

    close(to.read)
    to.read <- file(bgen_file, "rb")
    bgen_header <- load_bgen_header(to.read)
    
    L_H <- bgen_header$L_H ## length of header in bytes
    Layout <- bgen_header$Layout
    CompressedSNPBlocks <- bgen_header$CompressedSNPBlocks
    M <- bgen_header$M
    N <- bgen_header$N

    out <- load_bgen_sample_identifier_block(to.read, SampleIdentifiers, N)     
    seek(to.read, where =NA)    
    
    ##
    for(i_m in 1:M) {
        out <- load_variant_identifying_data_for_one_snp(to.read, Layout)
        variant_info <- out$variant_info
        num_K_alleles <- out$num_K_alleles
    }

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
    rawToChar(magic) ## yay
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
    }
    L_SI <- 8 + 2 * N + sum(L_si) ## length of sample identifier block
    return(
        list(
            sample_names = sample_names,
            L_SI = L_SI
        )
    )
}

load_variant_identifying_data_for_one_snp <- function(to.read, Layout) {
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
    num_K_alleles <- readBin(to.read, size = 2, "integer", n = 1, endian = "little")
    alleles <- array(NA, num_K_alleles)
    for(i_allele in 1:num_K_alleles) {
        L_ai <- readBin(to.read, size = 4, "integer", n = 1, endian = "little")
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
    return(
        list(
            variant_info = variant_info,
            num_K_alleles = num_K_alleles
        )
    )
}


load_genotypes_for_one_snp <- function(to.read, num_K_alleles, N, CompressedSNPBlocks) {
    ## Genotype data block
    C <- readBin(to.read, size = 4, "integer", n = 1, endian = "little") ## length of rest of the data for this variant
    if (CompressedSNPBlocks > 0) {
        D <- readBin(to.read, size = 4, "integer", n = 1, endian = "little") ## length after de-compression
    } else {
        D <- C
    }
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
    bit_prob <- as.integer(readBin(dataC, size = 1, "raw", n = 1, endian = "little"))
    if ((bit_prob < 1) | (32 < bit_prob)) {
        stop("bit prob is outside range")
    }
    if (phased_flag == 1) {
        ## haplotype probabilities
        stop("haplotype parsing not written")
    } else if (phased_flag == 0) {
        ##
        remaining_bytes <- D - 4 - 2 - 1 - 1 - N - 1 - 1
        dosage <- array(NA, N)
        gen_probs <- array(NA, c(3, N))
        if ((bit_prob %in% c(8, 16, 24, 32)) == FALSE) {
        stop("non multiple of 8 bit_prob not supported")
        }
        for(iSample in 1:N) {
            if (ploidy[iSample] != 2) {
                stop("this code does not support non-2 ploidy")
            }
            b_hom_ref <- readBin(dataC, size = 1, "raw", n = bit_prob / 8, endian = "little")
            b_het <- readBin(dataC, size = 1, "raw", n = bit_prob / 8, endian = "little")        
            if (is_missing[iSample]) {
                dosage[iSample] <- NA
                gen_probs[1:3, iSample] <- NA
            } else {
                x <- sum(as.integer(rawToBits(b_hom_ref)) * (2 ** (0:(bit_prob - 1)))    )
                p_hom_ref <- x / (2 ** bit_prob)
                x <- sum(as.integer(rawToBits(b_het)) * (2 ** (0:(bit_prob - 1)))    )
                p_het <- x / (2 ** bit_prob)
                p_hom_alt <- 1 - p_hom_ref - p_het
                dosage[iSample] <- p_het * 2 * p_hom_alt
                gen_probs[1:3, iSample] <- c(p_hom_ref, p_het, p_hom_alt)
            }
        }
    }
    return(
        list(
            dosage = dosage,
            gen_probs = gen_probs
        )
    )
}
