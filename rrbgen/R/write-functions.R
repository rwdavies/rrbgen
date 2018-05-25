write_bgen_header <- function(
    to.write,
    offset,
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
    if ((Layout %in% c(0, 1, 2)) == FALSE) {
        stop(paste0("Layout must be either 0, 1 or 2 but you have selected:", Layout))
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
    return(NULL)
}


write_bgen_sample_identifier_block <- function(to.write, binary_start, sample_names, SampleIdentifiers = 1) {
    ##
    seek(to.write, where = binary_start)
    if (SampleIdentifiers == 1) {
        ## first, convert and build rest
        N <- length(sample_names)        
        sample_names_as_raw <- lapply(sample_names, charToRaw)
        L_Si <- lapply(sample_names_as_raw, length) ## per-sample
        L_SI <- as.integer(8 + 2 * N + sum(unlist(L_Si))) ## total
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
