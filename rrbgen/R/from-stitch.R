## does not need to be the same as STITCH, but convenient
## is a way to specify how to spread N samples of nCores coresn
getSampleRange <- function(
  N,
  nCores
) {
    ## upon closer inspection, this might be slow? oh well
    if (nCores == 1) {
        x <- rep(1, N)
    } else {
        x <- as.integer(cut(1:N, nCores))
    }
    w <- which(diff(x) > 0)
    start <- c(1, w + 1)
    end <- c(w, N)
    sampleRange <- lapply(1:nCores, function(i_core) {
        x <- c(start[i_core], end[i_core])
        if (sum(is.na(x)) > 0)
            return(NULL)
        return(x)
    })
    sampleRange <- sampleRange[sapply(sampleRange, is.null) == FALSE]
    return(sampleRange)
}
