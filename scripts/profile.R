#!/usr/bin/env Rscript

library("proftools")
library("rrbgen")

## change directory to one up from scripts, no matter how this was called
args <- commandArgs(trailingOnly = FALSE)
for(key in c("--file=", "--f=")) {
    i <- substr(args, 1, nchar(key)) == key
    if (sum(i) == 1) {
        script_dir <- dirname(substr(args[i], nchar(key) + 1, 1000))
        setwd(file.path(script_dir, "../"))        
    }
}


setwd("~/Dropbox/rrbgen/rrbgen/R/"); source("read-functions.R") ;source("write-functions.R"); source("test-drivers.R");setwd("~/Dropbox/rrbgen/")



profout <- tempfile()
Rprof(file = profout, gc.profiling = TRUE, line.profiling = TRUE)

## what to profile
message(date())
external_bgen_dir <- "./external/bgen/"
bit <- 16
external_bgen_file <- file.path(external_bgen_dir, paste0("example.", bit, "bits.bgen"))
message(paste0(date(), ":Start"))
for(i in 1:5) {
    message(paste0(date(), ": Load"))
    out <- rrbgen_load(external_bgen_file)
    message(paste0(date(), ": Write"))
    ## now write as well!
    gp <- out$gp
    sample_names <- out$sample_names
    var_info <- out$var_info
    CompressedSNPBlocks <- 1
    bgen_file <- tempfile()
    B_bit_prob <- 16
    rrbgen_write(
        bgen_file,
        sample_names = sample_names,
        var_info = var_info,
        gp = gp,
        CompressedSNPBlocks = CompressedSNPBlocks,
        B_bit_prob = B_bit_prob
    )
}
message(paste0(date(), ":Done"))
## end of what to profile

Rprof(NULL)
pd <- readProfileData(profout)
title <- Sys.getenv("TITLE")

output_plot <- "profile.pdf"
pdf(output_plot, height = 24, width = 24)
par(mfrow = c(3, 1), oma=c(0, 0, 3, 0))
flameGraph(pd, order = "time", main = "Time")
flameGraph(pd, order = "alpha", main = "Alphabetically")
flameGraph(pd, order = "hot", main = "Hot path")
title(title, outer=TRUE)
dev.off()