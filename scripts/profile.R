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


setwd("~/Dropbox/rrbgen/rrbgen/R/")
source("read-functions.R")
source("write-functions.R")
source("test-drivers.R")
source("from-stitch.R")
setwd("~/Dropbox/rrbgen/")


nCores <- 1
CompressedSNPBlocks <- 1
B_bit_prob <- 16


## make really big file
sample_data_file <- "~/Downloads/profile_data.RData"
file <- tempfile()
if (file.exists(sample_data_file) == FALSE) {
    print_message("Prepare sample data")
    M <- 10000
    N <- 2000
    sample_names <- paste0("samp", 1:N)
    var_info <- make_fake_var_info(M)
    var_ids <- var_info[, "varid"]
    gp <- make_fake_gp(sample_names, var_ids, random_fraction = 0)
    list_of_gp_raw_t <- convert_gp_to_list_of_raw(
        gp,
        nCores = nCores,
        B_bit_prob = B_bit_prob
    )
    save(list_of_gp_raw_t, var_info, sample_names, file = sample_data_file)
} else {
    load(sample_data_file)
}

print_message("begin profile ")
profout <- tempfile()
Rprof(file = profout, gc.profiling = TRUE, line.profiling = TRUE)

########################### profile this
rrbgen_write(
    bgen_file = tempfile(),
    sample_names = sample_names,
    var_info = var_info,
    list_of_gp_raw_t = list_of_gp_raw_t,
    CompressedSNPBlocks = CompressedSNPBlocks,
    B_bit_prob = B_bit_prob,
    verbose = TRUE,
    nCores = 1
)
########################### end of what to profile
Rprof(NULL)
print_message("end profile")


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
