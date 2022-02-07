#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

source("benchmark_analysis.R")
batch <- as.numeric(args[1])
batch.size <- 1
run(batch, batch.size)
