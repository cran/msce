library(testthat)
library(RcppParallel)

cat("Starting tests for 'msce' package ...")

RcppParallel::setThreadOptions(numThreads = 2) # nice to CRAN

test_check("msce")

RcppParallel::setThreadOptions(numThreads = "auto")
