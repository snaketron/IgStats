source("R/Posterior.R")
source("R/Stats.R")
source("R/Util.R")


require(rstan)
rstan_options(auto_write = TRUE)


d <- read.csv(file = "cells.data.tsv", header = T,
              as.is = T, sep = "\t")

d <- d[, c("sample", "condition", "cdr3.aa")]
d$cdr3.sequence <- d$cdr3.aa
d$cdr3.aa <- NULL
cdr3.data <- d


stats.length <- computeCdr3LengthStats(cdr3.data = cdr3.data,
                                       hdi.level = 0.95,
                                       cores = 4)
