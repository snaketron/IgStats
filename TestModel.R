source("R/Posterior.R")
source("R/Stats.R")
source("R/Util.R")

require(Biostrings)
require(Peptides)
require(rstan)
rstan_options(auto_write = TRUE)
require(ggplot2)

d <- read.csv(file = "cells.data.tsv", header = T,
              as.is = T, sep = "\t")
d$replicate <- NA
d <- d[, c("replicate" ,"sample", "condition", "cdr3.aa")]
d$cdr3.sequence <- d$cdr3.aa
d$cdr3.aa <- NULL
cdr3.data <- d



# table(cdr3.data$sample, cdr3.data$condition)
# cdr3.data <- cdr3.data[cdr3.data$sample %in% c("11B_S3", "12B_S7",
#                                                "11E_S6", "12E_S8",
#                                                "1N_S15", "3N_S17"), ]


x <- parseCdr3Data(cdr3.data = cdr3.data)
x <- getStanFormattedCdr3Data(parsed.cdr3.data = x)$length

m <- rstan::stan_model(file = "inst/extdata/cdr3_length_s.stan",
                       auto_write = TRUE)

g <- rstan::sampling(object = m,
                     data = x,
                     chains = 2,
                     cores = 2,
                     iter = 2000,
                     warmup = 1000,
                     control = list(adapt_delta = 0.99,
                                    max_treedepth = 12),
                     algorithm = "NUTS")

