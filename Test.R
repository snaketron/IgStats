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
d <- d[, c("sample", "condition", "cdr3.aa")]
d$cdr3.sequence <- d$cdr3.aa
d$cdr3.aa <- NULL
cdr3.data <- d

table(cdr3.data$sample, cdr3.data$condition)
cdr3.data <- cdr3.data[cdr3.data$sample %in% c("11B_S3", "12B_S7",
                                               "11E_S6", "12E_S8",
                                               "1N_S15", "3N_S17"), ]


stats <- computeCdr3Stats(cdr3.data = cdr3.data,
                          hdi.level = 0.95,
                          cores = 2)

ggplot(data = stats$length$stats$summary.sample)+
  geom_point(aes(x = condition, y = mean, group = sample, col = condition),
             position = position_dodge(width = 0.65))+
  geom_errorbar(aes(x = condition, y = mean, group = sample,
                    col = condition, ymin = mean_L, ymax = mean_H),
             position = position_dodge(width = 0.65), width = 0.3)+
  theme_bw()+
  theme(legend.position = "top")


ggplot(data = stats$length$stats$summary.condition)+
  geom_point(aes(x = condition, y = mean, col = condition))+
  geom_errorbar(aes(x = condition, y = mean, col = condition,
                    ymin = mean_L, ymax = mean_H), width = 0.2)+
  theme_bw()+
  theme(legend.position = "top")
