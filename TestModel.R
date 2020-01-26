source("R/Posterior.R")
source("R/Stats.R")
source("R/Util.R")

require(Biostrings)
require(Peptides)
require(rstan)
rstan_options(auto_write = TRUE)


d <- read.csv(file = "cells.data.tsv", header = T,
              as.is = T, sep = "\t")
d <- d[, c("sample", "condition", "cdr3.aa")]
d$cdr3.sequence <- d$cdr3.aa
d$cdr3.aa <- NULL
cdr3.data <- d
table(cdr3.data$sample)

# cdr3.data <- cdr3.data[cdr3.data$condition == "E", ]
table(cdr3.data$sample)

x <- parseCdr3Data(cdr3.data = cdr3.data)
x <- getStanFormattedCdr3Data(parsed.cdr3.data = x)$length

m_nb <- rstan::stan_model(file = "inst/extdata/cdr3_length_nb.stan", auto_write = TRUE)
m_norm <- rstan::stan_model(file = "inst/extdata/cdr3_length.stan", auto_write = TRUE)


g_nb <- rstan::sampling(object = m_nb,
                        data = x,
                        chains = 2,
                        cores = 2,
                        iter = 2000,
                        warmup = 1000,
                        control = list(adapt_delta = 0.95,
                                       max_treedepth = 10),
                        algorithm = "NUTS")


g_norm <- rstan::sampling(object = m_norm,
                          data = x,
                          chains = 2,
                          cores = 2,
                          iter = 2000,
                          warmup = 1000,
                          control = list(adapt_delta = 0.99,
                                         max_treedepth = 12),
                          algorithm = "NUTS")


require(loo)
loo::loo(g_nb)
loo::loo(g_norm)


d_nb <- data.frame(summary(g_nb, par = "Yhat_sample")$summary)
d_normal <- data.frame(summary(g_nb, par = "Yhat_sample")$summary)
d_nb$obs <- aggregate(x$Y~x$G, FUN = mean)[, 2]
d_normal$obs <- aggregate(x$Y~x$G, FUN = mean)[, 2]

ggplot(data = d_nb)+
  geom_point(aes(y = mean, x = obs))+
  geom_errorbar(aes(y = mean, x = obs, ymin = X2.5., ymax = X97.5.))


ggplot(data = d_normal)+
  geom_point(aes(y = mean, x = obs))+
  geom_errorbar(aes(y = mean, x = obs, ymin = X2.5., ymax = X97.5.))

d <- data.frame(mean.nb = d_nb$mean,
                L.nb = d_nb$X2.5.,
                H.nb = d_nb$X97.5.,
                mean.normal = d_normal$mean,
                L.normal = d_normal$X2.5.,
                H.normal = d_normal$X97.5.,
                obs = aggregate(x$Y~x$G, FUN = mean)[, 2],
                G = aggregate(x$Y~x$G, FUN = mean)[, 1])

ggplot(data = d)+
  geom_point(aes(y = mean.nb, x = mean.normal))+
  geom_errorbar(aes(y = mean.nb, x = mean.normal, ymin = L.nb, ymax = H.nb))+
  geom_errorbarh(aes(y = mean.nb, x = mean.normal, xmin = L.normal, xmax = H.normal))


ggplot(data = d)+
  geom_point(aes(x = G, y = obs), col = "green")+
  geom_point(aes(x = G-0.2, y = mean.nb), col = "red")+
  geom_point(aes(x = G+0.2, y = mean.normal), col = "blue")+
  geom_errorbar(aes(x = G-0.2, y = mean.nb, ymin = L.nb, ymax = H.nb),
                col = "red", width = 0.1)+
  geom_errorbar(aes(x = G+0.2, y = mean.normal, ymin = L.normal, ymax = H.normal),
                col = "blue", width = 0.1)+
  scale_y_continuous(limits = c(12, 17))
