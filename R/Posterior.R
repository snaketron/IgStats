
getCdr3LengthChargePosterior <- function(glm.cdr3.length,
                                         hdi.level,
                                         cdr3.data,
                                         stan.data) {

  summary.sample <- rstan::summary(object = glm.cdr3.length,
                                   digits = 4, pars = "mu_sample",
                                   prob = c(0.5, (1-hdi.level)/2,
                                            1-(1-hdi.level)/2))
  summary.sample <- summary.sample$summary
  summary.sample <- data.frame(summary.sample)
  colnames(summary.sample) <- c("mean", "mean_se",
                                "mean_sd", "mean_median",
                                "mean_L", "mean_H",
                                "Neff", "Rhat")
  summary.sample$sample_id <- unique(stan.data$G)
  summary.sample$sample <- unique(stan.data$Gorg)
  summary.sample$condition <- stan.data$Corg

  summary.condition <- rstan::summary(object = glm.cdr3.length,
                                      digits = 4, pars = "mu_condition",
                                      prob = c(0.5, (1-hdi.level)/2,
                                               1-(1-hdi.level)/2))
  summary.condition <- summary.condition$summary
  summary.condition <- data.frame(summary.condition)
  colnames(summary.condition) <- c("mean", "mean_se",
                                   "mean_sd", "mean_median",
                                   "mean_L", "mean_H",
                                   "Neff", "Rhat")
  summary.condition$condition <- NA
  for(i in 1:nrow(summary.condition)) {
    summary.condition$condition[i] <- stan.data$Corg[stan.data$C == i][1]
  }

  # return
  return (list(summary.sample = summary.sample,
               summary.condition = summary.condition))
}



getCdr3AAPosterior <- function(glm.cdr3.aa, hdi.level, cdr3.data) {

  cdr3.data.unique <- cdr3.data[duplicated(cdr3.data$sample_id) == FALSE, ]

  summary.sample <- rstan::summary(object = glm.cdr3.aa,
                                   digits = 4, pars = "mu_sample",
                                   prob = c(0.5, (1-hdi.level)/2,
                                            1-(1-hdi.level)/2))
  summary.sample <- summary.sample$summary
  summary.sample <- data.frame(summary.sample)
  colnames(summary.sample) <- c("mean", "mean_se",
                                "mean_sd", "mean_median",
                                "mean_L", "mean_H",
                                "Neff", "Rhat")
  summary.sample$sample_id <- rep(x = cdr3.data.unique$sample_id, each = 20)
  summary.sample$sample <- rep(x = cdr3.data.unique$sample, each = 20)
  summary.sample$condition <- rep(x = cdr3.data.unique$condition, each = 20)
  summary.sample$AA <- rep(x = Biostrings::AA_STANDARD,
                           times = max(cdr3.data.unique$sample_id))


  summary.condition <- rstan::summary(object = glm.cdr3.aa,
                                      digits = 4, pars = "mu_condition",
                                      prob = c(0.5, (1-hdi.level)/2,
                                               1-(1-hdi.level)/2))
  summary.condition <- summary.condition$summary
  summary.condition <- data.frame(summary.condition)
  colnames(summary.condition) <- c("mean", "mean_se",
                                   "mean_sd", "mean_median",
                                   "mean_L", "mean_H",
                                   "Neff", "Rhat")
  cs <- NA
  for(i in 1:nrow(summary.condition)) {
    cs <- c(cs, rep(x = cdr3.data.unique$condition[cdr3.data.unique$C == i][1], each = 20))
  }
  summary.condition$condition <- cs
  summary.condition$AA <- rep(x = Biostrings::AA_STANDARD,
                              times = max(cdr3.data.unique$C))

  # return
  return (list(summary.sample = summary.sample,
               summary.condition = summary.condition))
}
