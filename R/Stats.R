

# Description:
computeCdr3LengthStats <- function(cdr3.data,
                                   hdi.level = 0.95,
                                   cores) {

  # check inputs
  # checkInput(usage.data = usage.data,
  #            mcmc.chains = as.integer(x = mcmc.chains),
  #            mcmc.cores = as.integer(x = mcmc.cores),
  #            mcmc.steps = as.integer(x = mcmc.steps),
  #            mcmc.warmup = as.integer(x = mcmc.warmup),
  #            hdi.level = hdi.level)


  # format input
  # usage.data.raw <- usage.data
  # usage.data <- getUsageData(usage = usage.data.raw)


  # infer DNA or
  cdr3.data <- parseCdr3Data(cdr3.data = cdr3.data)
  stan.dl <- getStanFormattedCdr3Data(parsed.cdr3.data = cdr3.data)
  browser()


  # model
  message("Compiling model ... \n")
  model.file <- system.file("extdata",
                            "cdr3_length.stan",
                            package = "IgStats")
  model <- rstan::stan_model(file = model.file,
                             auto_write = TRUE)


  # stan sampling
  pars.relevant <- c("mu_sample",
                     "sigma_sample",
                     "mu_condition",
                     "sigma_condition",
                     "Yhat_sample")
  glm <- rstan::sampling(object =  ,
                         data = stan.dl$length,
                         chains = 4,
                         cores = cores,
                         iter = 5000,
                         warmup = 1500,
                         refresh = 250,
                         control = list(adapt_delta = 0.99,
                                        max_treedepth = 12),
                         pars = pars.relevant,
                         algorithm = "NUTS")


  # get summary
  message("Computing summaries ... \n")
  stats <- getCdr3LengthPosterior(glm.cdr3.length = glm,
                                  hdi.level = hdi.level,
                                  cdr3.data = cdr3.data)

  return (stats)
}



# Description:
computeCdr3ChargeStats <- function(cdr3.data,
                                   hdi.level = 0.95,
                                   cores) {

  # check inputs
  # checkInput(usage.data = usage.data,
  #            mcmc.chains = as.integer(x = mcmc.chains),
  #            mcmc.cores = as.integer(x = mcmc.cores),
  #            mcmc.steps = as.integer(x = mcmc.steps),
  #            mcmc.warmup = as.integer(x = mcmc.warmup),
  #            hdi.level = hdi.level)


  # format input
  # usage.data.raw <- usage.data
  # usage.data <- getUsageData(usage = usage.data.raw)


  # infer DNA or
  cdr3.data <- parseCdr3Data(cdr3.data = cdr3.data)


  # model
  message("Compiling model ... \n")
  model.file <- system.file("extdata",
                            "cdr3_length.stan",
                            package = "IgStats")
  model <- rstan::stan_model(file = model.file,
                             auto_write = TRUE)


  # stan sampling
  pars.relevant <- c("mu_sample",
                     "sigma_sample",
                     "mu_condition",
                     "sigma_condition",
                     "Yhat_sample")
  glm <- rstan::sampling(object = model,
                         data = cdr3.data,
                         chains = 4,
                         cores = cores,
                         iter = 5000,
                         warmup = 1500,
                         refresh = 250,
                         control = list(adapt_delta = 0.99,
                                        max_treedepth = 12),
                         pars = pars.relevant,
                         algorithm = "NUTS")


  # get summary
  message("Computing summaries ... \n")
  stats <- getCdr3ChargePosterior(glm.cdr3.length = glm,
                                  hdi.level = hdi.level,
                                  cdr3.data = cdr3.data)

  return (stats)
}



# Description:
computeCdr3AAStats <- function(cdr3.data,
                               hdi.level = 0.95,
                               cores) {

  # check inputs
  # checkInput(usage.data = usage.data,
  #            mcmc.chains = as.integer(x = mcmc.chains),
  #            mcmc.cores = as.integer(x = mcmc.cores),
  #            mcmc.steps = as.integer(x = mcmc.steps),
  #            mcmc.warmup = as.integer(x = mcmc.warmup),
  #            hdi.level = hdi.level)


  # format input
  # usage.data.raw <- usage.data
  # usage.data <- getUsageData(usage = usage.data.raw)


  # infer DNA or
  cdr3.data <- parseCdr3Data(cdr3.data = cdr3.data)


  # model
  message("Compiling model ... \n")
  model.file <- system.file("extdata",
                            "cdr3_length.stan",
                            package = "IgStats")
  model <- rstan::stan_model(file = model.file,
                             auto_write = TRUE)


  # stan sampling
  pars.relevant <- c("mu_sample",
                     "sigma_sample",
                     "mu_condition",
                     "sigma_condition",
                     "Yhat_sample")
  glm <- rstan::sampling(object = model,
                         data = cdr3.data,
                         chains = 4,
                         cores = cores,
                         iter = 5000,
                         warmup = 1500,
                         refresh = 250,
                         control = list(adapt_delta = 0.99,
                                        max_treedepth = 12),
                         pars = pars.relevant,
                         algorithm = "NUTS")


  # get summary
  message("Computing summaries ... \n")
  stats <- getCdr3AAPosterior(glm.cdr3.length = glm,
                              hdi.level = hdi.level,
                              cdr3.data = cdr3.data)

  return (stats)
}

