

# Description:
computeCdr3Stats <- function(cdr3.data, hdi.level = 0.95, cores) {


  getLengthStats <- function(cdr3.data, stan.data,
                             hdi.level, cores) {

    # model
    message("Compiling model ... \n")
    model.file <- "inst/extdata/cdr3_length.stan"
    # model.file <- system.file("extdata",
    #                           "cdr3_length.stan",
    #                           package = "IgStats")
    model <- rstan::stan_model(file = model.file,
                               auto_write = TRUE)


    # stan sampling
    pars.relevant <- c("mu_sample", "sigma_sample",
                       "mu_condition", "sigma_condition",
                       "Yhat_sample")
    glm <- rstan::sampling(object = model,
                           data = stan.data,
                           chains = 4,
                           cores = cores,
                           iter = 5000,
                           warmup = 2500,
                           control = list(adapt_delta = 0.99,
                                          max_treedepth = 12),
                           pars = pars.relevant,
                           algorithm = "NUTS")

    # get summary
    message("Computing summaries ... \n")
    stats <- getCdr3LengthChargePosterior(glm.cdr3.length = glm,
                                          hdi.level = hdi.level,
                                          cdr3.data = cdr3.data,
                                          stan.data = stan.data)

    return (list(stats = stats, glm = glm, stan.data = stan.data))
  }


  getChargeStats <- function(cdr3.data, stan.data,
                             hdi.level, cores) {

    # model
    message("Compiling model ... \n")
    model.file <- "inst/extdata/cdr3_charge.stan"
    # model.file <- system.file("extdata",
    #                           "cdr3_length.stan",
    #                           package = "IgStats")
    model <- rstan::stan_model(file = model.file,
                               auto_write = TRUE)


    # stan sampling
    pars.relevant <- c("mu_sample", "sigma_sample",
                       "mu_condition", "sigma_condition",
                       "Yhat_sample")
    glm <- rstan::sampling(object = model,
                           data = stan.data,
                           chains = 4,
                           cores = cores,
                           iter = 5000,
                           warmup = 2500,
                           control = list(adapt_delta = 0.99,
                                          max_treedepth = 12),
                           pars = pars.relevant,
                           algorithm = "NUTS")

    # get summary
    message("Computing summaries ... \n")
    stats <- getCdr3LengthChargePosterior(glm.cdr3.length = glm,
                                          hdi.level = hdi.level,
                                          cdr3.data = cdr3.data,
                                          stan.data = stan.data)


    return (list(stats = stats, glm = glm, stan.data = stan.data))
  }



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


  # infer DNA
  cdr3.data <- parseCdr3Data(cdr3.data = cdr3.data)

  stan.data <- getStanFormattedCdr3Data(parsed.cdr3.data = cdr3.data)

  out.length <- getLengthStats(cdr3.data = cdr3.data,
                               stan.data = stan.data$length,
                               hdi.level = hdi.level,
                               cores = cores)

  out.charge.ph6 <- getChargeStats(cdr3.data = cdr3.data,
                                   stan.data = stan.data$charge.6,
                                   hdi.level = hdi.level,
                                   cores = cores)

  out.charge.ph7 <- getChargeStats(cdr3.data = cdr3.data,
                                   stan.data = stan.data$charge.7,
                                   hdi.level = hdi.level,
                                   cores = cores)

  return (list(length = out.length,
               charge.ph6 = out.charge.ph6,
               charge.ph7 = out.charge.ph7))
}

