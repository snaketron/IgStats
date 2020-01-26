

# Description:
# * cdr3.data:
#    -cdr3.sequence
#    -sample
#    -condition
# * Computes:
#    -length
#    -net-charge at pH = 6.0
#    -net-charge at pH = 7.0
#    -AA composition
parseCdr3Data <- function(cdr3.data) {

  # Group
  cdr3.data$sample_id <- paste(cdr3.data$replicate, cdr3.data$sample,
                               cdr3.data$condition, sep = '_')
  cdr3.data$sample_id <- as.numeric(as.factor(cdr3.data$sample_id))
  cdr3.data <- cdr3.data[order(cdr3.data$sample_id, decreasing = FALSE), ]

  # S
  cdr3.data$S <- as.numeric(as.factor(cdr3.data$sample))

  # C
  cdr3.data$C <- as.numeric(as.factor(cdr3.data$condition))



  # Length:
  cdr3.data$length <- nchar(x = cdr3.data$cdr3.sequence)

  # Net charge:
  cdr3.data$charge6 <- Peptides::charge(seq = cdr3.data$cdr3.sequence,
                                        pH = 6, pKscale = "Lehninger")
  cdr3.data$charge7 <- Peptides::charge(seq = cdr3.data$cdr3.sequence,
                                        pH = 7, pKscale = "Lehninger")

  # Counts AA:
  aa.counts <- data.frame(letterFrequency(
    x = Biostrings::AAStringSet(cdr3.data$cdr3.sequence),
    letters = Biostrings::AA_STANDARD, as.prob = FALSE))
  colnames(aa.counts) <- paste0("letter_", colnames(aa.counts))

  cdr3.data <- cbind(cdr3.data, aa.counts)

  return (cdr3.data)
}



# Description:
# * shm.data:
#    -Y: mutation count
#    -N: region
#    -region
#    -sample
#    -condition
#    -cell.type
#    -replicate
# * Computes:
#    -get G
parseShmData <- function(shm.data) {

  # Group
  shm.data$sample_id <- paste(shm.data$region, shm.data$sample,
                              shm.data$condition, sep = '_')
  shm.data$sample_id <- as.numeric(as.factor(shm.data$sample_id))
  shm.data <- shm.data[order(shm.data$sample_id, decreasing = FALSE), ]

  return (shm.data)
}



# Description:
# Compute pmax
getPmax <- function(glm.ext, par) {

  getPmaxGene <- function(x, par.data) {
    p <- par.data[,x]
    l <- length(p)
    o <- max(sum(p < 0)/l, sum(p>0)/l)
    return(o)
  }

  par.data <- glm.ext[[par]]
  pmax <- vapply(X = seq_len(length.out = ncol(par.data)),
                 FUN = getPmaxGene,
                 par.data = par.data,
                 FUN.VALUE = numeric(length = 1))
  pmax <- 2 * pmax - 1
  return(pmax)
}




# Description:
# Get stan
getStanCdr3Data_rsc <- function(parsed.cdr3.data) {
  # file: cdr3_length_scr.stan
  # int<lower=0> N;                       // number of data points
  # int<lower=0> Ng;                      // number of samples
  # int<lower=0> Ns;                      // number of subjects
  # int<lower=0> Nc;                      // number of conditions
  # int<lower=0> Gcount [Ng];             // BCRs per sample
  # real <lower=0> Y [N];                 // CDR3 lengths
  # int <lower=0> G [N];                  // sample IDs
  # int <lower=0> C [Ns];                 // condition IDs
  # int <lower=0> S [Ng];                 // subject IDs
  # real <lower=0> prior_mean;            // prior mean in condition (of normal pdf)
  # real <lower=0> prior_sigma;           // prior sigma in condition (of normal pdf)

  unique.cdr3.data <- parsed.cdr3.data[duplicated(parsed.cdr3.data$sample_id) == F, ]
  unique.cdr3.data <- unique.cdr3.data[order(unique.cdr3.data$sample_id, decreasing = F), ]

  unique.sample.data <- parsed.cdr3.data[duplicated(unique.sample.data$S) == F, ]
  unique.sample.data <- unique.sample.data[order(unique.sample.data$S, decreasing = F), ]

  length <- list(
    N = nrow(parsed.cdr3.data),
    Ng = max(parsed.cdr3.data$sample_id),
    Nc = max(unique.cdr3.data$C),
    Ns = max(unique.cdr3.data$S),
    Gcount = as.numeric(table(parsed.cdr3.data$sample_id)),
    Y = parsed.cdr3.data$length,
    C = unique.sample.data$C,
    Corg = unique.sample.data$condition,
    S = unique.sample.data$S,
    Sorg = unique.sample.data$sample,
    G = parsed.cdr3.data$sample_id,
    Gorg = parsed.cdr3.data$sample,
    prior_mean = 15,
    prior_sigma = 10)

  charge.6 <- list(
    N = nrow(parsed.cdr3.data),
    Ng = max(parsed.cdr3.data$sample_id),
    Nc = max(unique.cdr3.data$C),
    Ns = max(unique.cdr3.data$S),
    Gcount = as.numeric(table(parsed.cdr3.data$sample_id)),
    Y = parsed.cdr3.data$charge6,
    C = unique.sample.data$C,
    Corg = unique.sample.data$condition,
    S = unique.sample.data$S,
    Sorg = unique.sample.data$sample,
    G = parsed.cdr3.data$sample_id,
    Gorg = parsed.cdr3.data$sample,
    prior_mean = 0,
    prior_sigma = 20)

  charge.7 <- list(
    N = nrow(parsed.cdr3.data),
    Ng = max(parsed.cdr3.data$sample_id),
    Nc = max(unique.cdr3.data$C),
    Ns = max(unique.cdr3.data$S),
    Gcount = as.numeric(table(parsed.cdr3.data$sample_id)),
    Y = parsed.cdr3.data$charge7,
    C = unique.sample.data$C,
    Corg = unique.sample.data$condition,
    S = unique.sample.data$S,
    Sorg = unique.sample.data$sample,
    G = parsed.cdr3.data$sample_id,
    Gorg = parsed.cdr3.data$sample,
    prior_mean = 0,
    prior_sigma = 20)

  j <- which(regexpr(pattern = "letter", text = colnames(parsed.cdr3.data)) != -1)
  AA <- matrix(data = 0, nrow = max(parsed.cdr3.data$sample_id), ncol = 20)
  for(g in 1:max(parsed.cdr3.data$sample_id)) {
    AA[g, ] <- as.numeric(colSums(x = parsed.cdr3.data[parsed.cdr3.data$sample_id == g, j]))
  }
  colnames(AA) <- gsub(pattern = "letter\\_", replacement = '',
                       colnames(parsed.cdr3.data)[j])

  aa <- list(
    Na = 20,
    Ng = max(parsed.cdr3.data$sample_id),
    Nc = max(unique.cdr3.data$C),
    Y = AA,
    Yorg = AA,
    C = unique.sample.data$C,
    Corg = unique.sample.data$condition,
    S = unique.sample.data$S,
    Sorg = unique.sample.data$sample,
    G = parsed.cdr3.data$sample_id,
    Gorg = parsed.cdr3.data$sample)


  stan.data <- list(length = length,
                    charge.6 = charge.6,
                    charge.7 = charge.7,
                    aa = aa)
  return (stan.data)
}



# Description:
# Get stan
getStanCdr3Data_sc <- function(parsed.cdr3.data) {
  # file: cdr3_length_sc.stan
  # int<lower=0> N;                       // number of data points
  # int<lower=0> Ng;                      // number of all samples
  # int<lower=0> Nc;                      // number of conditions
  # int<lower=0> Gcount [Ng];             // number of groups
  # real <lower=0> Y [N];                 // CDR3 length
  # int <lower=0> C [Ng];                 // condition ID
  # int <lower=0> G [N];                  // sample IDs
  # real <lower=0> prior_mean;            // prior mean in condition (of normal pdf)
  # real <lower=0> prior_sigma;           // prior sigma in condition (of normal pdf)

  unique.cdr3.data <- parsed.cdr3.data[duplicated(parsed.cdr3.data$sample_id) == F, ]
  unique.cdr3.data <- unique.cdr3.data[order(unique.cdr3.data$sample_id, decreasing = F), ]

  length <- list(
    N = nrow(parsed.cdr3.data),
    Ng = max(parsed.cdr3.data$sample_id),
    Nc = max(unique.cdr3.data$C),
    Gcount = as.numeric(table(parsed.cdr3.data$sample_id)),
    Y = parsed.cdr3.data$length,
    C = unique.cdr3.data$C,
    Corg = unique.cdr3.data$condition,
    G = parsed.cdr3.data$sample_id,
    Gorg = parsed.cdr3.data$sample,
    prior_mean = 15,
    prior_sigma = 10)

  charge.6 <- list(
    N = nrow(parsed.cdr3.data),
    Ng = max(parsed.cdr3.data$sample_id),
    Nc = max(unique.cdr3.data$C),
    Gcount = as.numeric(table(parsed.cdr3.data$sample_id)),
    Y = parsed.cdr3.data$charge6,
    C = unique.cdr3.data$C,
    Corg = unique.cdr3.data$condition,
    G = parsed.cdr3.data$sample_id,
    Gorg = parsed.cdr3.data$sample,
    prior_mean = 0,
    prior_sigma = 20)

  charge.7 <- list(
    N = nrow(parsed.cdr3.data),
    Ng = max(parsed.cdr3.data$sample_id),
    Nc = max(unique.cdr3.data$C),
    Gcount = as.numeric(table(parsed.cdr3.data$sample_id)),
    Y = parsed.cdr3.data$charge7,
    C = unique.cdr3.data$C,
    Corg = unique.cdr3.data$condition,
    G = parsed.cdr3.data$sample_id,
    Gorg = parsed.cdr3.data$sample,
    prior_mean = 0,
    prior_sigma = 20)

  j <- which(regexpr(pattern = "letter", text = colnames(parsed.cdr3.data)) != -1)
  AA <- matrix(data = 0, nrow = max(parsed.cdr3.data$sample_id), ncol = 20)
  for(g in 1:max(parsed.cdr3.data$sample_id)) {
    AA[g, ] <- as.numeric(colSums(x = parsed.cdr3.data[parsed.cdr3.data$sample_id == g, j]))
  }
  colnames(AA) <- gsub(pattern = "letter\\_", replacement = '',
                       colnames(parsed.cdr3.data)[j])

  aa <- list(
    Na = 20,
    Ng = max(parsed.cdr3.data$sample_id),
    Nc = max(unique.cdr3.data$C),
    Y = AA,
    Yorg = AA,
    C = unique.cdr3.data$C,
    Corg = unique.cdr3.data$condition,
    G = parsed.cdr3.data$sample_id,
    Gorg = parsed.cdr3.data$sample)


  stan.data <- list(length = length,
                    charge.6 = charge.6,
                    charge.7 = charge.7,
                    aa = aa)
  return (stan.data)
}



# Description:
# Get stan
getStanCdr3Data_s <- function(parsed.cdr3.data) {
  # file: cdr3_length_s.stan
  # int<lower=0> N;                       // number of data points
  # int<lower=0> Ng;                      // number of all samples
  # int<lower=0> Gcount [Ng];             // number of groups
  # real <lower=0> Y [N];                 // CDR3 length
  # int <lower=0> G [N];                  // sample IDs
  # real <lower=0> prior_mean;            // prior mean in samples (of normal pdf)
  # real <lower=0> prior_sigma;           // prior sigma in samples (of normal pdf)

  unique.cdr3.data <- parsed.cdr3.data[duplicated(parsed.cdr3.data$sample_id) == F, ]
  unique.cdr3.data <- unique.cdr3.data[order(unique.cdr3.data$sample_id, decreasing = F), ]

  length <- list(
    N = nrow(parsed.cdr3.data),
    Ng = max(parsed.cdr3.data$sample_id),
    Gcount = as.numeric(table(parsed.cdr3.data$sample_id)),
    Y = parsed.cdr3.data$length,
    Corg = unique.cdr3.data$condition,
    G = parsed.cdr3.data$sample_id,
    Gorg = parsed.cdr3.data$sample,
    prior_mean = 15,
    prior_sigma = 10)

  charge.6 <- list(
    N = nrow(parsed.cdr3.data),
    Ng = max(parsed.cdr3.data$sample_id),
    Gcount = as.numeric(table(parsed.cdr3.data$sample_id)),
    Y = parsed.cdr3.data$charge6,
    Corg = unique.cdr3.data$condition,
    G = parsed.cdr3.data$sample_id,
    Gorg = parsed.cdr3.data$sample,
    prior_mean = 0,
    prior_sigma = 20)

  charge.7 <- list(
    N = nrow(parsed.cdr3.data),
    Ng = max(parsed.cdr3.data$sample_id),
    Gcount = as.numeric(table(parsed.cdr3.data$sample_id)),
    Y = parsed.cdr3.data$charge7,
    G = parsed.cdr3.data$sample_id,
    Gorg = parsed.cdr3.data$sample,
    prior_mean = 0,
    prior_sigma = 20)

  j <- which(regexpr(pattern = "letter", text = colnames(parsed.cdr3.data)) != -1)
  AA <- matrix(data = 0, nrow = max(parsed.cdr3.data$sample_id), ncol = 20)
  for(g in 1:max(parsed.cdr3.data$sample_id)) {
    AA[g, ] <- as.numeric(colSums(x = parsed.cdr3.data[parsed.cdr3.data$sample_id == g, j]))
  }
  colnames(AA) <- gsub(pattern = "letter\\_", replacement = '',
                       colnames(parsed.cdr3.data)[j])

  aa <- list(
    Na = 20,
    Ng = max(parsed.cdr3.data$sample_id),
    Y = AA,
    Yorg = AA,
    G = parsed.cdr3.data$sample_id,
    Gorg = parsed.cdr3.data$sample)


  stan.data <- list(length = length,
                    charge.6 = charge.6,
                    charge.7 = charge.7,
                    aa = aa)
  return (stan.data)
}



