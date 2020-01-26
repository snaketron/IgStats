

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
  cdr3.data$sample_id <- paste(cdr3.data$sample, cdr3.data$condition, sep = '_')
  cdr3.data$sample_id <- as.numeric(as.factor(cdr3.data$sample_id))
  cdr3.data <- cdr3.data[order(cdr3.data$sample_id, decreasing = FALSE), ]
  cdr3.data$sample_id <- cdr3.data$sample_id


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
getStanFormattedCdr3Data <- function(parsed.cdr3.data) {
  parsed.cdr3.data$C <- as.numeric(as.factor(parsed.cdr3.data$condition))
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
    prior_mean_condition = 15,
    prior_sigma_condition = 10)

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
    prior_mean_condition = 0,
    prior_sigma_condition = 20)

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
    prior_mean_condition = 0,
    prior_sigma_condition = 20)

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

