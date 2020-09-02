################################################################################
# Replication file for Cattaneo, Ma, Masatlioglu and Suleymanov (2019)
# 28-AUG-2019
################################################################################

# This replication file provides R code to test against 5 preferences for 1000 
#   simulations. Effective sample size is 200. 

################################################################################
# Load package
################################################################################

# The R code provided relies on the package "ramchoice", which is available on
#   CRAN. Run the following if it is not already installed.
#
#                 install.packages("ramchoice")

rm(list=ls(all=TRUE))
library("ramchoice")

################################################################################
# Generate Data
################################################################################

# This function simulates dataset under the Logit attention rule
#   n    : effective sample size
#   uSize: number of alternatives
#   mSize: size of menu/choice problem
#   par  : parameter

genData <- function(n, uSize, mSize, varsig) {
  # determine population choice rule
  prob_vec <- rep(0, mSize)
  for (i in 1:mSize) { # enumerate over alternatives
    for (j in 1:(mSize-i+1)) { # enumerate over consideration set size
      temp1 <- factorial(mSize-i) / factorial(j-1) / factorial(mSize-i-j+1)
      temp2 <- j^varsig / 
        sum((1:mSize)^varsig * factorial(mSize) / factorial(1:mSize) / factorial((mSize-1):0))
      prob_vec[i] = prob_vec[i] + temp1 * temp2
    }
  }
  # initialize
  allMenus <- t(combn(uSize, mSize))
  menu <- choice <- matrix(0, nrow=n*nrow(allMenus), ncol=uSize) 
  
  for (i in 1:nrow(allMenus)) {
    for (j in 1:n) {
      menu[j+(i-1)*n, allMenus[i, ]] <- 1
      choice[j+(i-1)*n, sort(allMenus[i, ])[rmultinom(1, 1, prob_vec) == 1]] <- 1
    }
  }
  return(list(menu=menu, choice=choice))
}

################################################################################
# Prepare for simulation
################################################################################

# Simulation parameters
n      <- 200    # Effective sample size
repe   <- 1000   # Number of Monte Carlo repetitions
varsig <- 2      # Model parameter

# Following is the list of preferences/null hypotheses. Each row represents one
#   preference ordering.
pref_list <- matrix(c(1, 2, 3, 4, 5,
                      2, 3, 4, 5, 1,
                      3, 4, 5, 2, 1,
                      4, 5, 3, 2, 1,
                      5, 4, 3, 2, 1), ncol=5, byrow=TRUE)

# Following is an empty matrix used to store inference result. Each row corresponds
#   to one Monte Carlo repetition, and each column represents one preference.
result <- matrix(0, ncol=5, nrow=repe)

################################################################################
# Inference
################################################################################

# The whole exercise takes 8 minutes (0.32s for one Monte Carlo repetition). 
# System: MacOS 10.13.1, R 3.4.1 
#   2017 Macbook Pro 13inch, 2.3 GHz Intel Core i5, 16 GB 2133 MHz LPDDR3 

# Final result:
#   --------------------------------
#   Hypothesis      | Empirical size
#   ----------------+---------------
#   Null            |
#     1, 2, 3, 4, 5 | 0.000
#     2, 3, 4, 5, 1 | 0.004
#   ----------------+---------------
#   Alternative     |
#     3, 4, 5, 2, 1 | 0.144
#     4, 5, 3, 2, 1 | 0.275
#     5, 4, 3, 2, 1 | 0.351
#   --------------------------------

set.seed(42) # Seed
ptm <- proc.time() # Track time used

for (i in 1:repe) { # iterate through 2000 simulated datasets

  # Generate data
  menu <- choice <- matrix(0, nrow=0, ncol=5)
  for (j in 5:2) {
    temp <- genData(n, 5, j, varsig)
    menu <- rbind(menu, temp$menu); choice <- rbind(choice, temp$choice)
  }

  # Inference result
  temp_ram <- rAtte(menu, choice, pref_list = pref_list, method = "GMS",
                    nCritSimu = 2000,
                    BARatio2MS = 0.1, BARatio2UB = 0.1, MNRatioGMS = NULL,
                    limDataCorr = FALSE,
                    attBinary = 1)

  # Add inference result to the "result" matrix. We use the proposed critical value,
  #   and nominal level 0.05
  result[i, ] <- 1 * (temp_ram$Tstat > temp_ram$critVal$GMS[, 2])

  # Tracking progress
  if (i %% 50 == 0) cat(paste(i/10, "% completed.\n", sep=""))
}

proc.time() - ptm # Track time used

# Final result
colnames(result) <- c("H1", "H2", "H3", "H4", "H5")
colMeans(result)
