################################################################################
# Replication file for Cattaneo, Ma, Masatlioglu and Suleymanov (2020)
# 20-OCT-2021
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

# The whole exercise takes 10-20 minutes

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
    temp <- logitSimu(n, 5, j, varsig)
    menu <- rbind(menu, temp$menu); choice <- rbind(choice, temp$choice)
  }

  # Inference result
  temp_ram <- revealPref(menu, choice, pref_list = pref_list, method = "GMS",
                         nCritSimu = 2000,
                         BARatio2MS = 0.1, BARatio2UB = 0.1, MNRatioGMS = NULL,
                         RAM = TRUE, AOM = FALSE,
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
