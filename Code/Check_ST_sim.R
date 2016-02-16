# clear environment
rm(list = ls())
gc()

if( Sys.getenv("USERNAME") %in% c("James.Thorson","xJames.Thorson") ) setwd( "C:/Users/James.Thorson/Desktop/Project_git/Trout_GRF/" )

#######################
# Load libraries
#######################
library(TMB)
library(minqa)
library(dplyr)
library(lubridate)
library(tidyr)
library(ggplot2)
source("Functions/Input_Functions.R")
source("Functions/simOUGMRF.R")
source("Functions/simST.R")
source("Functions/runOUGMRF.R")
source("Functions/summary_functions.R")

dir_out <- "Output"
if(!file.exists("Output/Power_Sim/Data/")) {
  dir.create("Output/Power_Sim/Data/", recursive = TRUE)
}

#######################
# Load data
#######################
load( "Data/White_River_Network.RData")
colnames(family)[1] = "child_name"
family = cbind( family, "child_b"=1:nrow(family) )

# Set Conditions
mean_N <- 10
n_years_vec <- c(20)
n_years <- max(n_years_vec)
sample_sites_vec <- c(nrow(family))
p <- c(0.75, 0.75, 0.75)
theta <- 1
SD <- 0.3
rhot <- 0.5
SD_t <- 1
theta_st <- c(0.25, 0.5, 0.75, 1)
SD_st <- c(0.05, 0.1, 0.2, 0.3)
rho <- c(0.25, 0.5, 0.75)

n_sim <- 25

# Covariates
# add spatially varying covariates constant in time
gamma_j <- c(0.5) # doesn't work if I use more than 1 coef
X_i <- matrix(rnorm(nrow(family), 0, 1), nrow(family), length(gamma_j))
# replicate by number of years
X_ij <- do.call("rbind", rep(list(X_i), n_years))

# Set TMB code
Version = "OU_GMRF_v1h"

# Sanity checks
if( TRUE ){
  network <- simST(family, theta = theta, SD = SD, rhot = rhot, SD_t = 0, theta_st = theta_st[1], SD_st = SD_st[1], mean_N = mean_N, n_years = n_years, rho = rho[1], gamma_j = gamma_j, X_ij=X_ij, p = p, spatial = FALSE, temporal = TRUE, spatiotemporal = TRUE)
  # Check average sample AR of log-density at each site
  # This is only equal to rhot, or rho, when the other process has variance (i.e., SD_t or SD_st) fixed at zero, and is otherwise some kind of weighted average of the two (I think)
  ar1 = function(vec) var(vec[-length(vec)],vec[-1]) / var(vec)
  mean(apply(network$log_Npred_bt, MARGIN=1, FUN=ar1))
}
