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
source("Functions/simST.R")
source("Functions/runOUGMRF.R")
source("Functions/summary_functions.R")


#######################
# Load data
#######################
load( "Data/White_River_Network.RData")
colnames(family)[1] = "child_name"
family = cbind( family, "child_b"=1:nrow(family) )

# Set TMB code
Version = "OU_GMRF_v1i"

# Number of simulations
n_sim <- 5

# Set Conditions
mean_N <- 10 
n_years_vec <- 20
n_years <- max(n_years_vec)
sample_sites_vec <- 300
p <- c(0.5, 0.5, 0.5)  # Detection probability for each of three passes
theta <- 0.5 # range for spatial variation
SD <- 0   #  Marginal SD of spatial variation
rhot <- 0.5 # 0.1 # Autocorrelation over time
SD_t <- 0.3 # 0  # Conditional SD for variation over time
theta_st <- 0.5 # Range for spatio-temporal variation
SD_st <- 0 # Marginal SD of spatial component of spatio-temporal variation
rho <- 0    # Correlation among years for spatio-temporal variation
gamma_j <- c(0.2) # doesn't work if I use more than 1 coef
X_i <- matrix(rnorm(nrow(family), 0, 1), nrow(family), length(gamma_j))
X_ij <- do.call("rbind", rep(list(X_i), n_years))

Options_vec = c("SpatialTF"=0, "TemporalTF"=1, "SpatiotemporalTF"=0, "DetectabilityTF"=0, "ObsModel"=1, "OverdispersedTF"=0, "abundTF"=0)

# Compile
if(FALSE) {
  dyn.unload(dynlib(paste0("Code/", Version)))
  file.remove( paste0("Code/", Version,c(".o",".dll")) )
}
compile( paste0("Code/", Version,".cpp") )

#######################
# Simulations
#######################
df_sims <- data.frame(iter=integer(),
                      theta=numeric(),
                      theta_hat=numeric(),
                      rhot=numeric(),
                      rhot_hat=numeric(),
                      sigmat=numeric(),
                      sigmat_hat=numeric(),
                      theta_st=numeric(),
                      theta_st_hat=numeric(),
                      SD = numeric(),
                      SD_hat = numeric(),
                      SD_st = numeric(),
                      SD_st_hat = numeric(),
                      SD_inf = numeric(),
                      SD_inf_hat = numeric(),
                      SD_st_inf = numeric(),
                      SD_st_inf_hat = numeric(),
                      rho_st=numeric(),
                      rho_st_hat=numeric(),
                      gamma_j=numeric(),
                      gamma_j_hat=numeric(),
                      stringsAsFactors=FALSE)

counter <- 0

# start loop
for(i in 1:n_sim) {
  # set conditions
  df_N <- matrix(NA, nrow(family), length(sample_sites_vec))
  df_coef <- list()
  df_sd <- list()
  df_rmse <- NA
  N_hat <- list()
  
  # simulate abundance and counts on network
  network <- simST(family, theta = theta, SD = SD, rhot = rhot, SD_t = SD_t, theta_st = theta_st, SD_st = SD_st, mean_N = mean_N, n_years = n_years, rho = rho, gamma_j = gamma_j, X_ij=X_ij, p = p, spatial = FALSE, temporal = TRUE, spatiotemporal = FALSE)
  
  # thin to sample sites and years
  counter <- counter + 1
  c_ip <- network$c_ip
  n <- length(network[["N_i"]])
  c_ip$t_i <- network$t_i
  c_ip$child_b <- rep(family$child_b, times = n_years)
  
  # Sample locations
  sites <- 1:max(family$child_b)
  remove_sites <- sample(sites, size = max(family$child_b) - sample_sites_vec[1], replace = FALSE) 
  c_ip[which(c_ip$child_b %in% remove_sites), 1:3] <- NA
  
  # Sample Years
  years <- 1:n_years_vec[1]
  c_ip[which(!(c_ip$t_i %in% years)), 1:3] <- NA
  c_ip <- as.matrix(c_ip[ , 1:length(p)])
  
  #--------- Fit Model in TMB ----------
  start <- 1
  end <- 2
  Calc_lambda_ip <- rep(NA, length.out = nrow(network$c_ip))
  Calc_lambda_ip[start:end] <- 1
  Calc_lambda_ip[is.na(Calc_lambda_ip)] <- 0
  
  # need to make df, family, and c_ip agree
  df_counts <- as.data.frame(c_ip)
  df_counts$child_b <- rep(1:nrow(family), times = n_years)
  df <- dplyr::left_join(df_counts, family)
  
  # Make inputs
  Inputs <- makeInput(family = family, df = df, c_ip = c_ip, options = Options_vec, X = X_ij, t_i = network$t_i, version = Version, CalcSD_lambda_ip = Calc_lambda_ip, spatial_equal = TRUE)
  
  # Run
  mod <- runOUGMRF(inputs = Inputs)
  
  # Debugging info -- added by Jim but turned off by default
  if(FALSE){
    # Diagnostics -- estimation model
    unlist(mod$Report[c('log_mean','detectrate','extradetectionSD','rho_st','rhot','SD_inf','SD_st_inf','SDinput','SDinput_st','sigmaIID','sigmat','theta','theta_st')])
    sapply( mod$ParHat[c('log_extradetectrate_i','Epsiloninput_d','Deltainput_t','Nu_dt')], FUN=function(vec){c(mean(vec),sd(vec))} )
    table(names(mod$obj$env$last.par.best))
    mod$Report$rho_t_b
    mod$Report$SDinput_t_b
    mod$ParHat$Nu_dt
    
    # Diagnostics -- simulation model
    sapply( network[c('x_t','x_b','x_bt')], FUN=function(vec){c(mean(vec),sd(vec))} )
    network$x_bt
  }
  
  #----------- summarize -----------
  
  if(class(mod) == "try-error") {
    df_sims[counter, "iter"] <- i
    df_sims[counter, "n_sites"] <- sample_sites_vec[1]
    df_sims[counter, "n_years"] <- n_years_vec[1]
    df_sims[counter, "theta"] <- theta
    df_sims[counter, "theta_hat"] <- NA_real_
    df_sims[counter, "rhot"] <- rhot
    df_sims[counter, "rhot_hat"] <- NA_real_
    df_sims[counter, "sigmat"] <- SD_t
    df_sims[counter, "sigmat_hat"] <- NA_real_
    df_sims[counter, "theta_st"] <- theta_st
    df_sims[counter, "theta_st_hat"] <- NA_real_
    df_sims[counter, "SD"] <- SD
    df_sims[counter, "SD_hat"] <- NA_real_
    df_sims[counter, "SD_st"] <- SD_st
    df_sims[counter, "SD_st_hat"] <- NA_real_
    df_sims[counter, "SD_inf"] <- SD / ((2 * theta) ^ 0.5)
    df_sims[counter, "SD_inf_hat"] <- NA_real_
    df_sims[counter, "SD_st_inf"] <- SD_st / ((2 * theta_st) ^ 0.5)
    df_sims[counter, "SD_st_inf_hat"] <- NA_real_
    df_sims[counter, "rho_st"] <- rho
    df_sims[counter, "rho_st_hat"] <- NA_real_
    df_sims[counter, "gamma_j"] <- gamma_j
    df_sims[counter, "gamma_j_hat"] <- NA_real_
    df_sims[counter, "converge"] <- FALSE
  } else {
    # check convergence
    converge <- FALSE
    try(converge <- mod$opt$convergence == 0)
    try(converge <- ifelse(converge == TRUE, !any(is.na(mod$SD$sd)), converge))
    large_sd <- FALSE
    try(large_sd <- max(mod$SD$sd, na.rm = T) > 100)
    
    if(converge == FALSE | large_sd == TRUE) {
      df_sims[counter, "iter"] <- i
      df_sims[counter, "n_sites"] <- sample_sites_vec[1]
      df_sims[counter, "n_years"] <- n_years_vec[1]
      df_sims[counter, "theta"] <- theta
      df_sims[counter, "theta_hat"] <- NA_real_
      df_sims[counter, "rhot"] <- rhot
      df_sims[counter, "rhot_hat"] <- NA_real_
      df_sims[counter, "sigmat"] <- SD_t
      df_sims[counter, "sigmat_hat"] <- NA_real_
      df_sims[counter, "theta_st"] <- theta_st
      df_sims[counter, "theta_st_hat"] <- NA_real_
      df_sims[counter, "SD"] <- SD
      df_sims[counter, "SD_hat"] <- NA_real_
      df_sims[counter, "SD_st"] <- SD_st
      df_sims[counter, "SD_st_hat"] <- NA_real_
      df_sims[counter, "SD_inf"] <- SD / ((2 * theta) ^ 0.5)
      df_sims[counter, "SD_inf_hat"] <- NA_real_
      df_sims[counter, "SD_st_inf"] <- SD_st / ((2 * theta_st) ^ 0.5)
      df_sims[counter, "SD_st_inf_hat"] <- NA_real_
      df_sims[counter, "rho_st"] <- rho
      df_sims[counter, "rho_st_hat"] <- NA_real_
      df_sims[counter, "gamma_j"] <- gamma_j
      df_sims[counter, "gamma_j_hat"] <- NA_real_
      df_sims[counter, "converge"] <- FALSE
    } else {
      df_sims[counter, "iter"] <- i
      df_sims[counter, "n_sites"] <- sample_sites_vec[1]
      df_sims[counter, "n_years"] <- n_years_vec[1]
      df_sims[counter, "theta"] <- theta
      df_sims[counter, "theta_hat"] <- mod$Report$theta
      df_sims[counter, "rhot"] <- rhot
      df_sims[counter, "rhot_hat"] <- mod$Report$rhot
      df_sims[counter, "sigmat"] <- SD_t
      df_sims[counter, "sigmat_hat"] <- mod$Report$sigmat
      df_sims[counter, "theta_st"] <- theta_st
      df_sims[counter, "theta_st_hat"] <- mod$Report$theta_st
      df_sims[counter, "SD"] <- SD
      df_sims[counter, "SD_hat"] <- mod$Report$SDinput
      df_sims[counter, "SD_st"] <- SD_st
      df_sims[counter, "SD_st_hat"] <- mod$Report$SDinput_st
      df_sims[counter, "SD_inf"] <- SD / ((2 * theta) ^ 0.5)
      df_sims[counter, "SD_inf_hat"] <- mod$Report$SD_inf
      df_sims[counter, "SD_st_inf"] <- SD_st / ((2 * theta_st) ^ 0.5)
      df_sims[counter, "SD_st_inf_hat"] <- mod$Report$SD_st_inf
      df_sims[counter, "rho_st"] <- rho
      df_sims[counter, "rho_st_hat"] <- mod$Report$rho_st
      df_sims[counter, "gamma_j"] <- gamma_j
      df_sims[counter, "gamma_j_hat"] <- mod$Report$gamma_j
      df_sims[counter, "converge"] <- TRUE
    }
  }
} # end sim iter

####################
# Results
####################

# rhot_hat should equal rhot consistently
df_sims %>%
  dplyr::select(iter, rhot, rhot_hat, sigmat, sigmat_hat, converge) %>%
  summary()

