# 1. Proof it works
# Spatial model varying theta and sigma. I think it makes sense to compare these to a non-spatial model (sorry a sort of 3rd axis). Do you think we need to vary sigma since it affects the spatial correlation but as a constant and not in relation to distance? I don't think we need to test temporal or spatiotemporal components here.
# 
# 2. Power analysis
# spatiotemporal model varying the number of years and sites with data
# 
# 3. Performance on axis
# 
# I'm unsure if this is necessary for this paper but I could vary the detection rate as you suggest holding everything else constant in a spatial model. This would be relevant for other fish species, other taxa (stream salamanders), and *maybe* YOY vs adults.

if( Sys.getenv("USERNAME") %in% c("James.Thorson","xJames.Thorson") ) setwd( "C:/Users/James.Thorson/Desktop/Project_git/Trout_GRF/" )

# clear environment
rm(list = ls())
gc()

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

#######################
# Load data
#######################
load( "Data/White_River_Network.RData")
colnames(family)[1] = "child_name"
family = cbind( family, "child_b"=1:nrow(family) )

# Set Conditions
mean_N <- 50
n_years_vec <- c(4, 8, 10, 15, 20)
n_years <- max(n_years_vec)
sample_sites_vec <- c(25, 50, 100, 200, nrow(family))
p <- c(0.75, 0.75, 0.75)
theta <- 0.1
SD <- 0.05
rhot <- 0.5
SD_t <- 0.1
theta_st <- 0.1
SD_st <- 0.05
rho <- 0.8

n_sim <- 200

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
  network <- simST(family, theta = theta, SD = SD, rhot = rhot, SD_t = SD_t, theta_st = theta_st, SD_st = SD_st, mean_N = mean_N, n_years = n_years, rho = rho, gamma_j = gamma_j, X_ij=X_ij, p = p)
  # Check average sample AR of log-density at each site
  # This is only equal to rhot, or rho, when the other process has variance (i.e., SD_t or SD_st) fixed at zero, and is otherwise some kind of weighted average of the two (I think)
  ar1 = function(vec) var(vec[-length(vec)],vec[-1]) / var(vec)
  mean(apply(network$log_Npred_bt, MARGIN=1, FUN=ar1))
}

# Compile
if(FALSE) {
  dyn.unload(dynlib(paste0("Code/", Version)))
  file.remove( paste0("Code/", Version,c(".o",".dll")) )
}
compile( paste0("Code/", Version,".cpp") )

#---------- simulations --------------
dat <- data.frame(iter=integer(),
                  n_sites=integer(),
                  n_years=integer(),
                  spatialTF=integer(),
                  mean_N=numeric(),
                  min_N=numeric(),
                  max_N=numeric(),
                  mean_N_est=numeric(),
                  N_se=numeric(),
                  RMSE=numeric(),
                  theta=numeric(),
                  theta_hat=numeric(),
                  rhot=numeric(),
                  rhot_hat=numeric(),
                  sigmat=numeric(),
                  sigmat_hat=numeric(),
                  theta_st=numeric(),
                  theta_st_hat=numeric(),
                  rho_st=numeric(),
                  rho_st_hat=numeric(),
                  gamma_j=numeric(),
                  gamma_j_hat=numeric(),
                  stringsAsFactors=FALSE)

df_sims <- dat

counter <- 0
########## Set Up Parallel Processing ##########
library(foreach)
library(doParallel)

# set up parallel backend & make database connection available to all workers
nc <- min(c(detectCores()-1, 15)) #
cl <- makeCluster(nc, type = "PSOCK")
registerDoParallel(cl)

# # setup to write out to monitor progress
# logFile = paste0(data_dir, "/log_file.txt")
# logFile_Finish = paste0(data_dir, "/log_file_finish.txt")
# cat("Monitoring progress of prediction loop in parallel", file=logFile, append=FALSE, sep = "\n")
# cat("Monitoring the finish of each loop", file=logFile_Finish, append=FALSE, sep = "\n")

########## Run Parallel Loop ########## 
# start loop
df_sims <- foreach(i = 1:n_sim, 
                        .inorder=FALSE, 
                        .combine = rbind,
                        .packages=c("TMB",
                                    "dplyr",
                                    "minqa",
                                    "lubridate",
                                    "tidyr")
                        #.export = c("indexDeployments", "deriveMetrics") # shouldn't be needed after update package
                        #.export = c("derive_metrics_par")#,
                        #.export=ls(envir=globalenv(),
                        #          "indexDeployments")# shouldn't be needed after update package
) %dopar% {
  #for(i in 1:2) {
  source("Functions/Input_Functions.R")
  source("Functions/simOUGMRF.R")
  source("Functions/simST.R")
  source("Functions/runOUGMRF.R")
  source("Functions/summary_functions.R")
  
  # for(i in 1:n_sim) {
  # set conditions
  df_N <- matrix(NA, nrow(family), length(sample_sites_vec))
  df_coef <- list()
  df_sd <- list()
  df_rmse <- NA
  N_hat <- list()
  
  # simulate abundance and counts on network
  #set.seed(723750)
  network <- simST(family, theta = theta, SD = SD, rhot = rhot, SD_t = SD_t, theta_st = theta_st, SD_st = SD_st, mean_N = mean_N, n_years = n_years, rho = rho, gamma_j = gamma_j, X_ij=X_ij, p = p)
  str(network)
  summary(network$N_i)
  
  # thin to sample sites and years
  for(b in 1:length(sample_sites_vec)) {
    for(ti in 1:length(n_years_vec)) {
      for(s in 1:2) {
        counter <- counter + 1
        c_ip <- network$c_ip
        n <- length(network[["N_i"]])
        
        c_ip$t_i <- network$t_i
        c_ip$child_b <- rep(family$child_b, times = n_years)
        
        # Sample locations
        sites <- 1:max(family$child_b)
        remove_sites <- sample(sites, size = max(family$child_b) - sample_sites_vec[b], replace = FALSE) 
        c_ip[which(c_ip$child_b %in% remove_sites), 1:3] <- NA
        
        # Sample Years
        years <- 1:n_years_vec[ti]
        c_ip[which(!(c_ip$t_i %in% years)), 1:3] <- NA
        c_ip <- as.matrix(c_ip[ , 1:length(p)])
        
        #--------- Fit Model in TMB ----------
        start <- 1
        end <- 2
        Calc_lambda_ip <- rep(NA, length.out = nrow(network$c_ip))
        Calc_lambda_ip[start:end] <- 1
        Calc_lambda_ip[is.na(Calc_lambda_ip)] <- 0
        
        if(s == 1) {
          Options_vec = c("SpatialTF"=0, "TemporalTF"=1, "SpatiotemporalTF"=0, "DetectabilityTF"=1, "ObsModel"=1, "OverdispersedTF"=0, "abundTF"=0)
        }
        if(s == 2) {
          Options_vec = c("SpatialTF"=1, "TemporalTF"=1, "SpatiotemporalTF"=1, "DetectabilityTF"=1, "ObsModel"=1, "OverdispersedTF"=0, "abundTF"=0)
        }
        
        # need to make df, family, and c_ip agree
        df_counts <- as.data.frame(c_ip)
        df_counts$child_b <- rep(1:nrow(family), times = n_years)
        df <- dplyr::left_join(df_counts, family)
        
        # Make inputs
        Inputs <- makeInput(family = family, df = df, c_ip = c_ip, options = Options_vec, X = X_ij, t_i = network$t_i, version = Version, CalcSD_lambda_ip = Calc_lambda_ip)
        
        mod <- runOUGMRF(inputs = Inputs)
        
        #----------- summarize -----------
        
        if(class(mod) == "try-error") {
          dat[counter, "iter"] <- i
          dat[counter, "n_sites"] <- sample_sites_vec[b]
          dat[counter, "n_years"] <- n_years_vec[ti]
          dat[counter, "spatialTF"] <- s - 1
          dat[counter, "mean_N"] <- mean(network$N_i)
          dat[counter, "min_N"] <- min_N = min(network$N_i, na.rm = T)
          dat[counter, "max_N"] <- max_N = max(network$N_i, na.rm = T)
          dat[counter, "mean_N_est"] <- NA_real_
          dat[counter, "N_se"] <- NA_real_
          dat[counter, "RMSE"] <- NA_real_
          dat[counter, "theta"] <- theta
          dat[counter, "theta_hat"] <- NA_real_
          dat[counter, "rhot"] <- rhot
          dat[counter, "rhot_hat"] <- NA_real_
          dat[counter, "sigmat"] <- SD_t
          dat[counter, "sigmat_hat"] <- NA_real_
          dat[counter, "theta_st"] <- theta_st
          dat[counter, "theta_st_hat"] <- NA_real_
          dat[counter, "rho_st"] <- rho
          dat[counter, "rho_st_hat"] <- NA_real_
          dat[counter, "gamma_j"] <- gamma_j
          dat[counter, "gamma_j_hat"] <- NA_real_
          dat[counter, "converge"] <- FALSE
        } else {
          # check convergence
          converge <- FALSE
          try(converge <- mod$opt$convergence == 0 & !(mean(mod_out$SD$sd) == "NaN") & !any(is.na(mod_out$SD$sd)) & max(mod_out$SD$sd, na.rm = T) < 100)
          
          if(converge == FALSE) {
            dat[counter, "iter"] <- i
            dat[counter, "n_sites"] <- sample_sites_vec[b]
            dat[counter, "n_years"] <- n_years_vec[ti]
            dat[counter, "spatialTF"] <- s - 1
            dat[counter, "mean_N"] <- mean(network$N_i)
            dat[counter, "min_N"] <- min_N = min(network$N_i, na.rm = T)
            dat[counter, "max_N"] <- max_N = max(network$N_i, na.rm = T)
            dat[counter, "mean_N_est"] <- NA_real_
            dat[counter, "N_se"] <- NA_real_
            dat[counter, "RMSE"] <- NA_real_
            dat[counter, "theta"] <- theta
            dat[counter, "theta_hat"] <- NA_real_
            dat[counter, "rhot"] <- rhot
            dat[counter, "rhot_hat"] <- NA_real_
            dat[counter, "sigmat"] <- SD_t
            dat[counter, "sigmat_hat"] <- NA_real_
            dat[counter, "theta_st"] <- theta_st
            dat[counter, "theta_st_hat"] <- NA_real_
            dat[counter, "rho_st"] <- rho
            dat[counter, "rho_st_hat"] <- NA_real_
            dat[counter, "gamma_j"] <- gamma_j
            dat[counter, "gamma_j_hat"] <- NA_real_
            dat[counter, "converge"] <- FALSE
          } else {
          try(N_se <- mod$SD$sd[which(names(mod$SD$value) == "mean_N")])
          N_se <- ifelse(is.null(N_se), NA_real_, N_se)
          if(!is.na(N_se)) {
          N_se <- ifelse(N_se == "NaN", NA_real_, N_se)
          }
          df_N <- data.frame(N_i = network$N_i, N_hat = NA_real_)
          try(df_N <- data.frame(N_i = network$N_i, N_hat = mod$Report$N_ip[,1]))
          
          dat[counter, "iter"] <- i
          dat[counter, "n_sites"] <- sample_sites_vec[b]
          dat[counter, "n_years"] <- n_years_vec[ti]
          dat[counter, "spatialTF"] <- s - 1
          dat[counter, "mean_N"] <- mean(network$N_i)
          dat[counter, "min_N"] <- min_N = min(network$N_i, na.rm = T)
          dat[counter, "max_N"] <- max_N = max(network$N_i, na.rm = T)
          dat[counter, "mean_N_est"] <- mean(mod$Report$N_ip[ , 1])
          dat[counter, "N_se"] <- N_se
          dat[counter, "RMSE"] <- rmse(df_N$N_i - df_N$N_hat)
          dat[counter, "theta"] <- theta
          dat[counter, "theta_hat"] <- mod$Report$theta
          dat[counter, "rhot"] <- rhot
          dat[counter, "rhot_hat"] <- mod$Report$rhot
          dat[counter, "sigmat"] <- SD_t
          dat[counter, "sigmat_hat"] <- mod$Report$sigmat
          dat[counter, "theta_st"] <- theta_st
          dat[counter, "theta_st_hat"] <- mod$Report$theta_st
          dat[counter, "rho_st"] <- rho
          dat[counter, "rho_st_hat"] <- mod$Report$rho_st
          dat[counter, "gamma_j"] <- gamma_j
          dat[counter, "gamma_j_hat"] <- mod$Report$gamma_j
          dat[counter, "converge"] <- FALSE
          
          }
        }
        #         mean(network$N_i)
        #         SD_means <- data.frame(param = names(mod$SD$value), 
        #                                est = as.numeric(mod$SD$value), 
        #                                sd = as.numeric(mod$SD$sd), stringsAsFactors = F)
        
        #-------- save sim iter output -------
        save(network, mod, file = paste0("Output/Power_Sim/Data/sim_", i, "_st_", s-1, "_sites_", sample_sites_vec[b], "_years_", n_years_vec[ti], ".RData")) 
      } # end spatial TF loop
    } # end site loop
  } # end year loop
  write.csv(dat, file = paste0("Output/Power_Sim/Data/summary_sim", i, ".csv"), row.names = FALSE)
  return(dat)
} # end sim iter
stopCluster(cl)
closeAllConnections()

save(df_sims, file = "Output/Power_Sim/STsim_Results.RData")
