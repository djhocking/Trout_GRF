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
source("Functions/runOUGMRF.R")
source("Functions/summary_functions.R")

dir_out <- "Output"

#######################
# Load data
#######################
load( "Data/White_River_Network.RData")
colnames(family)[1] = "child_name"
family = cbind( family, "child_b"=1:nrow(family) )

if( any(family$dist_b < (max(family$dist_b,na.rm=TRUE)/1e4),na.rm=TRUE) ) {
  family$dist_b = ifelse( family$dist_b<(max(family$dist_b,na.rm=TRUE)/1e4), (max(family$dist_b,na.rm=TRUE)/1e4), family$dist_b)
  warning("Negative distances fixed to lower bound")
}

Version = "OU_GMRF_v1h"

# Compile
if(FALSE) {
  dyn.unload(dynlib(paste0("Code/", Version)))
  file.remove( paste0("Code/", Version,c(".o",".dll")) )
}
compile( paste0("Code/", Version,".cpp") )

###############################
# vary theta, SD, spatial
###############################
n_sim <- 200
# vary theta, SD, sigmaIID, log_mean, and sample_pct for spatial/non-spatial
spatial_cor <- c("extremely high", "high", "medium", "low", "very low")
theta_vec <- c(0.5, 1, 3, 5)
SD_vec <- c(0.1, 0.2, 0.3, 0.4) # 0.01, 0.1, 0.25 work with theta = c(0.5, 1, 5) and mean_N=c(5,10,50)
sigmaIID_vec <- c(0)
mean_N <- c(10)
spatial_vec <- c(TRUE, FALSE)
sample_pct_vec <- c(1) #

# Detection probability for removal sampling
p <- c(0.75, 0.75, 0.75)

# Covariates
gamma_j <- c(0.5) # doesn't work if I use more than 1 coef
X_ij <- matrix(rnorm(nrow(family), 0, 1), nrow(family), length(gamma_j))

options_list <- list(c("SpatialTF"=0, "TemporalTF"=0, "SpatiotemporalTF"=0, "DetectabilityTF"=1, "ObsModel"=1, "OverdispersedTF"=0, "abundTF"=1),
                     c("SpatialTF"=1, "TemporalTF"=0, "SpatiotemporalTF"=0, "DetectabilityTF"=1, "ObsModel"=1, "OverdispersedTF"=0, "abundTF"=1))

#######################
# Simulate GMRF following O-U process
#######################

if(!file.exists("Output/Sim_Spatial/Figures/")) {
  dir.create("Output/Sim_Spatial/Figures/", recursive = TRUE)
}
if(!file.exists("Output/Sim_Spatial/Data/")) {
  dir.create("Output/Sim_Spatial/Data/", recursive = TRUE)
}

I <- length(theta_vec)
J <- length(SD_vec)
K <- length(mean_N)
L <- length(options_list)
M <- length(spatial_vec)
df_N <- matrix(NA, nrow(family), I*J*K*L + K*L)
df_coef <- list()
df_sd <- list()
df_rmse <- NA

counter <- 0
mod2 <- list()

library(foreach)
library(doParallel)

# set up parallel backend & make database connection available to all workers
nc <- min(c(detectCores()-1, 12)) #
cl <- makeCluster(nc, type = "PSOCK")
registerDoParallel(cl)

# # setup to write out to monitor progress
#logFile = "Output/Sim_Spatial/Data/sims_table.csv")
# logFile_Finish = paste0(data_dir, "/log_file_finish.txt")
# cat("Monitoring progress of prediction loop in parallel", file=logFile, append=FALSE, sep = "\n")
# cat("Monitoring the finish of each loop", file=logFile_Finish, append=FALSE, sep = "\n")

########## Run Parallel Loop ########## 
# start loop
df_sims <- foreach(sim = 1:n_sim, 
                   .inorder=FALSE, 
                   .combine = rbind,
                   .packages=c("TMB",
                               "dplyr",
                               "minqa",
                               "lubridate",
                               "tidyr")
) %dopar% {
  
  source("Functions/Input_Functions.R")
  source("Functions/simOUGMRF.R")
  source("Functions/runOUGMRF.R")
  source("Functions/summary_functions.R")
  
  table_full <- NULL
  table_full2 <- NULL
  table_combined <- NULL
  
  for(i in 1:length(theta_vec)) { 
    for(j in 1:length(SD_vec)) {
      for(k in 1:length(mean_N)) {
        m <- 1
        
        network <- simOUGMRF(family = family, 
                             theta = theta_vec[i],
                             SD = SD_vec[j],
                             mean_N = mean_N[k],
                             gamma = gamma_j,
                             X_ij = X_ij,
                             p = p,
                             sample_pct = sample_pct_vec[1],
                             spatial = spatial_vec[m]
        )
        
        for(l in 1:length(options_list)) {
          #counter <- counter + 1
          Options_vec = options_list[[l]]
          
          start <- 1
          end <- 2
          Calc_lambda_ip <- rep(NA, length.out = nrow(network$c_ip))
          Calc_lambda_ip[start:end] <- 1
          Calc_lambda_ip[is.na(Calc_lambda_ip)] <- 0
          
          # Make inputs
          Inputs <- makeInput(family = family, c_ip = network$c_ip, options = Options_vec, X = X_ij, t_i = network$t_i, version = Version, CalcSD_lambda_ip = Calc_lambda_ip)
          
          # run model
          mod_out <- NULL
          optim_err <- tryCatch(
            mod_out <- runOUGMRF(inputs = Inputs),
            error = function(e) e
          )
          
          error <- TRUE
          if(!inherits(optim_err, "try-error")) { # if not error check that model output produced
            if(!is.null(mod_out)) { # if output produced, check SD produced
              error <- FALSE
              converge <- FALSE
              if(!is.null(mod_out$SD)) { # if SD produced, check sd produced
                if(!is.null(mod_out$SD$sd)) { # if sd produced, check that converged
                  converge <- ifelse(mod_out$opt$convergence == 0 & !any(is.na(mod_out$SD$sd)), TRUE, FALSE)
                }
              }
            }
          }
          
          if(error == TRUE | converge == FALSE) {
            table_i <- data.frame(sim = sim,
                                  spatial = spatial_vec[m],
                                  sp_mod = Options_vec[["SpatialTF"]],
                                  theta = theta_vec[i], 
                                  SD_ou = SD_vec[j], 
                                  mean_N_set = mean_N[k], 
                                  mean_N = NA_real_,
                                  mean_N_hat = NA_real_,
                                  min_N = NA_real_,
                                  max_N = NA_real_,
                                  N_se = NA_real_,
                                  AIC = NA_real_,
                                  error = error,
                                  converge = converge, 
                                  rmse = NA_real_, 
                                  resid_mean = NA_real_, 
                                  theta_hat = NA_real_, 
                                  SD_ou_hat = NA_real_, 
                                  coef_hat = NA_real_, 
                                  coef_sd = NA_real_,
                                  coef_true = gamma_j,
                                  stringsAsFactors = FALSE)
          } else { # if model output produced with no error and no convergence problems
            # calculate expected conditional N_i
            N_hat <- data.frame(N_hat = mod_out$Report$N_ip[,1], lambda_hat = mod_out$Report$lambda_ip[ , 1], pass_1 = network$c_ip[ , 1])
            N_hat <- N_hat %>%
              dplyr::mutate(N_hat = ifelse(is.na(pass_1), lambda_hat, N_hat))
            
            # Get SE of mean N across the network
            N_se <- NULL
            try(N_se <- mod_out$SD$sd[which(names(mod_out$SD$value) == "mean_N")])
            N_se <- ifelse(is.null(N_se), NA_real_, N_se)
            if(!is.na(N_se)) {
              N_se <- ifelse(N_se == "NaN", NA_real_, N_se)
            }
            
            # make output table for specific conditions and simulation
            table_i <- data.frame(sim = sim,
                                  spatial = spatial_vec[m],
                                  sp_mod = Options_vec[["SpatialTF"]],
                                  theta = theta_vec[i], 
                                  SD_ou = SD_vec[j], 
                                  mean_N_set = mean_N[k], 
                                  mean_N = mean(network$N_i, na.rm = T),
                                  mean_N_hat = mean(mod_out$Report$N_ip[ , 1]),
                                  min_N = min(network$N_i, na.rm = T),
                                  max_N = max(network$N_i, na.rm = T), 
                                  N_se = N_se,
                                  error = error,
                                  converge = converge, 
                                  AIC = mod_out$opt$AIC,
                                  rmse = rmse(network$N_i - N_hat$N_hat), 
                                  resid_mean = mean(network$N_i - N_hat$N_hat), 
                                  theta_hat = mod_out$Report$theta, 
                                  SD_ou_hat = mod_out$Report$SDinput, 
                                  coef_hat = mod_out$Report$gamma_j, 
                                  coef_sd = mod_out$SD$sd[which(names(mod_out$SD$value) == "gamma_j")],
                                  coef_true = gamma_j,
                                  stringsAsFactors = FALSE)
          }
          
          # apend to existing table across conditions for entire simulation iteration
          if(is.null(table_full)) {
            table_full <- table_i
          } else {
            table_full <- dplyr::bind_rows(table_full, table_i)
          }
          save(network, mod_out, file = paste0("Output/Sim_Spatial/Data/sim_", sim, "_theta_", theta_vec[i], "_SD_", SD_vec[j], "_spMod_", options_list[[l]][1], "_meanN_", mean_N[k], "_spatial_data_", m, ".RData")) 
        } # end options_list (spatial on/off) loop
        } # end mean_N loop
      } # end SD_vec loop
    } # end theta_vec loop
    
    # repeat the process for non-spatial data (don't need to vary theta and SD)
    m <- 2
    network2 <- simOUGMRF(family = family, 
                         theta = theta_vec[1],
                         SD = SD_vec[1],
                         mean_N = mean_N[1],
                         gamma = gamma_j,
                         X_ij = X_ij,
                         p = p,
                         sample_pct = sample_pct_vec[1],
                         spatial = spatial_vec[m]
    )
    
    for(l in 1:length(options_list)) {
      #counter <- counter + 1
      Options_vec = options_list[[l]]
      
      start <- 1
      end <- 2
      Calc_lambda_ip <- rep(NA, length.out = nrow(network2$c_ip))
      Calc_lambda_ip[start:end] <- 1
      Calc_lambda_ip[is.na(Calc_lambda_ip)] <- 0
      
      # Make inputs
      Inputs <- makeInput(family = family, c_ip = network2$c_ip, options = Options_vec, X = X_ij, t_i = network2$t_i, version = Version, CalcSD_lambda_ip = Calc_lambda_ip)
      
      # run model
      mod_out2 <- NULL
      optim_err <- tryCatch(
        mod_out2 <- runOUGMRF(inputs = Inputs),
        error = function(e) e
      )
      
      error <- TRUE
      if(!inherits(optim_err, "try-error")) { # if not error check that model output produced
        if(!is.null(mod_out2)) { # if output produced, check SD produced
          error <- FALSE
          converge <- FALSE
          if(!is.null(mod_out2$SD)) { # if SD produced, check sd produced
            if(!is.null(mod_out2$SD$sd)) { # if sd produced, check that converged
              converge <- ifelse(mod_out2$opt$convergence == 0 & !any(is.na(mod_out2$SD$sd)), TRUE, FALSE)
            }
          }
        }
      }
      
      if(error == TRUE | converge == FALSE) {
        table2_i <- data.frame(sim = sim,
                              spatial = spatial_vec[m],
                              sp_mod = Options_vec[["SpatialTF"]],
                              theta = NA_real_, 
                              SD_ou = NA_real_, 
                              mean_N_set = mean_N[1], 
                              mean_N = NA_real_,
                              mean_N_hat = NA_real_,
                              min_N = NA_real_,
                              max_N = NA_real_,
                              N_se = NA_real_,
                              AIC = NA_real_,
                              error = error,
                              converge = converge, 
                              rmse = NA_real_, 
                              resid_mean = NA_real_, 
                              theta_hat = NA_real_, 
                              SD_ou_hat = NA_real_, 
                              coef_hat = NA_real_, 
                              coef_sd = NA_real_,
                              coef_true = gamma_j,
                              stringsAsFactors = FALSE)
      } else { # if model output produced with no error and no convergence problems
        # calculate expected conditional N_i
        N_hat <- data.frame(N_hat = mod_out2$Report$N_ip[,1], lambda_hat = mod_out2$Report$lambda_ip[ , 1], pass_1 = network2$c_ip[ , 1])
        N_hat <- N_hat %>%
          dplyr::mutate(N_hat = ifelse(is.na(pass_1), lambda_hat, N_hat))
        
        # Get SE of mean N across the network2
        N_se <- NULL
        try(N_se <- mod_out2$SD$sd[which(names(mod_out2$SD$value) == "mean_N")])
        N_se <- ifelse(is.null(N_se), NA_real_, N_se)
        if(!is.na(N_se)) {
          N_se <- ifelse(N_se == "NaN", NA_real_, N_se)
        }
        
        # make output table for specific conditions and simulation
        table2_i <- data.frame(sim = sim,
                              spatial = spatial_vec[m],
                              sp_mod = Options_vec[["SpatialTF"]],
                              theta = NA_real_, 
                              SD_ou = NA_real_, 
                              mean_N_set = mean_N[1], 
                              mean_N = mean(network2$N_i, na.rm = T),
                              mean_N_hat = mean(mod_out2$Report$N_ip[ , 1]),
                              min_N = min(network2$N_i, na.rm = T),
                              max_N = max(network2$N_i, na.rm = T), 
                              N_se = N_se,
                              error = error,
                              converge = converge, 
                              AIC = mod_out2$opt$AIC,
                              rmse = rmse(network2$N_i - N_hat$N_hat), 
                              resid_mean = mean(network2$N_i - N_hat$N_hat), 
                              theta_hat = mod_out2$Report$theta, 
                              SD_ou_hat = mod_out2$Report$SDinput, 
                              coef_hat = mod_out2$Report$gamma_j, 
                              coef_sd = mod_out2$SD$sd[which(names(mod_out2$SD$value) == "gamma_j")],
                              coef_true = gamma_j,
                              stringsAsFactors = FALSE)
      }
      
      # apend to existing table across conditions for entire simulation iteration
      if(is.null(table_full2)) {
        table_full2 <- table2_i
      } else {
        table_full2 <- dplyr::bind_rows(table_full2, table2_i)
      }
      save(network2, mod_out2, file = paste0("Output/Sim_Spatial/Data/sim_", sim, "_theta_", theta_vec[1], "_SD_", SD_vec[1], "_spMod_", options_list[[l]][1], "_meanN_", mean_N[1], "_spatial_data_", m, ".RData")) 
    } # end spatial model options loop
  table_combined <- dplyr::bind_rows(table_full, table_full2)
  return(table_combined)
} # end foreach dopar loop
stopCluster(cl)
closeAllConnections()

save(df_sims, file = "Output/Sim_Spatial/Sim_Results.RData")
write.csv(df_sims, file = "Output/Sim_Spatial/Sim_Results.csv", row.names = FALSE)   
          
          
          
          
             
              