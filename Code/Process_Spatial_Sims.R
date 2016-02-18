
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

n_sim <- 200
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

table_full <- NULL
for(sim in 1:n_sim) {
  for(i in 1:length(theta_vec)) { 
    for(j in 1:length(SD_vec)) { # theta must be larger than SD?
      for(k in 1:length(mean_N)) {
        m <- 1
          for(l in 1:length(options_list)) {
            mod_out <- NULL
            file_input <- paste0("Output/Sim_Spatial/Data/sim_", sim, "_theta_", theta_vec[i], "_SD_", SD_vec[j], "_spatialTF_", options_list[[l]][1], "_meanN_", mean_N[k], ".RData")
            if(file.exists(file_input)) {
            load(file = paste0("Output/Sim_Spatial/Data/sim_", sim, "_theta_", theta_vec[i], "_SD_", SD_vec[j], "_spatialTF_", options_list[[l]][1], "_meanN_", mean_N[k], ".RData")) 
            }
            
            #counter <- counter + 1
            Options_vec = options_list[[l]]
            
            if(is.null(mod_out)) {
              #mod2[[counter]] <- NA
              table_i <- data.frame(spatial = spatial_vec[m],
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
                                    converge = FALSE, 
                                    rmse = NA_real_, 
                                    resid_mean = NA_real_, 
                                    theta_hat = NA_real_, 
                                    SD_ou_hat = NA_real_, 
                                    coef_hat = NA_real_, 
                                    coef_sd = NA_real_,
                                    coef_true = gamma_j)
            } else {          
              # check convergence
              converge <- FALSE
              if(!is.null(mod_out$SD$sd)) {
              try(converge <- mod_out$opt$convergence == 0 & !any(is.na(mod_out$SD$sd)))
              }
              
              # organize table of output
              if(converge == FALSE) {
                table_i <- data.frame(spatial = spatial_vec[m],
                                      sp_mod = Options_vec[["SpatialTF"]],
                                      theta = theta_vec[i], 
                                      SD_ou = SD_vec[j], 
                                      mean_N_set = mean_N[k], 
                                      mean_N = mean(network$N_i, na.rm = T),
                                      mean_N_hat = NA_real_,
                                      min_N = min(network$N_i, na.rm = T),
                                      max_N = max(network$N_i, na.rm = T), 
                                      N_se = NA_real_,
                                      converge = converge, 
                                      AIC = NA_real_,
                                      rmse = NA_real_, 
                                      resid_mean = NA_real_, 
                                      theta_hat = NA_real_, 
                                      SD_ou_hat = NA_real_, 
                                      coef_hat = NA_real_, 
                                      coef_sd = NA_real_,
                                      coef_true = gamma_j)
              } else {
                # extract predicted abundance
                N_hat <- data.frame(N_hat = mod_out$Report$N_ip[,1], lambda_hat = mod_out$Report$lambda_ip[ , 1], pass_1 = network$c_ip[ , 1])
                N_hat <- N_hat %>%
                  dplyr::mutate(N_hat = ifelse(is.na(pass_1), lambda_hat, N_hat))
                
                N_se <- NULL
                try(N_se <- mod_out$SD$sd[which(names(mod_out$SD$value) == "mean_N")])
                N_se <- ifelse(is.null(N_se), NA_real_, N_se)
                if(!is.na(N_se)) {
                  N_se <- ifelse(N_se == "NaN", NA_real_, N_se)
                }
                
                table_i <- data.frame(spatial = spatial_vec[m],
                                      sp_mod = Options_vec[["SpatialTF"]],
                                      theta = theta_vec[i], 
                                      SD_ou = SD_vec[j], 
                                      mean_N_set = mean_N[k], 
                                      mean_N = mean(network$N_i, na.rm = T),
                                      mean_N_hat = mean(mod_out$Report$N_ip[ , 1]),
                                      min_N = min(network$N_i, na.rm = T),
                                      max_N = max(network$N_i, na.rm = T), 
                                      N_se = N_se,
                                      converge = converge, 
                                      AIC = mod_out$opt$AIC,
                                      rmse = rmse(network$N_i - N_hat$N_hat), 
                                      resid_mean = mean(network$N_i - N_hat$N_hat), 
                                      theta_hat = mod_out$Report$theta, 
                                      SD_ou_hat = mod_out$Report$SDinput, 
                                      coef_hat = mod_out$Report$gamma_j, 
                                      coef_sd = mod_out$SD$sd[which(names(mod_out$SD$value) == "gamma_j")],
                                      coef_true = gamma_j)
              }
            }
            table_i$sim <- sim
            
            # apend to existing table
            if(is.null(table_full)) {
              table_full <- table_i
            } else {
              table_full <- dplyr::bind_rows(table_full, table_i)
            }
            rm(network)
            rm(mod_out)
          }
        }
      }
    }
  }




# Need to redo so don't have separate value for every m

df_sims <- table_full %>%
  dplyr::mutate(theta = ifelse(spatial == FALSE, NA_real_, theta),
                SD_ou = ifelse(spatial == FALSE, NA_real_, SD_ou)) %>%
  dplyr::distinct() %>%
  dplyr::group_by(spatial, sp_mod, theta, SD_ou, mean_N_set, coef_true)

df_sims_spatialT <- dplyr::filter(df_sims, spatial == TRUE)
df_sims_spatialF <- dplyr::filter(df_sims, spatial == FALSE)

rows <- 1:(nrow(df_sims_spatialF)-3)
df_sims_spatialF$shift2 <- c(0, 0, 0, as.numeric(df_sims_spatialF[rows, ]$sim))
df_sims_spatialF <- df_sims_spatialF %>%
  dplyr::mutate(extra = ifelse(sim == shift2, TRUE, FALSE)) %>%
  dplyr::filter(extra == FALSE) %>%
  dplyr::select(-shift2, -extra)

df_sims <- bind_rows(df_sims_spatialT, df_sims_spatialF)

save(df_sims, file = "Output/Sim_Spatial/Sim_Results.RData")
write.csv(df_sims, file = "Output/Sim_Spatial/Sim_Results.csv", row.names = FALSE)

