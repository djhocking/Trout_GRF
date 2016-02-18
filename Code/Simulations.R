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

# 1. Proof it works
# Spatial model varying theta and sigma. I think it makes sense to compare these to a non-spatial model (sorry a sort of 3rd axis). Do you think we need to vary sigma since it affects the spatial correlation but as a constant and not in relation to distance? I don't think we need to test temporal or spatiotemporal components here.
# 
# 2. Power analysis
# spatiotemporal model varying the number of years and sites with data
# 
# 3. Performance on axis
# I'm unsure if this is necessary for this paper but I could vary the detection rate as you suggest holding everything else constant in a spatial model. This would be relevant for other fish species, other taxa (stream salamanders), and *maybe* YOY vs adults.

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


n_sim <- 200
counter <- 0
mod2 <- list()

library(foreach)
library(doParallel)

# set up parallel backend & make database connection available to all workers
nc <- min(c(detectCores()-1, 12)) #
cl <- makeCluster(nc, type = "PSOCK")
registerDoParallel(cl)

# # setup to write out to monitor progress
# logFile = paste0(data_dir, "/log_file.txt")
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
  #counter <- 0
  
#   saveRDS(list(spatial_cor = spatial_cor,
#                theta_vec = theta_vec,
#                SD_vec = SD_vec,
#                sigmaIID_vec = sigmaIID_vec,
#                mean_N = mean_N,
#                spatial_vec = spatial_vec,
#                p = p, 
#                gamma_j = gamma_j),
#           file = paste0("Output/Sim_Spatial/Data/spatial_sims_input_iter_", sim, ".RData"
#           ))
  
  #set.seed(1987654)
  # adjust so don't have a new random sample when comparing spatial and non-spatial models
  
  for(i in 1:length(theta_vec)) { 
    for(j in 1:length(SD_vec)) { # theta must be larger than SD?
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
          
          if(inherits(optim_err, "try-error")) {
            #mod2[[counter]] <- NA
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
            # compile all model output
            #mod2[[counter]] <- mod_out
            
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
              
              # compile output of interest
#               df_N[ , counter] <- N_hat$N_hat
#               df_coef[[counter]] <- as.numeric(mod_out$SD$value)
#               df_sd[[counter]] <- as.numeric(mod_out$SD$sd)
#               df_rmse[counter] <- rmse(network$N_i - N_hat$N_hat)
              
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
          # table_full[counter, "sim"] <- sim
          
          #           df_N_plot <- data.frame(N_hat = N_hat$N_hat, N_i = network$N_i, obsTF = ifelse(is.na(network$c_ip[,1]), FALSE, TRUE))
          #           g <- ggplot(df_N_plot, aes(N_i, N_hat, colour = obsTF)) + geom_point() + geom_abline(intercept = 0, slope = 1, colour = "blue") + theme_bw() + ggtitle(paste0("spatial=", as.integer(spatial_vec[m]), "; sp_mod=", Options_vec[["SpatialTF"]], "; theta=", theta_vec[i], "; SD_ou=", SD_vec[j], "; mean_N=", mean_N[k]))
          #           ggsave(filename = paste0("Output/Sim_Spatial/Figures/model_", counter[i], ".pdf"), plot = g)
          
         save(network, mod_out, file = paste0("Output/Sim_Spatial/Data/sim_", sim, "_theta_", theta_vec[i], "_SD_", SD_vec[j], "_spatialTF_", options_list[[l]][1], "_meanN_", mean_N[k], ".RData")) 
        }
      }
    }
  }
  
  for(k in 1:length(mean_N)) {
    m <- 2
    
    network <- simOUGMRF(family = family, 
                         theta = theta_vec[1],
                         SD = SD_vec[1],
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
      optim_err <- tryCatch(
        mod_out <- runOUGMRF(inputs = Inputs),
        error = function(e) e
      )
      
      if(inherits(optim_err, "try-error")) {
        #mod2[[counter]] <- NA
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
        # compile all model output
        #mod2[[counter]] <- mod_out
        
        # check convergence
        converge <- FALSE
        try(converge <- mod_out$opt$convergence == 0 & !any(is.na(mod_out$SD$sd)))
        
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
          
          # compile output of interest
#           df_N[ , counter] <- N_hat$N_hat
#           df_coef[[counter]] <- as.numeric(mod_out$SD$value)
#           df_sd[[counter]] <- as.numeric(mod_out$SD$sd)
#           df_rmse[counter] <- rmse(network$N_i - N_hat$N_hat)
          
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
      #table_full$sim <- sim
      
      #           df_N_plot <- data.frame(N_hat = N_hat$N_hat, N_i = network$N_i, obsTF = ifelse(is.na(network$c_ip[,1]), FALSE, TRUE))
      #           g <- ggplot(df_N_plot, aes(N_i, N_hat, colour = obsTF)) + geom_point() + geom_abline(intercept = 0, slope = 1, colour = "blue") + theme_bw() + ggtitle(paste0("spatial=", as.integer(spatial_vec[m]), "; sp_mod=", Options_vec[["SpatialTF"]], "; theta=", theta_vec[i], "; SD_ou=", SD_vec[j], "; mean_N=", mean_N[k]))
      #           ggsave(filename = paste0("Output/Sim_Spatial/Figures/model_", counter[i], ".pdf"), plot = g)
      
      save(network, mod_out, file = paste0("Output/Sim_Spatial/Data/sim_", sim, "_theta_", theta_vec[i], "_SD_", SD_vec[j], "_spatialTF_", options_list[[l]][1], "_meanN_", mean_N[k], ".RData")) 
    }
  }
  table_full <- as.data.frame(table_full)
  #write.csv(table_full, file = paste0("Output/Sim_Spatial/Data/summary_sim", i, ".csv"), row.names = FALSE)
  return(table_full)
} # end sim iter
stopCluster(cl)
closeAllConnections()

save(df_sims, file = "Output/Sim_Spatial/Sim_Results.RData")
write.csv(df_sims, file = "Output/Sim_Spatial/Sim_Results.csv", row.names = FALSE)


















#####################


for(k in 1:length(mean_N)) {
  for(l in 1:length(options_list)) {
    m <- 2
    counter <- counter + 1
    network <- simOUGMRF(family = family, 
                         theta = theta_vec[1],
                         SD = SD_vec[1],
                         mean_N = mean_N[k],
                         gamma = gamma_j,
                         X_ij = X_ij,
                         p = p,
                         sample_pct = sample_pct_vec[1],
                         spatial = spatial_vec[m]
    )
    
    Options_vec = options_list[[l]]
    
    # Make inputs
    Inputs <- makeInput(family = family, c_ip = network$c_ip, options = Options_vec, X = X_ij, t_i = network$t_i, version = Version, CalcSD_lambda_ip = Calc_lambda_ip)
    
    # run model
    optim_err <- tryCatch(
      mod_out <- runOUGMRF(inputs = Inputs),
      error = function(e) e
    )
    
    if(inherits(optim_err, "try-error")) {
      mod2[[counter]] <- NA
      table_i <- data.frame(spatial = spatial_vec[m],
                            sp_mod = Options_vec[["SpatialTF"]],
                            theta = theta_vec[1], 
                            SD_ou = SD_vec[1], 
                            mean_N = mean_N[k],  
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
      # compile all model output
      mod2[[counter]] <- mod_out
      
      # extract predicted abundance
      N_hat <- data.frame(N_hat = mod_out$Report$N_ip[,1], lambda_hat = mod_out$Report$lambda_ip[ , 1], pass_1 = network$c_ip[ , 1])
      N_hat <- N_hat %>%
        dplyr::mutate(N_hat = ifelse(is.na(pass_1), lambda_hat, N_hat))
      
      # compile output of interest
      df_N[ , counter] <- N_hat$N_hat
      df_coef[[counter]] <- as.numeric(mod_out$SD$value)
      df_sd[[counter]] <- as.numeric(mod_out$SD$sd)
      df_rmse[counter] <- rmse(network$N_i - N_hat$N_hat)
      
      # check convergence
      converge <- mod_out$opt$convergence == 0 & !(mean(mod_out$SD$sd) == "NaN") & max(mod_out$SD$sd, na.rm = T) < 100
      
      # organize table of output
      if(converge == FALSE) {
        table_i <- data.frame(spatial = spatial_vec[m],
                              sp_mod = Options_vec[["SpatialTF"]],
                              theta = theta_vec[1], 
                              SD_ou = SD_vec[1], 
                              mean_N = mean_N[k],  
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
        table_i <- data.frame(spatial = spatial_vec[m],
                              sp_mod = Options_vec[["SpatialTF"]],
                              theta = theta_vec[1], 
                              SD_ou = SD_vec[1], 
                              mean_N = mean_N[k], 
                              converge = converge, 
                              AIC = mod_out$opt$AIC,
                              rmse = df_rmse[counter], 
                              resid_mean = mean(network$N_i - N_hat$N_hat), 
                              theta_hat = mod_out$Report$theta, 
                              SD_ou_hat = mod_out$Report$SDinput, 
                              coef_hat = mod_out$Report$gamma_j, 
                              coef_sd = mod_out$SD$sd[which(names(mod_out$SD$value) == "gamma_j")],
                              coef_true = gamma_j)
      }
    }
    
    # apend to existing table
    if(counter == 1) {
      table_full <- table_i
    } else {
      table_full <- dplyr::bind_rows(table_full, table_i)
    }
    
    df_N_plot <- data.frame(N_hat = N_hat$N_hat, N_i = network$N_i, obsTF = ifelse(is.na(network$c_ip[,1]), FALSE, TRUE))
    g <- ggplot(df_N_plot, aes(N_i, N_hat)) + geom_point() + geom_abline(aes(0,1), colour = "blue") + theme_bw() + ggtitle(paste0("spatial=", as.integer(spatial_vec[m]), "; sp_mod=", Options_vec[["SpatialTF"]], "; theta=", theta_vec[i], "; SD_ou=", SD_vec[j], "; mean_N=", mean_N[k]))
    ggsave(filename = paste0("Output/model_", counter[i], ".pdf"), plot = g)
  }
}
}
}
########################

table_full$model <- 1:nrow(table_full)


# format full table
table_full <- table_full %>%
  dplyr::mutate(theta_hat = ifelse(sp_mod == 1, theta_hat, NA_real_),
                SD_ou_hat = ifelse(sp_mod == 1, SD_ou_hat, NA_real_),
                theta = ifelse(spatial == TRUE, theta, NA_real_),
                SD_ou = ifelse(spatial == TRUE, SD_ou, NA_real_)) %>%
  dplyr::arrange(-spatial)
#data.frame(unclass(table_full))
format(table_full, digits = 3)

table_full %>%
  dplyr::group_by(spatial, sp_mod) %>%
  dplyr::summarise(mean_coef = mean(coef_hat, na.rm = T),
                   sd_coef = mean(coef_sd, na.rm = T))

table_full %>%
  dplyr::group_by(spatial, sp_mod) %>%
  dplyr::filter(rmse < 100) %>% # seems to be convergence problem not picked up
  dplyr::summarise(mean_resid = mean(resid_mean, na.rm = T),
                   mean_rmse = mean(rmse, na.rm = T))

saveRDS(table_full, file = file.path(dir_out, "Sim_Table.RData"))

save.image(file = file.path(dir_out, "Sim_Results_Full.RData"))



# relation of sigma_b and rho_b

df_sigma <- data.frame(SD_b_hat = mod2[[35]]$Report$SDinput_b, 
                       N_i = mod2[[35]]$Report$N_ip[, 1], 
                       dist = family$dist_b,
                       rho_b = mod2[[35]]$Report$rho_b)

ggplot(df_sigma, aes(dist, SD_b_hat)) + geom_point()

ggplot(df_sigma, aes(dist, rho_b)) + geom_point()

ggplot(df_sigma, aes(SD_b_hat, rho_b)) + geom_point()

# ability to recover SDinput_b - since combo of SDinput and theta which are not recovered well
df_sigma$SD_b <- ((mod2[[35]]$Report$SDinput^2)/(2*mod2[[35]]$Report$theta) * (1-exp(-2*mod2[[35]]$Report$theta*family$dist_b)))^0.5

SD <- 0.25
theta <- 0.5
df_sigma$SD_b_input <- ((SD^2)/(2*theta) * (1-exp(-2*theta*family$dist_b)))^0.5

ggplot(df_sigma, aes(SD_b, SD_b_hat)) + geom_point()

ggplot(df_sigma, aes(SD_b_input, SD_b_hat)) + geom_abline(intercept = 0, slope = 1, colour = "blue") + geom_point()




n_iters <- 100
theta_hat <- rep(NA, length.out = n_iters)
SD_hat <- rep(NA, length.out = n_iters)
for(i in 1:n_iters) {
  network <- simOUGMRF(family = family, 
                       theta = 0.2, # adult estimate from sp + t model
                       SD = 0.1, # adult estimate from sp + t model
                       mean_N = 40,
                       gamma = gamma_j,
                       X_ij = X_ij,
                       p = p,
                       sample_pct = sample_pct_vec[1],
                       spatial = TRUE
  )
  
  Options_vec = c("SpatialTF"=1, "TemporalTF"=0, "SpatiotemporalTF"=0, "DetectabilityTF"=1, "ObsModel"=1, "OverdispersedTF"=0, "abundTF"=1)
  
  Calc_lambda_ip <- rep(0, length.out = nrow(network$c_ip))
  
  # Make inputs
  Inputs <- makeInput(family = family, c_ip = network$c_ip, options = Options_vec, X = X_ij, t_i = network$t_i, version = Version, CalcSD_lambda_ip = Calc_lambda_ip)
  
  try({
    # run model
    mod_out <- runOUGMRF(inputs = Inputs)
    
    #str(mod_out$Report)
    
    theta_hat[i] <- mod_out$Report$theta
    SD_hat[i] <- mod_out$Report$SDinput
  })
  
}
df_reps <- data.frame(theta_hat, SD_hat)
df_reps[which(df_reps$theta > 100), ] <- NA
format(df_reps, digits = 2, scientific = FALSE)

summary(df_reps)



# what if the distance is on a different scale?
family2 <- family
faimly2$dist <- family$dist/10
n_iters <- 100
theta_hat <- rep(NA, length.out = n_iters)
SD_hat <- rep(NA, length.out = n_iters)
for(i in 1:n_iters) {
  network <- simOUGMRF(family = family2, 
                       theta = 0.2, # adult estimate from sp + t model
                       SD = 0.1, # adult estimate from sp + t model
                       mean_N = 40,
                       gamma = gamma_j,
                       X_ij = X_ij,
                       p = p,
                       sample_pct = sample_pct_vec[1],
                       spatial = TRUE
  )
  
  Options_vec = c("SpatialTF"=1, "TemporalTF"=0, "SpatiotemporalTF"=0, "DetectabilityTF"=1, "ObsModel"=1, "OverdispersedTF"=0, "abundTF"=1)
  
  Calc_lambda_ip <- rep(0, length.out = nrow(network$c_ip))
  
  # Make inputs
  Inputs <- makeInput(family = family2, c_ip = network$c_ip, options = Options_vec, X = X_ij, t_i = network$t_i, version = Version, CalcSD_lambda_ip = Calc_lambda_ip)
  
  try({
    # run model
    mod_out <- runOUGMRF(inputs = Inputs)
    
    #str(mod_out$Report)
    
    theta_hat[i] <- mod_out$Report$theta
    SD_hat[i] <- mod_out$Report$SDinput
  })
  
}
df_reps2 <- data.frame(theta_hat, SD_hat)
df_reps2[which(df_reps2$theta > 100), ] <- NA
format(df_reps2, digits = 2, scientific = FALSE)

summary(df_reps2)

# multiple by 10
family3 <- family
faimly3$dist <- family$dist*10
n_iters <- 100
theta_hat <- rep(NA, length.out = n_iters)
SD_hat <- rep(NA, length.out = n_iters)
for(i in 1:n_iters) {
  network <- simOUGMRF(family = family2, 
                       theta = 0.2, # adult estimate from sp + t model
                       SD = 0.1, # adult estimate from sp + t model
                       mean_N = 40,
                       gamma = gamma_j,
                       X_ij = X_ij,
                       p = p,
                       sample_pct = sample_pct_vec[1],
                       spatial = TRUE
  )
  
  Options_vec = c("SpatialTF"=1, "TemporalTF"=0, "SpatiotemporalTF"=0, "DetectabilityTF"=1, "ObsModel"=1, "OverdispersedTF"=0, "abundTF"=1)
  
  Calc_lambda_ip <- rep(0, length.out = nrow(network$c_ip))
  
  # Make inputs
  Inputs <- makeInput(family = family3, c_ip = network$c_ip, options = Options_vec, X = X_ij, t_i = network$t_i, version = Version, CalcSD_lambda_ip = Calc_lambda_ip)
  
  try({
    # run model
    mod_out <- runOUGMRF(inputs = Inputs)
    
    #str(mod_out$Report)
    
    theta_hat[i] <- mod_out$Report$theta
    SD_hat[i] <- mod_out$Report$SDinput
  })
  
}
df_reps3 <- data.frame(theta_hat, SD_hat)
df_reps3[which(df_reps3$theta > 100), ] <- NA
format(df_reps3, digits = 2, scientific = FALSE)

summary(df_reps3)




# multiple by 10
family4 <- family
family4$dist <- family$dist
n_iters <- 100
theta_hat <- rep(NA, length.out = n_iters)
SD_hat <- rep(NA, length.out = n_iters)
sig_b_rmse <- rep(NA, length.out = n_iters)
N_rmse <- rep(NA, length.out = n_iters)
theta <- 1
SD <- 0.5
for(i in 1:n_iters) {
  network <- simOUGMRF(family = family4, 
                       theta = theta, # adult estimate from sp + t model
                       SD = SD, # adult estimate from sp + t model
                       mean_N = 40,
                       gamma = gamma_j,
                       X_ij = X_ij,
                       p = p,
                       sample_pct = sample_pct_vec[1],
                       spatial = TRUE
  )
  
  Options_vec = c("SpatialTF"=1, "TemporalTF"=0, "SpatiotemporalTF"=0, "DetectabilityTF"=1, "ObsModel"=1, "OverdispersedTF"=0, "abundTF"=1)
  
  Calc_lambda_ip <- rep(0, length.out = nrow(network$c_ip))
  
  # Make inputs
  Inputs <- makeInput(family = family4, c_ip = network$c_ip, options = Options_vec, X = X_ij, t_i = network$t_i, version = Version, CalcSD_lambda_ip = Calc_lambda_ip)
  
  try({
    # run model
    mod_out <- runOUGMRF(inputs = Inputs)
    
    #str(mod_out$Report)
    
    sig_b <- sqrt((SD^2)/(2 * theta) * (1 - exp(-1 * theta * family$dist_b)))
    
    if(mod_out$opt$convergence == 1) {
      theta_hat[i] <- NA
      SD_hat[i] <- NA
      N_rmse[i] <- NA
      sig_b_rmse[i] <- NA
      sigma_b <- data.frame(iter = i, sig_b = NA, sig_b_hat = NA)
    } else {
      theta_hat[i] <- mod_out$Report$theta
      SD_hat[i] <- mod_out$Report$SDinput
      N_rmse[i] <- rmse(network$N_i - mod_out$Report$N_ip[,1])
      sig_b_rmse[i] <- rmse(sig_b - mod_out$Report$SDinput_b)
      sigma_b <- data.frame(iter = i, sig_b = sig_b, sig_b_hat = mod_out$Report$SDinput_b)
      #plot(sig_b, mod_out$Report$SDinput_b, type = "p")
    }
    if(i == 1) {
      df_sig_4 <- sigma_b
    } else {
      df_sig_4 <- dplyr::bind_rows(df_sig, sigma_b)
    }
  })
  
}
df_reps4 <- data.frame(theta_hat, SD_hat, N_rmse, sig_b_rmse)
df_reps4[which(df_reps4$theta > 100), ] <- NA
format(df_reps4, digits = 2, scientific = FALSE)

summary(df_reps4)

ggplot(df_reps4, aes(theta_hat, SD_hat)) + geom_point() + geom_point(aes(x=theta, y=SD), colour = "red", size = 4)

ggplot(df_sig, aes(sig_b, sig_b_hat)) + geom_point(aes(colour = iter)) + scale_color_gradient() + geom_abline(intercept = 0, slope = 1, colour = "red")


# multiple by 10
family5 <- family
family5$dist <- family$dist/100
n_iters <- 100
theta_hat <- rep(NA, length.out = n_iters)
SD_hat <- rep(NA, length.out = n_iters)
sig_b_rmse <- rep(NA, length.out = n_iters)
N_rmse <- rep(NA, length.out = n_iters)
theta <- 10
SD <- 5
for(i in 1:n_iters) {
  network <- simOUGMRF(family = family5, 
                       theta = theta, # adult estimate from sp + t model
                       SD = SD, # adult estimate from sp + t model
                       mean_N = 40,
                       gamma = gamma_j,
                       X_ij = X_ij,
                       p = p,
                       sample_pct = sample_pct_vec[1],
                       spatial = TRUE
  )
  
  Options_vec = c("SpatialTF"=1, "TemporalTF"=0, "SpatiotemporalTF"=0, "DetectabilityTF"=1, "ObsModel"=1, "OverdispersedTF"=0, "abundTF"=1)
  
  Calc_lambda_ip <- rep(0, length.out = nrow(network$c_ip))
  
  # Make inputs
  Inputs <- makeInput(family = family5, c_ip = network$c_ip, options = Options_vec, X = X_ij, t_i = network$t_i, version = Version, CalcSD_lambda_ip = Calc_lambda_ip)
  
  try({
    # run model
    mod_out <- runOUGMRF(inputs = Inputs)
    
    #str(mod_out$Report)
    
    sig_b <- sqrt((SD^2)/(2 * theta) * (1 - exp(-1 * theta * family$dist_b)))
    
    if(mod_out$opt$convergence == 1) {
      theta_hat[i] <- NA
      SD_hat[i] <- NA
      N_rmse[i] <- NA
      sig_b_rmse[i] <- NA
      sigma_b <- data.frame(iter = i, sig_b = NA, sig_b_hat = NA)
    } else {
      theta_hat[i] <- mod_out$Report$theta
      SD_hat[i] <- mod_out$Report$SDinput
      N_rmse[i] <- rmse(network$N_i - mod_out$Report$N_ip[,1])
      sig_b_rmse[i] <- rmse(sig_b - mod_out$Report$SDinput_b)
      sigma_b <- data.frame(iter = i, sig_b = sig_b, sig_b_hat = mod_out$Report$SDinput_b)
      #plot(sig_b, mod_out$Report$SDinput_b, type = "p")
    }
    if(i == 1) {
      df_sig_5 <- sigma_b
    } else {
      df_sig_5 <- dplyr::bind_rows(df_sig_5, sigma_b)
    }
  })
  
}
df_reps5 <- data.frame(theta_hat, SD_hat, N_rmse, sig_b_rmse)
df_reps5[which(df_reps5$theta > 100), ] <- NA
format(df_reps5, digits = 2, scientific = FALSE)

summary(df_reps5)

ggplot(df_reps5, aes(theta_hat, SD_hat)) + geom_point() + geom_point(aes(x=theta, y=SD), colour = "red", size = 5)

ggplot(df_sig, aes(sig_b, sig_b_hat)) + geom_point(aes(colour = iter)) + scale_color_gradient() + geom_abline(intercept = 0, slope = 1, colour = "red")

