# 1. Proof it works
# Spatial model varying theta and sigma. I think it makes sense to compare these to a non-spatial model (sorry a sort of 3rd axis). Do you think we need to vary sigma since it affects the spatial correlation but as a constant and not in relation to distance? I don't think we need to test temporal or spatiotemporal components here.
# 
# 2. Power analysis
# spatiotemporal model varying the number of years and sites with data
# 
# 3. Performance on axis
# 
# I'm unsure if this is necessary for this paper but I could vary the detection rate as you suggest holding everything else constant in a spatial model. This would be relevant for other fish species, other taxa (stream salamanders), and *maybe* YOY vs adults.

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
n_years_vec <- c(4, 8, 10, 15, 20)
n_years <- max(n_years_vec)
sample_sites_vec <- c(25, 50, 100, 200, nrow(family))
p <- c(0.5, 0.5, 0.5)  # Detection probability for each of three passes
theta <- 0.3 # range for spatial variation (works with 1)
SD <- 0.5   #  Marginal SD of spatial variation
rhot <- 0.6 # 0.1 # Autocorrelation over time
SD_t <- 0.2 # 0  # Conditional SD for variation over time
theta_st <- 0.3 # Range for spatio-temporal variation (0.5 works)
SD_st <- 0.4 # 0.15   # Marginal SD of spatial component of spatio-temporal variation
rho <- 0.7    # Correlation among years for spatio-temporal variation

n_sim <- 300

# Covariates
# add spatially varying covariates constant in time
gamma_j <- c(0.2) # doesn't work if I use more than 1 coef
X_i <- matrix(rnorm(nrow(family), 0, 1), nrow(family), length(gamma_j))
# replicate by number of years
X_ij <- do.call("rbind", rep(list(X_i), n_years))

save(mean_N, family, n_years_vec, n_years, sample_sites_vec, p, theta, SD, rhot, SD_t, theta_st, SD_st, rho, n_sim, gamma_j, X_ij, file = "Output/Power_Sim/Data/ST_Conditions.RData")


# Set TMB code
Version = "OU_GMRF_v1i"

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
                  SD = numeric(),
                  SD_hat = numeric(),
                  SD_st = numeric(),
                  SD_st_hat = numeric(),
                  var_sum_sp = numeric(),
                  var_sum_sp_hat = numeric(),
                  SD_inf = numeric(),
                  SD_inf_hat = numeric(),
                  SD_st_inf = numeric(),
                  SD_st_inf_hat = numeric(),
                  rho_st=numeric(),
                  rho_st_hat=numeric(),
                  gamma_j=numeric(),
                  gamma_j_hat=numeric(),
                  converge = logical(),
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
logFile = paste0("Output/Power_Sim/log_file.txt")
# logFile_Finish = paste0(data_dir, "/log_file_finish.txt")
cat("Monitoring progress of prediction loop in parallel", file=logFile, append=FALSE, sep = "\n")
# cat("Monitoring the finish of each loop", file=logFile_Finish, append=FALSE, sep = "\n")

########## Run Parallel Loop ########## 

# problem no file locking when writing in parallel so things seem to get messed up and written in the wrong places
# write.table(dat, file = paste0("Output/Power_Sim/Data/summary.csv"), sep = ",", row.names = FALSE, append = FALSE, col.names = TRUE)

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
  network <- simST(family, theta = theta, SD = SD, rhot = rhot, SD_t = SD_t, theta_st = theta_st, SD_st = SD_st, mean_N = mean_N, n_years = n_years, rho = rho, gamma_j = gamma_j, X_ij=X_ij, p = p, spatial = TRUE, temporal = TRUE, spatiotemporal = TRUE)
  # str(network)
  # summary(network$N_i)
  
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
          Options_vec = c("SpatialTF"=0, "TemporalTF"=1, "SpatiotemporalTF"=0, "DetectabilityTF"=0, "ObsModel"=1, "OverdispersedTF"=0, "abundTF"=0)
        }
        if(s == 2) {
          Options_vec = c("SpatialTF"=1, "TemporalTF"=1, "SpatiotemporalTF"=1, "DetectabilityTF"=0, "ObsModel"=1, "OverdispersedTF"=0, "abundTF"=0)
        }
        
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
          dat[counter, "iter"] <- i
          dat[counter, "n_sites"] <- sample_sites_vec[b]
          dat[counter, "n_years"] <- n_years_vec[ti]
          dat[counter, "spatialTF"] <- s - 1
          dat[counter, "mean_N"] <- mean(network$N_i)
          dat[counter, "min_N"] <- min(network$N_i, na.rm = T)
          dat[counter, "max_N"] <- max(network$N_i, na.rm = T)
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
          dat[counter, "SD"] <- SD
          dat[counter, "SD_hat"] <- NA_real_
          dat[counter, "SD_st"] <- SD_st
          dat[counter, "SD_st_hat"] <- NA_real_
          dat[counter, "var_sum_sp"] <- SD^2 + SD_st^2
          dat[counter, "var_sum_sp_hat"] <- NA_real_
          dat[counter, "SD_inf"] <- SD / ((2 * theta) ^ 0.5)
          dat[counter, "SD_inf_hat"] <- NA_real_
          dat[counter, "SD_st_inf"] <- SD_st / ((2 * theta_st) ^ 0.5)
          dat[counter, "SD_st_inf_hat"] <- NA_real_
          dat[counter, "rho_st"] <- rho
          dat[counter, "rho_st_hat"] <- NA_real_
          dat[counter, "gamma_j"] <- gamma_j
          dat[counter, "gamma_j_hat"] <- NA_real_
          dat[counter, "converge"] <- FALSE
        } else {
          # check convergence
          converge <- FALSE
          try(converge <- mod$opt$convergence == 0)
          try(converge <- ifelse(any(mod$opt[["final_gradient"]] > 0.001), FALSE, converge))
          try(converge <- ifelse(converge == TRUE, !any(is.na(mod$SD$sd)), converge))
          large_sd <- FALSE
          # try(large_sd <- max(mod$SD$sd, na.rm = T) > 100)
          
          if(converge == FALSE | large_sd == TRUE) {
            dat[counter, "iter"] <- i
            dat[counter, "n_sites"] <- sample_sites_vec[b]
            dat[counter, "n_years"] <- n_years_vec[ti]
            dat[counter, "spatialTF"] <- s - 1
            dat[counter, "mean_N"] <- mean(network$N_i)
            dat[counter, "min_N"] <- min(network$N_i, na.rm = T)
            dat[counter, "max_N"] <- max(network$N_i, na.rm = T)
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
            dat[counter, "SD"] <- SD
            dat[counter, "SD_hat"] <- NA_real_
            dat[counter, "SD_st"] <- SD_st
            dat[counter, "SD_st_hat"] <- NA_real_
            dat[counter, "var_sum_sp"] <- SD^2 + SD_st^2
            dat[counter, "var_sum_sp_hat"] <- NA_real_
            dat[counter, "SD_inf"] <- SD / ((2 * theta) ^ 0.5)
            dat[counter, "SD_inf_hat"] <- NA_real_
            dat[counter, "SD_st_inf"] <- SD_st / ((2 * theta_st) ^ 0.5)
            dat[counter, "SD_st_inf_hat"] <- NA_real_
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
          sd_sum <- summary(mod$SD)
          sd_sum <- data.frame(parameter = rownames(sd_sum), sd_sum, stringsAsFactors = FALSE)
          try(df_N <- data.frame(N_i = network$N_i, N_hat = sd_sum[sd_sum$parameter == "N_i", 4]))
          
          dat[counter, "iter"] <- i
          dat[counter, "n_sites"] <- sample_sites_vec[b]
          dat[counter, "n_years"] <- n_years_vec[ti]
          dat[counter, "spatialTF"] <- s - 1
          dat[counter, "mean_N"] <- mean(network$N_i)
          dat[counter, "min_N"] <- min(network$N_i, na.rm = T)
          dat[counter, "max_N"] <- max(network$N_i, na.rm = T)
          dat[counter, "mean_N_est"] <- sd_sum[sd_sum$parameter == "mean_N", 4]
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
          dat[counter, "SD"] <- SD
          dat[counter, "SD_hat"] <- mod$Report$SDinput
          dat[counter, "SD_st"] <- SD_st
          dat[counter, "SD_st_hat"] <- mod$Report$SDinput_st
          dat[counter, "var_sum_sp"] <- mod$Report$SDinput^2 + mod$Report$SDinput_st^2
          dat[counter, "var_sum_sp_hat"] <- SD^2 + SD_st^2
          dat[counter, "SD_inf"] <- SD / ((2 * theta) ^ 0.5)
          dat[counter, "SD_inf_hat"] <- mod$Report$SD_inf
          dat[counter, "SD_st_inf"] <- SD_st / ((2 * theta_st) ^ 0.5)
          dat[counter, "SD_st_inf_hat"] <- mod$Report$SD_st_inf
          dat[counter, "rho_st"] <- rho
          dat[counter, "rho_st_hat"] <- mod$Report$rho_st
          dat[counter, "gamma_j"] <- gamma_j
          dat[counter, "gamma_j_hat"] <- mod$Report$gamma_j
          dat[counter, "converge"] <- TRUE
          }
        }
        dat$rmse <- as.numeric(dat$RMSE)
        #         mean(network$N_i)
        #         SD_means <- data.frame(param = names(mod$SD$value), 
        #                                est = as.numeric(mod$SD$value), 
        #                                sd = as.numeric(mod$SD$sd), stringsAsFactors = F)
        
        #-------- save sim iter output -------
       # save(network, mod, file = paste0("Output/Power_Sim/Data/sim_", i, "_st_", s-1, "_sites_", sample_sites_vec[b], "_years_", n_years_vec[ti], ".RData"))
        
      } # end spatial TF loop
    } # end site loop
  } # end year loop
  
  cat(paste0("Finished sim ", i, " of ", n_sim, ": ", timestamp()), file=logFile, append=TRUE, sep = "\n")
  
  # problem no file locking when writing in parallel so things seem to get messed up and written in the wrong places
  write.table(dat, file = paste0("Output/Power_Sim/Data/summary_", i, ".csv"), sep = ",", row.names = FALSE)
   # save(dat, file = paste0("Output/Power_Sim/Data/summary_sim", i, ".RData"))
  return(dat) # if save before the return nothing gets returned
} # end sim iter
stopCluster(cl)
closeAllConnections()

save(df_sims, file = "Output/Power_Sim/STsim_Results.RData")
write.csv(df_sims, file = "Output/Power_Sim/STsim_Results.csv", row.names = FALSE)


library(readr)
# df_sims <- read_csv(paste0("Output/Power_Sim/Data/summary.csv")) #, header = TRUE, stringsAsFactors = FALSE)
df_sims <- read_csv(paste0("Output/Power_Sim/Data/summary_", 1, ".csv"))
for(sim in 2:n_sim) {
  if(file.exists(paste0("Output/Power_Sim/Data/summary_", sim, ".csv"))) {
    try(foo <- read_csv(paste0("Output/Power_Sim/Data/summary_", sim, ".csv")))
    if(class(foo)[[1]] != "try-error") {
      df_sims <- dplyr::bind_rows(df_sims, foo)
    }
  }
}
write.table(dat, file = paste0("Output/Power_Sim/Data/summary.csv"), sep = ",", row.names = FALSE)

dplyr::filter(df_sims, spatialTF == TRUE) %>% summary()
dplyr::filter(df_sims, spatialTF == TRUE & converge == TRUE) %>% summary()
dplyr::filter(df_sims, spatialTF == FALSE) %>% summary()

save(df_sims, file = "Output/Power_Sim/STsim_Results.RData")
write.csv(df_sims, file = "Output/Power_Sim/STsim_Results.csv", row.names = FALSE)
