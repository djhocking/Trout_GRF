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

#######################
# Simulate GMRF following O-U process
#######################

# vary theta, SD, sigmaIID, log_mean, and sample_pct for spatial/non-spatial
spatial_cor <- c("high", "med", "low")
theta_vec <- c(0.5, 1, 5)
SD_vec <- c(0.1, 0.25, 0.5) # 0.01, 0.1, 0.25 work with theta = c(0.5, 1, 5) and mean_N=c(5,10,50)
sigmaIID_vec <- c(0, 1, 2)
mean_N <- c(5, 10, 100)
spatial_vec <- c(TRUE, FALSE)

sample_pct_vec <- c(1, 0.5, 0.25, 0.1, 0.05) # only do this for a couple of above scenarios

# Detection probability for removal sampling
p <- c(0.75, 0.1875, 0.05)

# Covariates
gamma_j <- c(0.2) # doesn't work if I use more than 1 coef
X_ij <- matrix(rnorm(nrow(family), 0, 1), nrow(family), length(gamma_j))

#############################
# Effect of survey density
#############################

# simulate abundance on network
df_N <- matrix(NA, nrow(family), length(sample_pct_vec))
df_coef <- list()
df_sd <- list()
df_rmse <- NA
N_hat <- list()
set.seed(45641)
network <- simOUGMRF(family = family, 
                     theta = theta_vec[2],
                     SD = SD_vec[1],
                     mean_N = mean_N[2],
                     gamma = gamma_j,
                     X_ij = X_ij,
                     p = p,
                     sample_pct = 1
)

saveRDS(list(network=network, 
             sample_pct_vec=sample_pct_vec, 
             theta = theta_vec[2],
             SD = SD_vec[1],
             mean_N = mean_N[2],
             gamma = gamma_j,
             p = p), file = file.path(dir_out, "survey_density_input.RData"))

# hold out percent of data at random
for(i in 1:length(sample_pct_vec)) {
  c_ip <- network$c_ip
  n <- length(network[["N_i"]])
  exclude_rows <- sample(1:n, size = round((1-sample_pct_vec[i])*n, digits = 0), replace = FALSE)
  c_ip[exclude_rows, ] <- NA
  # initial site-years to get SD for lambda (cycle through just for best model)
  start <- 1
  end <- 2
  Calc_lambda_ip <- rep(NA, length.out = nrow(network$c_ip))
  Calc_lambda_ip[start:end] <- 1
  Calc_lambda_ip[is.na(Calc_lambda_ip)] <- 0
  
  Options_vec = c("SpatialTF"=1, "TemporalTF"=0, "SpatiotemporalTF"=0, "DetectabilityTF"=1, "ObsModel"=1, "OverdispersedTF"=0, "abundTF"=1)
  
  # Make inputs
  Inputs <- makeInput(family = family, c_ip = c_ip, options = Options_vec, X = X_ij, t_i = network$t_i, version = Version, CalcSD_lambda_ip = Calc_lambda_ip)
  
  mod <- runOUGMRF(inputs = Inputs)
  
  SD_means <- data.frame(param = names(mod$SD$value), 
                         est = as.numeric(mod$SD$value), 
                         sd = as.numeric(mod$SD$sd), stringsAsFactors = F)
  
  N_est <- SD_means %>%
    dplyr::filter(grepl("N", param))
  
  N_hat[[i]] <- data.frame(N_hat = mod$Report$N_ip[ , 1], N_SD = N_est$sd, lambda_hat = mod$Report$lambda_ip[ , 1], pass_1 = network$c_ip[ , 1], stringsAsFactors = F)
  N_hat[[i]] <- N_hat[[i]] %>%
    dplyr::mutate(N_hat = ifelse(is.na(pass_1), lambda_hat, N_hat))
  
  df_N[ , i] <- N_hat[[i]]$N_hat
  df_coef[[i]] <- as.numeric(mod$SD$value)
  df_sd[[i]] <- as.numeric(mod$SD$sd)
  df_rmse[i] <- rmse(network$N_i - N_hat[[i]]$N_hat)
  
  df_N_plot <- data.frame(N_hat = N_hat[[i]]$N_hat, 
                          N_i = network$N_i,
                          lambda_hat = N_hat[[i]]$lambda_hat,
                          obsTF = ifelse(is.na(network$c_ip[,1]), FALSE, TRUE))
  g <- ggplot(df_N_plot, aes(N_i, N_hat, colour = obsTF)) + geom_point() + geom_abline(aes(0,1), colour = "blue") + ggtitle(paste0(sample_pct_vec[i]*100, " Percent Surveyed")) + theme_bw() + xlim(0, 25) + ylim(0, 25)
  print(g)
  ggsave(filename = paste0("Output/pct", sample_pct_vec[i], ".pdf"), plot = g)
}

df_err <- data.frame(Percent = sample_pct_vec, RMSE = df_rmse, stringsAsFactors = F)
format(df_err, digits = 3)

coef_means <- data.frame(param = names(mod$SD$value), true = c(gamma_j, NA, NA, log(theta_vec[2]), SD_vec[1], NA, NA, network$N_i), matrix(unlist(df_coef), length(df_sd[[1]]), length(sample_pct_vec)), stringsAsFactors = F)
names(coef_means) <- c("param",
                  "true",
                  paste0("pct", sample_pct_vec[1]),
                  paste0("pct", sample_pct_vec[2]),
                  paste0("pct", sample_pct_vec[3]),
                  paste0("pct", sample_pct_vec[4]),
                  paste0("pct", sample_pct_vec[5])
)
format(coef_means, digits = 2, scientific = FALSE)

sd_pct <- data.frame(param = names(mod$SD$value), matrix(unlist(df_sd), length(df_sd[[1]]), length(df_sd)), stringsAsFactors = FALSE)
names(sd_pct) <- c("param",
                   paste0("pct", sample_pct_vec[1]),
                   paste0("pct", sample_pct_vec[2]),
                   paste0("pct", sample_pct_vec[3]),
                   paste0("pct", sample_pct_vec[4]),
                   paste0("pct", sample_pct_vec[5]))

sd_pct_coef <- sd_pct %>%
  dplyr::filter(!grepl("N_", param))
sd_pct_coef

sd_pct_N <- sd_pct %>%
  dplyr::filter(grepl("N_", param))
colMeans(sd_pct_N[ , !(names(sd_pct) %in% c("param"))], na.rm = T) # precision improves with sample size

saveRDS(list(coef_means = coef_means,
             sd_pct = sd_pct,
             df_err = df_err),
        file = file.path(dir_out, "survey_density_results.RData"))


###############################
# vary theta, SD, mean, spatial
###############################
# vary theta, SD, sigmaIID, log_mean, and sample_pct for spatial/non-spatial
spatial_cor <- c("high", "med", "low")
theta_vec <- c(0.5, 1, 5)
SD_vec <- c(0.1, 0.1, 0.25) # 0.01, 0.1, 0.25 work with theta = c(0.5, 1, 5) and mean_N=c(5,10,50)
sigmaIID_vec <- c(0, 1, 2)
mean_N <- c(5, 10, 100)
spatial_vec <- c(TRUE, FALSE)

sample_pct_vec <- c(1, 0.5, 0.25, 0.1, 0.05)

# Detection probability for removal sampling
p <- c(0.75, 0.1875, 0.05)

# Covariates
gamma_j <- c(0.2) # doesn't work if I use more than 1 coef
X_ij <- matrix(rnorm(nrow(family), 0, 1), nrow(family), length(gamma_j))

options_list <- list(c("SpatialTF"=0, "TemporalTF"=0, "SpatiotemporalTF"=0, "DetectabilityTF"=1, "ObsModel"=1, "OverdispersedTF"=0, "abundTF"=1),
                     c("SpatialTF"=1, "TemporalTF"=0, "SpatiotemporalTF"=0, "DetectabilityTF"=1, "ObsModel"=1, "OverdispersedTF"=0, "abundTF"=1))
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

saveRDS(list(spatial_cor = spatial_cor,
             theta_vec = theta_vec,
             SD_vec = SD_vec,
             sigmaIID_vec = sigmaIID_vec,
             mean_N = mean_N,
             spatial_vec = spatial_vec,
             p = p, 
             gamma_j = gamma_j),
        file = "Output/spatial_sims_input.RData"
        )

set.seed(1987654)
# adjust so don't have a new random sample when comparing spatial and non-spatial models

for(i in 1:length(theta_vec)) { 
  for(j in 1:length(SD_vec)) { # theta must be larger than SD?
    for(k in 1:length(mean_N)) {
      for(l in 1:length(options_list)) {
        m <- 1
          counter <- counter + 1
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
                                  theta = theta_vec[i], 
                                  SD_ou = SD_vec[j], 
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
                                    theta = theta_vec[i], 
                                    SD_ou = SD_vec[j], 
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
                                    theta = theta_vec[i], 
                                    SD_ou = SD_vec[j], 
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
          g <- ggplot(df_N_plot, aes(N_i, N_hat, colour = obsTF)) + geom_point() + geom_abline(aes(0,1), colour = "blue") + theme_bw() + ggtitle(paste0("spatial=", as.integer(spatial_vec[m]), "; sp_mod=", Options_vec[["SpatialTF"]], "; theta=", theta_vec[i], "; SD_ou=", SD_vec[j], "; mean_N=", mean_N[k]))
          ggsave(filename = paste0("Output/Figures/model_", counter[i], ".pdf"), plot = g)
        }
      }
    }
  }


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

