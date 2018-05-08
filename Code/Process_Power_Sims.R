#file_list <- list.files("Output/Power_Sim/Data/")

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

load(file = "Output/Power_Sim/Data/ST_Conditions.RData")

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
for(i in 1:n_sim) {
  for(b in 1:length(sample_sites_vec)) {
    for(ti in 1:length(n_years_vec)) {
      for(s in 1:2) {
        counter <- counter + 1
        if(file.exists(paste0("Output/Power_Sim/Data/sim_", i, "_st_", s-1, "_sites_", sample_sites_vec[b], "_years_", n_years_vec[ti], ".RData"))) {
          load(file = paste0("Output/Power_Sim/Data/sim_", i, "_st_", s-1, "_sites_", sample_sites_vec[b], "_years_", n_years_vec[ti], ".RData"))
          if(any(class(mod) == "try-error" | is.null(mod[[1]]) | is.null(mod[[2]]) | is.null(mod[[2]]))) {
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
    }
  }
      }
    }
  }
}

df_sims_extra <- dat
df_sims <- df_sims_extra %>%
  distinct()

df_sims <- df_sims %>%
  dplyr::filter(!is.na(iter),
                !is.na(n_sites),
                !is.na(n_years),
                !is.na(spatialTF),
                !is.na(rhot)) %>%
  dplyr::distinct()

df_sims <- df_sims %>%
  dplyr::mutate(converged = ifelse(is.na(N_se) | N_se > mean_N_est | is.na(mean_N_est) | converge == FALSE, FALSE, TRUE))

df_sims <- df_sims %>%
  dplyr::mutate(converged = ifelse(theta_hat > 10 | theta_st_hat > 10, FALSE, converged))

summary(df_sims)
summary(df_sims[which(df_sims$converged == TRUE), ])

saveRDS(df_sims, file = "Output/Power_Sim/STsim_Results.RData")
write.csv(df_sims, file = "Output/Power_Sim/STsim_Results.csv", row.names = FALSE)

str(df_sims)
















# what predicts prob of convergence
library(lme4)
df_sims_s <- scale(df_sims[ ,2:ncol(df_sims)])
df_sims_s <- data.frame(iter = df_sims$iter, df_sims_s, converge = df_sims$converged)
converge_lm <- glmer(converge ~ n_sites + n_years + spatialTF + n_years * spatialTF + (1|iter), data = df_sims_s, family = "binomial")
summary(converge_lm)

df_sims %>%
  dplyr::group_by(n_sites, spatialTF) %>%
  dplyr::select(-iter) %>%
  dplyr::summarise_each(funs(mean(., na.rm = T)))

df_sims %>%
  dplyr::group_by(n_years, spatialTF) %>%
  dplyr::select(-iter) %>%
  dplyr::summarise_each(funs(mean(., na.rm = T)))







mean_N <- 50
n_years_vec <- c(2, 4, 8, 10, 15, 20)
n_years <- max(n_years_vec)
sample_sites_vec <- c(10, 25, 50, 100, nrow(family))
p <- c(0.75, 0.75, 0.75)
theta <- 0.1
SD <- 0.05
rhot <- 0.8
SD_t <- 0.1
theta_st <- 0.1
SD_st <- 0.05
rho <- 0.5

n_sim <- 100
counter <- 0
for(i in 1:n_sim) {
  for(b in 1:length(sample_sites_vec)) {
    for(t in 1:length(n_years_vec)) {
      for(s in 1:2) {
        counter <- counter + 1
        if(file.exists(paste0("Output/Power_Sim/Data/sim_", i, "_st_", s-1, "_sites_", sample_sites_vec[b], "_years_", n_years_vec[t], ".RData"))) {
          load(file = paste0("Output/Power_Sim/Data/sim_", i, "_st_", s-1, "_sites_", sample_sites_vec[b], "_years_", n_years_vec[t], ".RData"))
          df_N <- data.frame(N_i = network$N_i, N_hat = mod$Report$N_ip[,1])
          if(class(mod) == "try-error") {
            dat[counter, "iter"] <- i
            dat[counter, "n_sites"] <- sample_sites_vec[b]
            dat[counter, "n_years"] <- n_years_vec[t]
            dat[counter, "spatialTF"] <- s - 1
            dat[counter, "mean_N"] <- mean(network$N_i)
            dat[counter, "mean_N_est"] <- NA_real_
            dat[counter, "N_se"] <- NA_real_
            dat[counter, "RMSE"] <- NA_real_
            dat[counter, "theta"] <- theta
            dat[counter, "theta_hat"] <- NA_real_
            dat[counter, "rhot"] <- rhot
            dat[counter, "rhot_hat"] <- NA_real_
            dat[counter, "theta_st"] <- theta_st
            dat[counter, "theta_st_hat"] <- NA_real_
          } else {
            N_se <- mod$SD$sd[which(names(mod$SD$value) == "mean_N")]
            
            dat[counter, "iter"] <- i
            dat[counter, "n_sites"] <- sample_sites_vec[b]
            dat[counter, "n_years"] <- n_years_vec[t]
            dat[counter, "spatialTF"] <- s - 1
            dat[counter, "mean_N"] <- mean(network$N_i)
            dat[counter, "mean_N_est"] <- mod$Report$mean_N
            dat[counter, "N_se"] <- N_se
            dat[counter, "RMSE"] <- rmse(df_N$N_i - df_N$N_hat)
            dat[counter, "theta"] <- theta
            dat[counter, "theta_hat"] <- mod$Report$theta
            dat[counter, "rhot"] <- rhot
            dat[counter, "rhot_hat"] <- mod$Report$rhot
            dat[counter, "theta_st"] <- theta_st
            dat[counter, "theta_st_hat"] <- mod$Report$theta_st
          }
        } else {
          dat[counter, "iter"] <- i
          dat[counter, "n_sites"] <- sample_sites_vec[b]
          dat[counter, "n_years"] <- n_years_vec[t]
          dat[counter, "spatialTF"] <- s - 1
          dat[counter, "mean_N"] <- NA_real_
          dat[counter, "mean_N_est"] <- NA_real_
          dat[counter, "N_se"] <- NA_real_
          dat[counter, "RMSE"] <- NA_real_
          dat[counter, "theta"] <- theta
          dat[counter, "theta_hat"] <- NA_real_
          dat[counter, "rhot"] <- rhot
          dat[counter, "rhot_hat"] <- NA_real_
          dat[counter, "theta_st"] <- theta_st
          dat[counter, "theta_st_hat"] <- NA_real_
        }# end if file exists statement
      }
    }
  }
}

df_sims <- dat %>%
  dplyr::mutate(converged = ifelse(is.na(N_se) | N_se > mean_N_est | is.na(mean_N_est) | converge == FALSE, FALSE, TRUE))

saveRDS(df_sims, file = "Output/Power_Sim/STsim_Results.RData")
write.csv(df_sims, file = "Output/Power_Sim/STsim_Results.RData", row.names = FALSE)

str(df_sims)
df_sim_summary <- df_sims %>%
  dplyr::group_by(n_sites, n_years, spatialTF) %>%
  dplyr::select(-iter) %>%
  dplyr::summarise_each(funs(mean(., na.rm = T)))
as.data.frame(df_sim_summary %>% dplyr::filter(spatialTF == 1)) 

# what predicts prob of convergence
library(lme4)
df_sims_s <- scale(df_sims[ ,2:ncol(df_sims)])
df_sims_s <- data.frame(iter = df_sims$iter, df_sims_s, converge = df_sims$converged)
converge_lm <- glmer(converge.1 ~ n_sites + n_years + spatialTF + n_years * spatialTF + (1|iter), data = df_sims_s, family = "binomial")
summary(converge_lm)

df_sims %>%
  dplyr::group_by(n_sites, spatialTF) %>%
  dplyr::select(-iter) %>%
  dplyr::summarise_each(funs(mean(., na.rm = T)))

df_sims %>%
  dplyr::group_by(n_years, spatialTF) %>%
  dplyr::select(-iter) %>%
  dplyr::summarise_each(funs(mean(., na.rm = T)))

SD_array <- array(NA, c(6, 3, 100))


SD_means <- data.frame(param = names(mod$SD$value), 
                       est = as.numeric(mod$SD$value), 
                       sd = as.numeric(mod$SD$sd), stringsAsFactors = F)

# N_est <- SD_means %>%
#   dplyr::filter(grepl("N", param))

N_est <- mod$Report$N_ip[ , 1]

plot(network$N_i, N_est)
abline(0,1, col = "blue")


N_hat[[i]] <- data.frame(N_hat = mod$Report$N_ip[ , 1], N_SD = N_est$sd, lambda_hat = mod$Report$lambda_ip[ , 1], pass_1 = network$c_ip[ , 1], stringsAsFactors = F)
N_hat[[i]] <- N_hat[[i]] %>%
  dplyr::mutate(N_hat = ifelse(is.na(pass_1), lambda_hat, N_hat))
