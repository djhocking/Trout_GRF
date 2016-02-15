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
n_years_vec <- c(10, 10)
n_years <- max(n_years_vec)
sample_sites_vec <- c(200)
p <- c(0.75, 0.75, 0.75)
theta <- 0.1
SD <- 0.05
rhot <- 0.8
SD_t <- 0.1
theta_st <- 0.1
SD_st <- 0.05
rho <- 0.5
s <- 2 # 1 = nonspatial, 2 = spatial

# Covariates
# add spatially varying covariates constant in time
gamma_j <- c(0.5) # doesn't work if I use more than 1 coef
X_i <- matrix(rnorm(nrow(family), 0, 1), nrow(family), length(gamma_j))
# replicate by number of years
X_ij <- do.call("rbind", rep(list(X_i), n_years))

# Set TMB code
Version = "OU_GMRF_v1h"

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
                  mean_N_est=numeric(),
                  N_se=numeric(),
                  RMSE=numeric(),
                  theta=numeric(),
                  theta_hat=numeric(),
                  rhot=numeric(),
                  rhot_hat=numeric(),
                  theta_st=numeric(),
                  theta_st_hat=numeric(),
                  stringsAsFactors=FALSE)

counter <- 0
# set conditions
df_N <- matrix(NA, nrow(family), length(sample_sites_vec))
df_coef <- list()
df_sd <- list()
df_rmse <- NA
N_hat <- list()

# simulate abundance and counts on network
set.seed(723750)
network <- simST(family, theta = theta, SD = SD, rhot = rhot, SD_t = SD_t, theta_st = theta_st, SD_st = SD_st, mean_N = mean_N, n_years = n_years, rho = rho, gamma_j = gamma_j, X_ij=X_ij, p = p) #, sample_n = sample_sites_vec[s], sample_years = NULL) # 14 minutes for 1 sim of 11,120 node network over 10 years, 2.6 sec for network of 359 nodes over 10 years
str(network)
summary(network$N_i)

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
  dat[counter, "iter"] <- 1
  dat[counter, "n_sites"] <- sample_sites_vec[1]
  dat[counter, "n_years"] <- n_years_vec[1]
  dat[counter, "spatialTF"] <- s - 1
  dat[counter, "mean_N"] <- mean(network$N_i)
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
} else {
  try(N_se <- mod$SD$sd[which(names(mod$SD$value) == "mean_N")])
  N_se <- ifelse(is.null(N_se), NA_real_, N_se)
  if(!is.na(N_se)) {
    N_se <- ifelse(N_se == "NaN", NA_real_, N_se)
  }
  df_N <- data.frame(N_i = network$N_i, N_hat = NA_real_)
  try(df_N <- data.frame(N_i = network$N_i, N_hat = mod$Report$N_ip[,1]))
  
  dat[counter, "iter"] <- 1
  dat[counter, "n_sites"] <- sample_sites_vec[1]
  dat[counter, "n_years"] <- n_years_vec[1]
  dat[counter, "spatialTF"] <- s - 1
  dat[counter, "mean_N"] <- mean(network$N_i)
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
}

N_est <- mod$Report$N_ip[ , 1]
plot(network$N_i, N_est)
abline(0,1, col = "blue")

c_est <- mod$Report$chat_ip[ , 1]
plot(c_ip[ , 1], c_est)
abline(0,1, col = "blue")



