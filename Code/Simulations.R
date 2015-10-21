# clear environment
rm(list = ls())
gc()

#######################
# Load libraries
#######################
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

#######################
# Simulate GMRF following O-U process
#######################

# vary theta, SD, sigmaIID, log_mean, and sample_pct for spatial/non-spatial
spatial_cor <- c("high", "med", "low")
theta_vec <- c(0.5, 1, 5)
SD_vec <- c(0.1, 0.5, 1)
sigmaIID_vec <- c(0, 1, 2)
mean_N <- c(5, 10, 50)

sample_pct_vec <- c(1, 0.5, 0.25, 0.1, 0.05) # only do this for a couple of above scenarios

# Detection probability for removal sampling
p <- c(0.75, 0.1875, 0.05)

# Covariates
gamma_j <- c(0.2) # doesn't work if I use more than 1 coef
X_ij <- matrix(rnorm(nrow(family), 0, 1), nrow(family), length(gamma_j))

#############################
# Effect of survey density
#############################
df_N <- matrix(NA, nrow(family), length(sample_pct_vec))
df_coef <- list()
df_sd <- list()
df_rmse <- NA
for(i in 1:length(sample_pct_vec)) {
  set.seed(45641)
  network <- simOUGMRF(family = family, 
                       theta = theta_vec[2],
                       SD = SD_vec[1],
                       mean_N = mean_N[2],
                       gamma = gamma_j,
                       X_ij = X_ij,
                       p = p,
                       sample_pct = sample_pct_vec[i]
  )
  
  Options_vec = c("SpatialTF"=1, "TemporalTF"=0, "SpatiotemporalTF"=0, "DetectabilityTF"=1, "ObsModel"=1, "OverdispersedTF"=0)
  
  # Make inputs
  Inputs <- makeInput(family = family, c_ip = network$c_ip, Options_vec = Options_vec, X = X_ij, t_i = network$t_i, version = Version)
  
  mod <- runOUGMRF(inputs = Inputs)
  
  N_hat <- data.frame(N_hat = mod$Report$N_ip[,1], lambda_hat = mod$Report$lambda_ip[ , 1], pass_1 = network$c_ip[ , 1])
  N_hat <- N_hat %>%
    dplyr::mutate(N_hat = ifelse(is.na(pass_1), lambda_hat, N_hat))
  
  df_N[ , i] <- N_hat$N_hat
  df_coef[[i]] <- as.numeric(mod$SD$value)
  df_sd[[i]] <- as.numeric(mod$SD$sd)
  df_rmse[i] <- rmse(network$N_i - N_hat$N_hat)
  
  df_N_plot <- data.frame(N_hat = N_hat$N_hat, N_i = network$N_i, obsTF = ifelse(is.na(network$c_ip[,1]), FALSE, TRUE))
  g <- ggplot(df_N_plot, aes(N_i, N_hat, colour = obsTF)) + geom_point() + geom_abline(aes(0,1), colour = "blue") + ggtitle(paste0(sample_pct_vec[i]*100, " Percent Surveyed")) + theme_bw()
  ggsave(filename = paste0("Output/Figures/pct", sample_pct_vec[i], ".pdf"), plot = g)
}

coefs <- data.frame(parameter = names(mod$SD$value), true = c(gamma_j, theta_vec[2], 1, 0, NA, SD_vec[1]), matrix(unlist(df_coef), length(df_sd[[1]]), length(sample_pct_vec)))
names(coefs) <- c("parameter",
                  "true",
                  paste0("pct", sample_pct_vec[1]),
                  paste0("pct", sample_pct_vec[2]),
                  paste0("pct", sample_pct_vec[3]),
                  paste0("pct", sample_pct_vec[4]),
                  paste0("pct", sample_pct_vec[5])
)
coefs

data.frame(parameter = names(mod$SD$value), matrix(unlist(df_sd), length(df_sd[[1]]), length(sample_pct_vec)))


###############################
# vary theta, SD, mean, spatial
###############################
# vary theta, SD, sigmaIID, log_mean, and sample_pct for spatial/non-spatial
spatial_cor <- c("high", "med", "low")
theta_vec <- c(0.5, 1, 5)
SD_vec <- c(0.1, 0.25, 0.5) # 0.01, 0.1, 0.25 work with theta = c(0.5, 1, 5) and mean_N=c(5,10,50)
sigmaIID_vec <- c(0, 1, 2)
mean_N <- c(5, 10, 100)
spatial_vec <- c(TRUE, FALSE)

sample_pct_vec <- c(1, 0.5, 0.25, 0.1, 0.05)

# Detection probability for removal sampling
p <- c(0.75, 0.1875, 0.05)

# Covariates
gamma_j <- c(0.2) # doesn't work if I use more than 1 coef
X_ij <- matrix(rnorm(nrow(family), 0, 1), nrow(family), length(gamma_j))

options_list <- list(c("SpatialTF"=0, "TemporalTF"=0, "SpatiotemporalTF"=0, "DetectabilityTF"=1, "ObsModel"=1, "OverdispersedTF"=0),
                     c("SpatialTF"=1, "TemporalTF"=0, "SpatiotemporalTF"=0, "DetectabilityTF"=1, "ObsModel"=1, "OverdispersedTF"=0))
I <- length(theta_vec)
J <- length(SD_vec)
K <- length(mean_N)
L <- length(options_list)
M <- length(spatial_vec)
df_N <- matrix(NA, nrow(family), I*J*K*L*M)
df_coef <- list()
df_sd <- list()
df_rmse <- NA
counter <- 0
mod2 <- list()
set.seed(19876544)
for(i in 1:length(theta_vec)) { 
  for(j in 1:length(SD_vec)) { # theta must be larger than SD?
    for(k in 1:length(mean_N)) {
      for(l in 1:length(options_list)) {
        for(m in 1:length(spatial_vec)) {
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
          Inputs <- makeInput(family = family, c_ip = network$c_ip, Options_vec = Options_vec, X = X_ij, t_i = network$t_i, version = Version)
          
          # run model
          optim_err <- tryCatch(
            mod_out <- runOUGMRF(inputs = Inputs),
            error = function(e) e
          )
          
          if(inherits(optim_err, "error")) {
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
}

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










###############################
# vary theta, SD, mean, spatial
###############################
# vary theta, SD, sigmaIID, log_mean, and sample_pct for spatial/non-spatial
spatial_cor <- c("high", "med", "low")
theta_vec <- c(0.5, 1, 5)
SD_vec <- theta_vec/2
sigmaIID_vec <- c(0, 1, 2)
mean_N <- c(5, 10, 50)

sample_pct_vec <- c(1, 0.5, 0.25, 0.1, 0.05)

# Detection probability for removal sampling
p <- c(0.75, 0.1875, 0.05)

# Covariates
gamma_j <- c(0.2) # doesn't work if I use more than 1 coef
X_ij <- matrix(rnorm(nrow(family), 0, 1), nrow(family), length(gamma_j))

options_list <- list(c("SpatialTF"=0, "TemporalTF"=0, "SpatiotemporalTF"=0, "DetectabilityTF"=1, "ObsModel"=1, "OverdispersedTF"=0),
                     c("SpatialTF"=1, "TemporalTF"=0, "SpatiotemporalTF"=0, "DetectabilityTF"=1, "ObsModel"=1, "OverdispersedTF"=0))
I <- length(theta_vec)
#J <- length(SD_vec)
K <- length(mean_N)
L <- length(options_list)
df_N <- matrix(NA, nrow(family), I*K*L) #*J
df_coef <- list()
df_sd <- list()
df_rmse <- NA
counter <- 0
mod2 <- list()
set.seed(9876544)
for(i in 1:length(theta_vec)) { 
  # for(j in 1:length(SD_vec)) { # theta must be larger than SD?
  for(k in 1:length(mean_N)) {
    for(l in 1:length(options_list)) {
      counter <- counter + 1
      network <- simOUGMRF(family = family, 
                           theta = theta_vec[i],
                           SD = SD_vec[i],
                           mean_N = mean_N[k],
                           gamma = gamma_j,
                           X_ij = X_ij,
                           p = p,
                           sample_pct = sample_pct_vec[2]
      )
      
      Options_vec = options_list[[l]]
      
      # Make inputs
      Inputs <- makeInput(family = family, c_ip = network$c_ip, Options_vec = Options_vec, X = X_ij, t_i = network$t_i, version = Version)
      
      mod2[[counter]] <- runOUGMRF(inputs = Inputs)
      
      N_hat <- data.frame(N_hat = mod2[[counter]]$Report$N_ip[,1], lambda_hat = mod2[[counter]]$Report$lambda_ip[ , 1], pass_1 = network$c_ip[ , 1])
      N_hat <- N_hat %>%
        dplyr::mutate(N_hat = ifelse(is.na(pass_1), lambda_hat, N_hat))
      
      df_N[ , counter] <- N_hat$N_hat
      df_coef[[counter]] <- as.numeric(mod2[[counter]]$SD$value)
      df_sd[[counter]] <- as.numeric(mod2[[counter]]$SD$sd)
      df_rmse[counter] <- rmse(network$N_i - N_hat$N_hat)
      
      #       df_N_plot <- data.frame(N_hat = N_hat$N_hat, N_i = network$N_i, obsTF = ifelse(is.na(network$c_ip[,1]), FALSE, TRUE))
      #       g <- ggplot(df_N_plot, aes(N_i, N_hat, colour = obsTF)) + geom_point() + geom_abline(aes(0,1), colour = "blue") + theme_bw() #+ ggtitle(paste0(sample_pct_vec[i]*100, " Percent Surveyed")) 
      #       ggsave(filename = paste0("Output/Figures/pct", sample_pct_vec[i], ".pdf"), plot = g)
    }
  }
  #  }
}







for(d in 1:length(theta_vec)) { 
  theta <- theta_vec[d]
  for(k in 1:length(SD_vec)) {
    SD <- SD_vec[k]
    
    
    # parameters
    theta = 10
    log_theta <- log(theta)
    SD = 0.125
    log_mean = 3
    sigmaIID <- 0.7
    detectrate <- 1.77
    
    # object
    condSD_b = x_b = rep(NA, nrow(family))
    
    # seed at top of network
    WhichRoot = which( is.na(family[,'parent_b']) )
    condSD_b[WhichRoot] = sqrt( SD^2 / 2*theta )
    set.seed(53476)
    x_b[WhichRoot] = rnorm(1, mean=0, sd=condSD_b[WhichRoot])
    
    # Loop through network
    set.seed(1234)
    while( TRUE ){
      for(i in 1:nrow(family)){
        if( is.na(x_b[i]) ){
          SimulatedNodes = which(!is.na(x_b))
          Match = match( family[i,'parent_b'], SimulatedNodes ) # Which
          if(length(Match)==1){
            condSD_b[i] = sqrt( SD^2/(2*theta) * (1-exp(-2*theta*family[i,'dist_b'])) )
            x_b[i] = x_b[SimulatedNodes[Match]] + rnorm(1, mean=0, sd=condSD_b[i])
          }
        }
      }
      # Stopping condition
      if( all(!is.na(x_b)) ) break()
    }
    
    # Covariates
    gamma_j <- c(0.5, 1, -0.2)
    X_ij <- matrix(rnorm(length(x_b), 0, 1), length(x_b), length(gamma_j))
    eta_i <- gamma_j * X_ij
    
    # Simulate poisson state process
    set.seed(79543)
    N_i = rpois( length(x_b), lambda=exp(x_b+log_mean+eta_i))
    
    # Simulate binomial observation (count) process
    p <- c(0.75, 0.1875, 0.05)
    set.seed(13435)
    c_ip <- matrix(NA, length(N_i), length(p))
    for(a in 1:length(N_i)) {
      for(b in 1:length(p)) {
        c_ip[a,b] <- rbinom(1, N_i[a], p[b])
      }
    }
    
    # temporal - no temporal variability currently
    t_i <- rep(2000, times = length(N_i))
    
    ################
    # Thin observations on the network
    ################
    
    # thin observations at random (to answer how does sampling density within network affect ability to estimate theta)
    sample_pct <- 0.5
    remove_pct <- 1 - sample_pct
    rows <- 1:nrow(c_ip)
    set.seed(18354)
    remove_rows <- sample(rows, size = trunc(length(rows)*remove_pct), replace = FALSE)
    c_ip_reduced <- c_ip
    c_ip_reduced[remove_rows, ] <- NA
    
    #######################
    # Fit in TMB
    #######################
    Version = "OU_GMRF_v1h"
    
    # Compile
    if(FALSE) {
      dyn.unload(dynlib(paste0("Code/", Version)))
      file.remove( paste0("Code/", Version,c(".o",".dll")) )
    }
    compile( paste0("Code/", Version,".cpp") )
    
    #----------------- Observation-Detection Only ------------------
    # Turn off random effects in v1f (0 means exclude a component, except for ObsModel)
    Options_vec = c("SpatialTF"=0, "TemporalTF"=0, "SpatiotemporalTF"=0, "DetectabilityTF"=1, "ObsModel"=1, "OverdispersedTF"=0)
    
    # Make inputs
    Inputs <- makeInput(family = family, c_ip = c_ip_reduced, Options_vec = Options_vec, X = X_ij, t_i = t_i, version = Version)
    
    # Make object
    dyn.load( dynlib(paste0("Code/", Version )))
    obj1 <- MakeADFun(data=Inputs$Data, parameters=Inputs$Params, random=Inputs$Random, map=Inputs$Map, hessian=FALSE, inner.control=list(maxit=1000) )
    
    # First run
    obj1$fn( obj1$par )
    # Check for parameters that don't do anything
    Which = which( obj1$gr( obj1$par )==0 )
    
    # Run model
    opt1 = nlminb(start=obj1$env$last.par.best[-c(obj1$env$random)], objective=obj1$fn, gradient=obj1$gr, control=list(eval.max=1e4, iter.max=1e4, trace=1, rel.tol=1e-14) )
    opt1[["final_gradient"]] = obj1$gr( opt1$par )
    opt1[["AIC"]] = 2*opt1$objective + 2*length(opt1$par)
    Report1 = obj1$report()
    SD1 = sdreport( obj1, bias.correct=FALSE )
    
    # look at theta estimate and SD
    df_coef_1 <- data.frame(Parameter = names(SD1$value), Estimate = as.numeric(SD1$value), SD = SD1$sd)
    
    # compare true and predicted abundance
    N_hat = data.frame(N_hat = Report1$N_ip[,1], lambda_hat = Report1$lambda_ip[ , 1], pass_1 = c_ip_reduced[ , 1])
    N_hat <- N_hat %>%
      dplyr::mutate(N_hat = ifelse(is.na(pass_1), lambda_hat, N_hat))
    df_N1 <- data.frame(N_hat = N_hat$N_hat, N_i, obsTF = ifelse(is.na(c_ip_reduced[,1]), FALSE, TRUE))
    ggplot(df_N1, aes(N_i, N_hat, colour = obsTF)) + geom_point() + geom_abline(aes(0,1), colour = "blue")
    
    rmse(df_N1$N_i - df_N1$N_hat)
    #--------------------------------------------------
    
    
    #----------------- Spatial Only ------------------
    # Turn off random effects in v1f (0 means exclude a component, except for ObsModel)
    Options_vec = c("SpatialTF"=1, "TemporalTF"=0, "SpatiotemporalTF"=0, "DetectabilityTF"=1, "ObsModel"=1, "OverdispersedTF"=0)
    
    # Make inputs
    Inputs <- makeInput(family = family, c_ip = c_ip, Options_vec = Options_vec, X = X_ij, t_i = t_i, version = Version)
    
    # Make object
    dyn.load( dynlib(paste0("Code/", Version )))
    obj3 <- MakeADFun(data=Inputs$Data, parameters=Inputs$Params, random=Inputs$Random, map=Inputs$Map, hessian=FALSE, inner.control=list(maxit=1000) )
    
    # First run
    obj3$fn( obj3$par )
    # Check for parameters that don't do anything
    Which = which( obj3$gr( obj3$par )==0 )
    
    # Run model
    opt3 = nlminb(start=obj3$env$last.par.best[-c(obj3$env$random)], objective=obj3$fn, gradient=obj3$gr, control=list(eval.max=1e4, iter.max=1e4, trace=1, rel.tol=1e-14) )
    opt3[["final_gradient"]] = obj3$gr( opt3$par )
    opt3[["AIC"]] = 2*opt3$objective + 2*length(opt3$par)
    Report3 = obj3$report()
    SD3 = sdreport( obj3, bias.correct=FALSE )
    
    opt3b <- bobyqa(par = obj3$env$last.par.best[-c(obj3$env$random)], fn = obj3$fn)
    Report3b = obj3$report()
    opt3b[["AIC"]] = 2*opt3b$fval + 2*length(opt3b$par)
    SD3b <- sdreport(obj3, bias.correct=FALSE )
    
    # look at theta estimate and SD
    df_coef_3 <- data.frame(Parameter = names(SD3b$value), Estimate = as.numeric(SD3b$value), SD = SD3b$sd)
    
    # compare true and predicted abundance
    N_hat = data.frame(N_hat = Report3b$N_ip[,1], lambda_hat = Report3b$lambda_ip[ , 1], pass_1 = c_ip_reduced[ , 1])
    N_hat <- N_hat %>%
      dplyr::mutate(N_hat = ifelse(is.na(pass_1), lambda_hat, N_hat))
    df_N3 <- data.frame(N_hat = N_hat$N_hat, N_i, obsTF = ifelse(is.na(c_ip_reduced[,1]), FALSE, TRUE))
    ggplot(df_N3, aes(N_i, N_hat, colour = obsTF)) + geom_point() + geom_abline(aes(0,1), colour = "blue")
    
    rmse <- function(error, na.rm = T) {
      sqrt(mean(error^2, na.rm = T))
    }
    rmse(df_N3$N_i - df_N3$N_hat) # fit is better when add spatial even though theta not recovered well
    #--------------------------------------------------
    
    df_coef_1
    df_coef_3
    
    rmse(df_N1$N_i - df_N1$N_hat)
    rmse(df_N3$N_i - df_N3$N_hat)
    
    
    
    
    
    
    
    
    
    