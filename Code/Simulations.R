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

#######################
# Simulate GMRF following O-U process
#######################

# parameters
theta = 0.5
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
gamma_j <- c(0.5)
X_ij <- matrix(rnorm(length(x_b), 0, 1), length(x_b), length(gamma_j))
eta_i <- gamma_j * X_ij

# Simulate poisson state process
set.seed(79543)
N_i = rpois( length(x_b), lambda=exp(x_b+log_mean+eta_i))

# Simulate binomial observation (count) process
p <- c(0.75, 0.1875, 0.05)
set.seed(13435)
c_ip <- matrix(NA, length(N_i), length(p))
for(i in 1:length(N_i)) {
  for(j in 1:length(p)) {
    c_ip[i,j] <- rbinom(1, N_i[i], p[j])
  }
}

# temporal
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

rmse <- function(error, na.rm = T) {
  sqrt(mean(error^2, na.rm = T))
}
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









