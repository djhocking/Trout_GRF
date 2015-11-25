# clear environment
rm(list = ls())
gc()

#######################
# Load libraries
#######################
library(TMB)
library(dplyr)
library(minqa)
source("Functions/Input_Functions.R")
source("Functions/runOUGMRF.R")

# Diagnose problems of convergence and SD estimation
DiagnosticDir <- "Diagnostics/"
# create code directory if doesn't exist
if (!file.exists(DiagnosticDir)) {
  dir.create(DiagnosticDir)
}

#######################
# Load data
#######################
load("Data/Prepared_Data_White_River.RData")

# remove year from X_ij now so it doesn't mess with testing of temporal and temporal-spatial mdoels
covs <- X_ij
X_ij <- as.matrix(dplyr::select(covs, length_std, temp_summer_std, prcp_winter_std))

# df = dataframe with all data including sites with multiple passes (multiple instances of each child)
# family = dataframe with unique child rows. Other columns are parents of each child, lat, lon, and other data associated with each child node.
# C_ip matrix of counts at site-year i on electrofish survey pass p
# X_ij matrix of covariates (j) for each site-year i as a design matrix NOT including the intercept (column of 1s)
# t_i vector of length i indicating the survey year for each site-year visit
# df_stds dataframe with three columns of the parameter name, means, stds used for z-score standardization of the continuously-distributed independent variables in X_ij

#######################
# Fit in TMB
#######################

# v1a -- Original version
# v1b -- added covariates matrix X_ij
# v1c- adds linear predictors to SD output and multinomial count process (HAS A BUG!)
# v1d- adds makes random variation in detection probability random, and fixed bug in detectprob calculation
# v1e- adds Temporal variation as AR1 process
# v1f- adds Spatiotemporal variation as kronecker product of O-U and AR1 process, plus interface to turn off components
# v1g- add IID lognormal variation (representing micro-variation)
#setwd( TmbFile )

################### Compare models with version g ################
Version = "OU_GMRF_v1h"
# v1a -- Original version
# v1b -- added covariates matrix X_ij
# v1c- adds linear predictors to SD output and multinomial count process (HAS A BUG!)
# v1d- adds makes random variation in detection probability random, and fixed bug in detectprob calculation
# v1e- adds Temporal variation as AR1 process
# v1f- adds Spatiotemporal variation as kronecker product of O-U and AR1 process, plus interface to turn off components
# v1g- add IID lognormal variation (representing micro-variation)
#setwd( TmbFile )

# Compile
if(FALSE) {
  dyn.unload(dynlib(paste0("Code/", Version)))
  file.remove( paste0("Code/", Version,c(".o",".dll")) )
}
compile( paste0("Code/", Version,".cpp") )

#----------------- Observation-Detection Only ------------------
# Turn off random effects in v1f (0 means exclude a component, except for ObsModel)
Options_vec = c("SpatialTF"=0, "TemporalTF"=0, "SpatiotemporalTF"=0, "DetectabilityTF"=0, "ObsModel"=1, "OverdispersedTF"=1)

start <- 1
end <- 10
Calc_lambda_ip <- rep(NA, length.out = nrow(c_ip))
Calc_lambda_ip[start:end] <- 1
Calc_lambda_ip[is.na(Calc_lambda_ip)] <- 0

# Make inputs
Inputs <- makeInput(family = family, c_ip = c_ip, options = Options_vec, X = X_ij, t_i = t_i, version = Version, CalcSD_lambda_ip = Calc_lambda_ip)

inputs = Inputs
dyn.load( dynlib(paste0("Code/", Version )))
obj1 <- MakeADFun(data=inputs$Data, parameters=inputs$Params, random=inputs$Random, map=inputs$Map, hessian=TRUE, inner.control=list(maxit=1000) )

# First run
obj1$fn( obj1$par )

# set initial so return and booleans don't fail
opt1 <- NULL
Report1 <- NULL
SD1 <- NULL

try({
  # Run model
  opt1 = nlminb(start=obj1$env$last.par.best[-c(obj1$env$random)], objective=obj1$fn, gradient=obj1$gr, control=list(eval.max=1e4, iter.max=1e4, trace=1, rel.tol=1e-14) )
  opt1[["Param"]] = names( opt1$par )
  opt1[["final_gradient"]] = obj1$gr( opt1$par )
  opt1[["AIC"]] = 2*opt1$objective + 2*length(opt1$par)
  Report1 = obj1$report()
  Report1[["Optimizer"]] <- "nlminb"
  SD1 = sdreport( obj1, bias.correct=TRUE )
  
  # Use BOBYQA optimization if PORTS fails
  #if(is.null(SD1)) {
  if(opt1$convergence == 1 | max(abs(opt1$final_gradient)) > 0.001 | is.null(SD1)) {
    opt1 <- bobyqa(par = obj1$env$last.par.best[-c(obj1$env$random)], fn = obj1$fn)
    opt1[["convergence"]] <- opt1$ierr
    Report1 = obj1$report()
    Report1[["Optimizer"]] <- "bobyqa"
    opt1[["Param"]] = names( opt1$par )
    opt1[["AIC"]] = 2*opt1$fval + 2*length(opt1$par)
    SD1 <- sdreport(obj1, bias.correct=TRUE )
  }
})

foo <- data.frame(N = Report1$N_ip[,1], lambda = Report1$lambda_ip[,1], N_min = rowSums(c_ip))
foo
  
N_sd <- rep(NA, length.out = nrow(c_ip))
SD <- data.frame(var = names(SD1$value), sd = SD1$sd)
SD <- dplyr::filter(SD, var == "SD_report")
N_sd[which(Calc_lambda_ip == 1)] <- SD$sd

foo <- data.frame(N = Report1$N_ip[,1], lambda = Report1$lambda_ip[,1], N_min = rowSums(c_ip), N_sd = N_sd, c_i1 = c_ip[,1], chat_i1 = Report1$chat_ip[,1])
foo

# get next set of SD
start <- 391
end <- 398
Calc_lambda_ip <- rep(NA, length.out = nrow(c_ip))
Calc_lambda_ip[start:end] <- 1
Calc_lambda_ip[is.na(Calc_lambda_ip)] <- 0

Inputs <- makeInput(family = family, c_ip = c_ip, options = Options_vec, X = X_ij, t_i = t_i, version = Version, CalcSD_lambda_ip = Calc_lambda_ip)
inputs = Inputs
dyn.load( dynlib(paste0("Code/", Version )))

obj1 <- MakeADFun(data=inputs$Data, parameters=inputs$Params, random=inputs$Random, map=inputs$Map, hessian=FALSE, inner.control=list(maxit=1000) )

# First run
obj1$fn( obj1$par )

obj1$env$data$calcSD_Lambda_ip <- Calc_lambda_ip

SD1 <- sdreport(obj1, bias.correct=TRUE ) # seems have to rerun the whole model each time

SD <- data.frame(var = names(SD1$value), sd = SD1$sd)
SD <- dplyr::filter(SD, var == "SD_report")
N_sd[which(Calc_lambda_ip == 1)] <- SD$sd

foo <- data.frame(N = Report1$N_ip[,1], lambda = Report1$lambda_ip[,1], N_min = rowSums(c_ip), N_sd = N_sd, c_i1 = c_ip[,1], chat_i1 = Report1$chat_ip[,1])
foo


  mod1 <- runOUGMRF(inputs = Inputs)
#--------------------------------------------------


#----------------- Temporal Only ------------------
# Turn off random effects in v1f (0 means exclude a component, except for ObsModel)
Options_vec = c("SpatialTF"=0, "TemporalTF"=1, "SpatiotemporalTF"=0, "DetectabilityTF"=1, "ObsModel"=1, "OverdispersedTF"=1)

# Make inputs
Inputs <- makeInput(family = family, c_ip = c_ip, options = Options_vec, X = X_ij, t_i = t_i, version = Version)

mod2 <- runOUGMRF(inputs = Inputs)
#--------------------------------------------------

#----------------- Spatial Only ------------------
# Turn off random effects in v1f (0 means exclude a component, except for ObsModel)
Options_vec = c("SpatialTF"=1, "TemporalTF"=0, "SpatiotemporalTF"=0, "DetectabilityTF"=1, "ObsModel"=1, "OverdispersedTF"=1)

# Make inputs
Inputs <- makeInput(family = family, c_ip = c_ip, options = Options_vec, X = X_ij, t_i = t_i, version = Version, CalcSD_lambda_ip = Calc_lambda_ip)

mod3 <- runOUGMRF(inputs = Inputs)
#--------------------------------------------------


#----------------- Spatiotemporal Only ------------------
# Turn off random effects in v1f (0 means exclude a component, except for ObsModel)
Options_vec = c("SpatialTF"=0, "TemporalTF"=0, "SpatiotemporalTF"=1, "DetectabilityTF"=1, "ObsModel"=1, "OverdispersedTF"=1)

# Make inputs
Inputs <- makeInput(family = family, c_ip = c_ip, options = Options_vec, X = X_ij, t_i = t_i, version = Version)

mod4 <- runOUGMRF(inputs = Inputs)
#--------------------------------------------------

#----------------- Temporal + Spatiotemporal ------------------
# Turn off random effects in v1f (0 means exclude a component, except for ObsModel)
Options_vec = c("SpatialTF"=0, "TemporalTF"=1, "SpatiotemporalTF"=1, "DetectabilityTF"=1, "ObsModel"=1)

# Make inputs
Inputs <- makeInput(family = family, c_ip = c_ip, options = Options_vec, X = X_ij, t_i = t_i, version = Version)

mod5 <- runOUGMRF(inputs = Inputs)
#--------------------------------------------------


#----------------- Spatial Temporal Spatiotemporal ------------------
# Turn off random effects in v1f (0 means exclude a component, except for ObsModel)
Options_vec = c("SpatialTF"=1, "TemporalTF"=1, "SpatiotemporalTF"=1, "DetectabilityTF"=1, "ObsModel"=1)

# Make inputs
Inputs <- makeInput(family = family, c_ip = c_ip, options = Options_vec, X = X_ij, t_i = t_i, version = Version)

mod6 <- runOUGMRF(inputs = Inputs)
#--------------------------------------------------

#----------------- Spatial + Temporal ------------------
# Turn off random effects in v1f (0 means exclude a component, except for ObsModel)
Options_vec = c("SpatialTF"=1, "TemporalTF"=1, "SpatiotemporalTF"=0, "DetectabilityTF"=1, "ObsModel"=1, "OverdispersedTF"=1)

# Make inputs
Inputs <- makeInput(family = family, c_ip = c_ip, options = Options_vec, X = X_ij, t_i = t_i, version = Version)

mod7 <- runOUGMRF(inputs = Inputs)
#--------------------------------------------------


#--------------- AIC -------------
# Models successful: 1 & 3 (standard and spatial)
aic_table <- data.frame(mod = c("Obs", "Spatial"), AIC = c(mod1$opt$AIC, mod3$opt$AIC)) 
aic_table <- dplyr::arrange(aic_table, AIC)
aic_table$delta_AIC <- 0
for(i in 2:nrow(aic_table)) {
  aic_table$delta_AIC[i] <- aic_table$AIC[i] - aic_table$AIC[1]
}
aic_table # spatial much better

str(mod3)
mod3$opt$convergence # 0 equal converged successfully
data.frame(Parameter = names(mod3$SD$value), Estimate = mod3$SD$value, SD = mod3$SD$sd)
mod3$SD$par.fixed

save(mod3_White_River = mod3, mod1_White_River = mod1, aic_table_White_River = aic_table, file = "Output/White_River_Summary.RData")


Model <- c("Obs", "Temporal", "Spatial", "Spatiotemporal", "Temporal + ST", "S+T+ST", "Spatial + Temporal")
M_num <- 1:length(Model)
AIC <- c(opt1$AIC, opt2$AIC, opt3$AIC, opt4$AIC, opt5$AIC, opt6$AIC, opt7$AIC)
aic_table <- data.frame(M_num, Model, AIC, stringsAsFactors = FALSE)
names(aic_table) <- c("M_num", "Model", "AIC")
aic_table <- dplyr::arrange(aic_table, AIC)
aic_table$delta_AIC <- 0
for(i in 2:nrow(aic_table)) {
  aic_table$delta_AIC[i] <- aic_table$AIC[i] - aic_table$AIC[1]
}
aic_table
#------------------------------------------------

# compare coefficient estimates
LCI <- SD4$value - (1.96 * SD4$sd) # lower CI rough estimate for best model
UCI <- SD4$value + (1.96 * SD4$sd)

coef_table <- data.frame(Parameter = names(SD4$value), Estimate = SD4$value, SD = SD4$sd, LCI, UCI, stringsAsFactors = FALSE)
for(i in 1:ncol(as.matrix(X_ij))) {
  coef_table$Parameter[i] <- colnames(as.matrix(X_ij))[i]
}
format(coef_table, digits = 2, scientific = 5)

SD_table <- data.frame(Parameter = names(SD6$value), 
                       SD1 = SD1$sd, 
                       SD2 = SD2$sd, 
                       SD3 = SD3$sd, 
                       SD4 = SD4$sd, 
                       SD5 = SD5$sd, 
                       SD6 = SD6$sd, 
                       SD7 = SD7$sd,
                       stringsAsFactors = FALSE)
for(i in 1:ncol(as.matrix(X_ij))) {
  SD_table$Parameter[i] <- colnames(as.matrix(X_ij))[i]
}
format(SD_table, digits = 2, scientific = 5)

# save
save.image(file = "Output/White_River.RData")

n = nrow(df)
data.frame(v1 = SD$unbiased$value[1:n], v2 = SD$unbiased$value[(n+1):(2*n)], v3 = SD$unbiased$value[(2*n+1):(3*n)])

N_unbias <- SD$unbiased$value[1:n]
N_sd <- SD$sd[1:n]

if(!exists(file.path("Output", Version))) {
  dir.create(path = paste0("Output/", Version, "/"), recursive = TRUE)
} 
capture.output( Sdreport, file=paste0("Output/", Version, "SD.txt"))

# compare predicted vs. observed on original scale
c_est <- Report$lambda_ip * Report$detectprob_ip

plot(c_ip[ , 1], c_est[ , 1], xlab = "Count Observed", ylab = "Count Estimated")
points(c_ip[ , 2], c_est[ , 2], pch = 2)
points(c_ip[ , 3], c_est[ , 3], pch = 3)
abline( a=0, b=1, lty="dotted")


# get outpts for maping rho and N
N_i <- Report$lambda_ip[ , 1]

df_N <- data.frame(child_name = df$child_name, year = df$year, c_sum = rowSums(df[ , c("pass_1", "pass_2", "pass_3")]), N_i, N_unbias, N_sd, p = Report$detectprob_ip, pass_1 = df$pass_1, pass_2 = df$pass_2, pass_3 = df$pass_3)

# only 1 rho per node - add child_name and join back to df_N if rho doesn't change by year
rho_b = data.frame(child_name = family$child_name, rho_b = Report$rho_b)
df_N <- left_join(df_N, rho_b) 

# check if predictions of N are at least as large as the number of individuals caught (assume: complete closure during removal passes and no mis-identification)
df_N <- df_N %>%
  dplyr::mutate(p_miss_total = (1-p.1)^3,
               # problem = ifelse(c_sum > N_i, TRUE, FALSE),
                N_100 = N_i / df$length_sample * 100)
df_N

df_N %>%
  dplyr::select(c_sum, N_i, N_unbias, N_sd, pass_1, pass_2, pass_3, p_miss_total, N_100) %>%
  dplyr::filter(!is.na(c_sum)) 

# site getting the same values at a given site across years!!!!
df_N %>%
  dplyr::select(child_name, year, c_sum, N_i, N_unbias, N_sd, pass_1, pass_2, pass_3, p_miss_total, N_100) %>%
  dplyr::filter(!is.na(c_sum))

# save output for Kyle to map
df_kyle <- df_N %>%
  dplyr::select(child_name, N_100, rho_b) # if rho doesn't change filter to unique child

write.csv(df_kyle, file = "Output/N_correlation.csv")


