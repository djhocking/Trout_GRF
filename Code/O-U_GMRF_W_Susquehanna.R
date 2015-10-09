# clear environment
rm(list = ls())
gc()
# setwd( "C:/Users/James.Thorson/Desktop/Project_git/Trout_GRF/" )

#######################
# Load libraries
#######################
# install.packages("minqa")
library(TMB)
library(dplyr)
library(minqa)
library(optimx)
source("Functions/Input_Functions.R")

#######################
# Load data
#######################
load("Data/Prepared_Data_W_Susquehanna.RData")

# remove year from X_ij now so it doesn't mess with testing of temporal and temporal-spatial mdoels
#X_ij <- X_ij[ , c("length_std", "effort_std")]
covs <- X_ij
X_ij <- as.matrix(dplyr::select(covs, length_std, effort_std, forest_std, surfcoarse_std))

# df = dataframe with all data including sites with multiple passes (multiple instances of each child)
# family = dataframe with unique child rows. Other columns are parents of each child, lat, lon, and other data associated with each child node.
# C_ip matrix of counts at site-year i on electrofish survey pass p
# X_ij matrix of covariates (j) for each site-year i as a design matrix including the intercept (column of 1s)
# t_i vector of length i indicating the survey year for each site-year visit
# df_stds dataframe with three columns of the parameter name, means, stds used for z-score standardization of the continuously-distributed independent variables in X_ij

#######################
# Fit in TMB
#######################

Version = "OU_GMRF_v1h"
# v1a -- Original version
# v1b -- added covariates matrix X_ij
# v1c- adds linear predictors to SD output and multinomial count process (HAS A BUG!)
# v1d- adds makes random variation in detection probability random, and fixed bug in detectprob calculation
# v1e- adds Temporal variation as AR1 process
# v1f- adds Spatiotemporal variation as kronecker product of O-U and AR1 process, plus interface to turn off components
# v1g- add IID lognormal variation (representing micro-variation)
#setwd( TmbFile )
# v1h - adds option to turn of IID lognormal variation

# Compile
if(FALSE) {
  dyn.unload(dynlib(paste0("Code/", Version)))
  file.remove( paste0("Code/", Version,c(".o",".dll")) )
}
compile( paste0("Code/", Version,".cpp") )

#----------------- Observation-Detection Only ------------------
# Turn off random effects in v1f (0 means exclude a component, except for ObsModel)
Options_vec = c("SpatialTF"=0, "TemporalTF"=0, "SpatiotemporalTF"=0, "DetectabilityTF"=0, "ObsModel"=1, "OverdispersedTF"=1)

# Make inputs
Inputs <- makeInput(family = family, c_ip = c_ip, options = Options_vec, X = X_ij, t_i = t_i, version = Version)

# Make object
dyn.load( dynlib(paste0("Code/", Version )))
obj1 <- MakeADFun(data=Inputs$Data, parameters=Inputs$Params, random=Inputs$Random, map=Inputs$Map, hessian=FALSE, inner.control=list(maxit=1000) )
Report1 = obj1$report()

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
# 
# opt1 <- bobyqa(par = obj1$env$last.par.best[-c(obj1$env$random)], fn = obj1$fn)
# Report1 = obj1$report()
# opt1[["AIC"]] = 2*opt1$fval + 2*length(opt1$par)
# SD1 <- sdreport(obj1, bias.correct=FALSE )

#--------------------------------------------------


#----------------- Temporal Only ------------------
# Turn off random effects in v1f (0 means exclude a component, except for ObsModel)
options_vec = c("SpatialTF"=0, "TemporalTF"=1, "SpatiotemporalTF"=0, "DetectabilityTF"=1, "ObsModel"=1, "OverdispersedTF"=1)

# Make inputs
Inputs <- makeInput(family = family, c_ip = c_ip, options = options_vec, X = X_ij, t_i = t_i, version = Version)

# Make object
dyn.load( dynlib(paste0("Code/", Version )))
obj2 <- MakeADFun(data=Inputs$Data, parameters=Inputs$Params, random=Inputs$Random, map=Inputs$Map, hessian=FALSE, inner.control=list(maxit=1000) )
Report = obj2$report()

# First run
obj2$fn( obj2$par )
# Check for parameters that don't do anything
Which = which( obj2$gr( obj2$par )==0 )

# Run model
opt2 = nlminb(start=obj2$env$last.par.best[-c(obj2$env$random)], objective=obj2$fn, gradient=obj2$gr, control=list(eval.max=1e4, iter.max=1e4, trace=1, rel.tol=1e-14) )
opt2[["final_gradient"]] = obj2$gr( opt2$par )
opt2[["AIC"]] = 2*opt2$objective + 2*length(opt2$par)
Report2 = obj2$report()
SD2 = sdreport( obj2, bias.correct=FALSE )

# 
# opt2 <- bobyqa(par = obj2$env$last.par.best, fn = obj2$fn)
# Report2 = obj2$report( par = obj2$env$last.par.best)
# opt2[["AIC"]] = 2*opt2$fval + 2*length(opt2$par)
# SD2 <- sdreport(obj2, bias.correct=FALSE )
#--------------------------------------------------

#----------------- Spatial Only ------------------
# Turn off random effects in v1f (0 means exclude a component, except for ObsModel)
Options_vec = c("SpatialTF"=1, "TemporalTF"=0, "SpatiotemporalTF"=0, "DetectabilityTF"=1, "ObsModel"=1, "OverdispersedTF"=1)

# Make inputs
Inputs <- makeInput(family = family, c_ip = c_ip, options = Options_vec, X = X_ij, t_i = t_i, version = Version)

# Make object
dyn.load( dynlib(paste0("Code/", Version )))
obj3 <- MakeADFun(data=Inputs$Data, parameters=Inputs$Params, random=Inputs$Random, map=Inputs$Map, hessian=FALSE, inner.control=list(maxit=1000) )
Report = obj3$report()

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

# opt3 <- bobyqa(par = obj3$env$last.par.best[-c(obj3$env$random)], fn = obj3$fn)
# Report3 = obj3$report()
# opt3[["AIC"]] = 2*opt3$fval + 2*length(opt3$par)
# SD3 <- sdreport(obj3, bias.correct=FALSE )
#--------------------------------------------------


#----------------- Spatiotemporal Only ------------------
# Turn off random effects in v1f (0 means exclude a component, except for ObsModel)
Options_vec = c("SpatialTF"=0, "TemporalTF"=0, "SpatiotemporalTF"=1, "DetectabilityTF"=0, "ObsModel"=1)

# Make inputs
Inputs <- makeInput(family = family, c_ip = c_ip, options = Options_vec, X = X_ij, t_i = t_i, version = Version)

# Make object
dyn.load( dynlib(paste0("Code/", Version )))
obj4 <- MakeADFun(data=Inputs$Data, parameters=Inputs$Params, random=Inputs$Random, map=Inputs$Map, hessian=FALSE, inner.control=list(maxit=1000) )
Report = obj4$report()

# First run
obj4$fn( obj4$par )
# Check for parameters that don't do anything
Which = which( obj4$gr( obj4$par )==0 )

# Run model (nlminb is slow)
opt4 = nlminb(start=obj4$env$last.par.best[-c(obj4$env$random)], objective=obj4$fn, gradient=obj4$gr, control=list(eval.max=4e4, iter.max=4e4, trace=4, rel.tol=4e-44) )
opt4[["final_gradient"]] = obj4$gr( opt4$par )
opt4[["AIC"]] = 2*opt4$objective + 2*length(opt4$par)
Report4 = obj4$report()
SD4 = sdreport( obj4, bias.correct=FALSE )

# opt4 <- bobyqa(par = obj4$env$last.par.best[-c(obj4$env$random)], fn = obj4$fn)
# Report4 = obj4$report()
# opt4[["AIC"]] = 2*opt4$fval + 2*length(opt4$par)
# SD4 <- sdreport(obj4, bias.correct=FALSE )

#--------------------------------------------------

#----------------- Temporal + Spatiotemporal ------------------
# Turn off random effects in v1f (0 means exclude a component, except for ObsModel)
Options_vec = c("SpatialTF"=0, "TemporalTF"=1, "SpatiotemporalTF"=1, "DetectabilityTF"=1, "ObsModel"=1)

# Make inputs
Inputs <- makeInput(family = family, c_ip = c_ip, options = Options_vec, X = X_ij, t_i = t_i, version = Version)

# Make object
dyn.load( dynlib(paste0("Code/", Version )))
obj5 <- MakeADFun(data=Inputs$Data, parameters=Inputs$Params, random=Inputs$Random, map=Inputs$Map, hessian=FALSE, inner.control=list(maxit=1000) )
Report = obj5$report()

# First run
obj5$fn( obj5$par )
#fn_test < obj5$fn(obj5$par)

# Check for parameters that don't do anything
Which = which( obj5$gr( obj5$par )==0 )

# Run model
opt5 = nlminb(start=obj5$env$last.par.best[-c(obj5$env$random)], objective=obj5$fn, gradient=obj5$gr, control=list(eval.max=5e4, iter.max=5e4, trace=5, rel.tol=5e-54) )
opt5[["final_gradient"]] = obj5$gr( opt5$par )
opt5[["AIC"]] = 2*opt5$objective + 2*length(opt5$par)
Report5 = obj5$report()
SD5 = sdreport( obj5, bias.correct=FALSE )

# opt5 <- bobyqa(par = obj5$env$last.par.best[-c(obj5$env$random)], fn = obj5$fn)
# Report5 = obj5$report()
# opt5[["AIC"]] = 2*opt5$fval + 2*length(opt5$par)
# SD5 <- sdreport(obj5, bias.correct=FALSE )

#--------------------------------------------------


#----------------- Spatial Temporal Spatiotemporal ------------------
# Turn off random effects in v1f (0 means exclude a component, except for ObsModel)
Options_vec = c("SpatialTF"=1, "TemporalTF"=1, "SpatiotemporalTF"=1, "DetectabilityTF"=1, "ObsModel"=1)

# Make inputs
Inputs <- makeInput(family = family, c_ip = c_ip, options = Options_vec, X = X_ij, t_i = t_i, version = Version)

# Make object
dyn.load( dynlib(paste0("Code/", Version )))
obj6 <- MakeADFun(data=Inputs$Data, parameters=Inputs$Params, random=Inputs$Random, map=Inputs$Map, hessian=FALSE, inner.control=list(maxit=1000) )
Report = obj6$report()

# First run
obj6$fn( obj6$par )
# Check for parameters that don't do anything
Which = which( obj6$gr( obj6$par )==0 )

# Run model
opt6 = nlminb(start=obj6$env$last.par.best[-c(obj6$env$random)], objective=obj6$fn, gradient=obj6$gr, control=list(eval.max=6e4, iter.max=6e4, trace=6, rel.tol=6e-64) )
opt6[["final_gradient"]] = obj6$gr( opt6$par )
opt6[["AIC"]] = 2*opt6$objective + 2*length(opt6$par)
Report6 = obj6$report()
SD6 = sdreport( obj6, bias.correct=FALSE )

# opt6 <- bobyqa(par = obj6$env$last.par.best[-c(obj6$env$random)], fn = obj6$fn)
# Report6 = obj6$report()
# opt6[["AIC"]] = 2*opt6$fval + 2*length(opt6$par)
# SD6 <- sdreport(obj6, bias.correct=FALSE )

#--------------------------------------------------

#----------------- Spatial + Temporal ------------------
# Turn off random effects in v1f (0 means exclude a component, except for ObsModel)
Options_vec = c("SpatialTF"=1, "TemporalTF"=1, "SpatiotemporalTF"=0, "DetectabilityTF"=1, "ObsModel"=1, "OverdispersedTF" = 1)

# Make inputs
Inputs <- makeInput(family = family, c_ip = c_ip, options = Options_vec, X = X_ij, t_i = t_i, version = Version)

# Make object
dyn.load( dynlib(paste0("Code/", Version )))
obj7 <- MakeADFun(data=Inputs$Data, parameters=Inputs$Params, random=Inputs$Random, map=Inputs$Map, hessian=FALSE, inner.control=list(maxit=1000) )
Report = obj7$report()

# First run
obj7$fn( obj7$par )
# Check for parameters that don't do anything
Which = which( obj7$gr( obj7$par )==0 )

# Run model
opt7 = nlminb(start=obj7$env$last.par.best[-c(obj7$env$random)], objective=obj7$fn, gradient=obj7$gr, control=list(eval.max=1e4, iter.max=1e4, trace=1, rel.tol=1e-14) )
opt7[["final_gradient"]] = obj7$gr( opt7$par )
opt7[["AIC"]] = 2*opt7$objective + 2*length(opt7$par)
Report7 = obj7$report()
SD7 = sdreport( obj7, bias.correct=FALSE )

# opt7 <- bobyqa(par = obj7$env$last.par.best[-c(obj7$env$random)], fn = obj7$fn)
# Report7 = obj7$report(obj7$env$last.par.best)
# opt7[["AIC"]] = 2*opt7$fval + 2*length(opt7$par)
# SD7 <- sdreport(obj7, bias.correct=FALSE )

#--------------------------------------------------


#----------------- Spatial + Spatiotemporal ------------------
# Turn off random effects in v1f (0 means exclude a component, except for ObsModel)
Options_vec = c("SpatialTF"=1, "TemporalTF"=0, "SpatiotemporalTF"=1, "DetectabilityTF"=1, "ObsModel"=1)

# Make inputs
Inputs <- makeInput(family = family, c_ip = c_ip, options = Options_vec, X = X_ij, t_i = t_i, version = Version)

# Make object
dyn.load( dynlib(paste0("Code/", Version )))
obj8 <- MakeADFun(data=Inputs$Data, parameters=Inputs$Params, random=Inputs$Random, map=Inputs$Map, hessian=FALSE, inner.control=list(maxit=1000) )
Report = obj8$report()

# First run
obj8$fn( obj8$par )
#fn_test < obj8$fn(obj8$par)

# Check for parameters that don't do anything
Which = which( obj8$gr( obj8$par )==0 )

# Run model
opt8 = nlminb(start=obj8$env$last.par.best[-c(obj8$env$random)], objective=obj8$fn, gradient=obj8$gr, control=list(eval.max=8e4, iter.max=8e4, trace=8, rel.tol=8e-84) )
opt8[["final_gradient"]] = obj8$gr( opt8$par )
opt8[["AIC"]] = 2*opt8$objective + 2*length(opt8$par)
Report8 = obj8$report()
SD8 = sdreport( obj8, bias.correct=FALSE )

# opt8 <- bobyqa(par = obj8$env$last.par.best[-c(obj8$env$random)], fn = obj8$fn)
# Report8 = obj8$report(obj8$env$last.par.best)
# opt8[["AIC"]] = 2*opt8$fval + 2*length(opt8$par)
# SD8 <- sdreport(obj8, bias.correct=FALSE )

#--------------------------------------------------

# Convergence 
c(opt1$message, opt2$message, opt3$message, opt4$message, opt5$message, opt6$message, opt7$message)

c(opt1$ierr, opt2$ierr, opt3$ierr, opt4$ierr, opt5$ierr, opt6$ierr, opt7$ierr, opt8$ierr)

#--------------- AIC -------------
Model <- c("Obs", "Temporal", "Spatial", "Spatiotemporal", "Temporal + ST", "S+T+ST", "Spatial + Temporal", "Spatial + ST") #
M_num <- 1:length(Model)
AIC <- c(opt1$AIC, opt2$AIC, opt3$AIC, opt4$AIC, opt5$AIC, opt6$AIC, opt7$AIC, opt8$AIC) # 
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
LCI <- SD7$value - (1.96 * SD7$sd) # lower CI rough estimate for best model
UCI <- SD7$value + (1.96 * SD7$sd)
  
coef_table <- data.frame(Parameter = names(SD7$value), Estimate = SD7$value, SD = SD7$sd, LCI, UCI, stringsAsFactors = FALSE)
for(i in 1:ncol(as.matrix(X_ij))) {
  coef_table$Parameter[i] <- colnames(as.matrix(X_ij))[i]
}
format(coef_table, digits = 2, scientific = 5)

SD_table <- data.frame(Parameter = names(SD3$value), 
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
save.image(file = "Output/W_Susquehanna.RData")

############### OLD STUFF Below ##############

# Get standard errors
Report6 = obj6$report()
Sdreport = sdreport( obj6 )
SD = sdreport( obj6, bias.correct=TRUE )
SD$unbiased$value

# values and predicitons at sites with data
df_out <- data.frame(df, lambda6 = Report6$lambda_ip[ , 1], lambda5 = Report5$lambda_ip[ , 1]) %>%
  dplyr::filter(!grepl("N_", child_name))
library(tidyr)
df_wide <- df_out %>%
  dplyr::select(child_name, year, lambda5) %>%
  spread(key = year, value = lambda5)
write.csv(df_wide, file = "Output/predictions_by_year.csv")


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


