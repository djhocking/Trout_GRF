# clear environment
rm(list = ls())
gc()

#######################
# Load libraries
#######################
library(TMB)
library(dplyr)
source("Functions/Input_Functions.R")

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
X_ij <- as.matrix(dplyr::select(covs, length_std))

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
Options_vec = c("SpatialTF"=0, "TemporalTF"=0, "SpatiotemporalTF"=1, "DetectabilityTF"=0, "ObsModel"=1, "OverdispersedTF"=1)

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
SD1 = sdreport( obj1, bias.correct=TRUE )
#--------------------------------------------------


#----------------- Temporal Only ------------------
# Turn off random effects in v1f (0 means exclude a component, except for ObsModel)
Options_vec = c("SpatialTF"=0, "TemporalTF"=1, "SpatiotemporalTF"=0, "DetectabilityTF"=1, "ObsModel"=1)

# Make inputs
Inputs <- makeInput(family = family, c_ip = c_ip, options = Options_vec, X = X_ij, t_i = t_i, version = Version)

# Make object
dyn.load( dynlib(paste0("Code/", Version )))
obj2 <- MakeADFun(data=Inputs$Data, parameters=Inputs$Params, random=Inputs$Random, map=Inputs$Map, hessian=FALSE, inner.control=list(maxit=1000) )
Report = obj2$report()

obj2$gr_orig = obj2$gr
obj2$fn_orig = obj2$fn
obj2$fn <- function( vec ){
  Fn = obj2$gr_orig(vec)
  if( any(is.na(Fn ))) capture.output( matrix(Fn,ncol=1,dimnames=list(names(obj2$par),NULL)), file=paste0(DiagnosticDir,"Fn2.txt"), append = TRUE )
  return( Fn )
}
obj2$gr = function( vec ){
  Gr = obj2$gr_orig(vec)
  if( any(is.na(Gr))) capture.output( matrix(Gr,ncol=1,dimnames=list(names(obj2$par),NULL)), file=paste0(DiagnosticDir,"gr2.txt"), append = TRUE )
  return( Gr )
}

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
#--------------------------------------------------

#----------------- Spatial Only ------------------
# Turn off random effects in v1f (0 means exclude a component, except for ObsModel)
Options_vec = c("SpatialTF"=1, "TemporalTF"=0, "SpatiotemporalTF"=0, "DetectabilityTF"=1, "ObsModel"=1)

# Make inputs
Inputs <- makeInput(family = family, c_ip = c_ip, options = Options_vec, X = X_ij, t_i = t_i, version = Version)

# Make object
dyn.load( dynlib(paste0("Code/", Version )))
obj3 <- MakeADFun(data=Inputs$Data, parameters=Inputs$Params, random=Inputs$Random, map=Inputs$Map, hessian=FALSE, inner.control=list(maxit=1000) )
Report = obj3$report()

obj3$gr_orig = obj3$gr
obj3$fn_orig = obj3$fn
obj3$fn <- function( vec ){
  Fn = obj3$gr_orig(vec)
  if( any(is.na(Fn ))) capture.output( matrix(Fn,ncol=1,dimnames=list(names(obj3$par),NULL)), file=paste0(DiagnosticDir,"Fn3.txt"), append = TRUE )
  return( Fn )
}
obj3$gr = function( vec ){
  Gr = obj3$gr_orig(vec)
  if( any(is.na(Gr))) capture.output( matrix(Gr,ncol=1,dimnames=list(names(obj3$par),NULL)), file=paste0(DiagnosticDir,"gr3.txt"), append = TRUE )
  return( Gr )
}

# First run
obj3$fn( obj3$par )
# Check for parameters that don't do anything
Which = which( obj3$gr( obj3$par )==0 )

# Run model
opt3 = nlminb(start=obj3$env$last.par.best[-c(obj3$env$random)], objective=obj3$fn, gradient=obj3$gr, control=list(eval.max=1e4, iter.max=1e4, trace=1, rel.tol=1e-14) )
opt3[["final_gradient"]] = obj3$gr( opt3$par )
opt3[["AIC"]] = 2*opt3$objective + 2*length(opt3$par)

Report3 = obj3$report()
SD3 = sdreport( obj3, bias.correct=TRUE )
#--------------------------------------------------


#----------------- Spatiotemporal Only ------------------
# Turn off random effects in v1f (0 means exclude a component, except for ObsModel)
Options_vec = c("SpatialTF"=0, "TemporalTF"=0, "SpatiotemporalTF"=1, "DetectabilityTF"=1, "ObsModel"=1)

# Make inputs
Inputs <- makeInput(family = family, c_ip = c_ip, options = Options_vec, X = X_ij, t_i = t_i, version = Version)

# Make object
dyn.load( dynlib(paste0("Code/", Version )))
obj4 <- MakeADFun(data=Inputs$Data, parameters=Inputs$Params, random=Inputs$Random, map=Inputs$Map, hessian=FALSE, inner.control=list(maxit=1000) )
Report = obj4$report()

obj4$gr_orig = obj4$gr
obj4$fn_orig = obj4$fn
obj4$fn <- function( vec ){
  Fn = obj4$gr_orig(vec)
  if( any(is.na(Fn ))) capture.output( matrix(Fn,ncol=1,dimnames=list(names(obj4$par),NULL)), file=paste0(DiagnosticDir,"Fn1.txt"), append = TRUE )
  return( Fn )
}
obj4$gr = function( vec ){
  Gr = obj4$gr_orig(vec)
  if( any(is.na(Gr))) capture.output( matrix(Gr,ncol=1,dimnames=list(names(obj4$par),NULL)), file=paste0(DiagnosticDir,"gr1.txt"), append = TRUE )
  return( Gr )
}
# First run
obj4$fn( obj4$par )
# Check for parameters that don't do anything
Which = which( obj4$gr( obj4$par )==0 )

# Run model
opt4 = nlminb(start=obj4$env$last.par.best[-c(obj4$env$random)], objective=obj4$fn, gradient=obj4$gr, control=list(eval.max=1e4, iter.max=1e4, trace=1, rel.tol=1e-14) )
opt4[["final_gradient"]] = obj4$gr( opt4$par )
opt4[["AIC"]] = 2*opt4$objective + 2*length(opt4$par)

Report4 = obj4$report()
SD4 = sdreport( obj4, bias.correct=TRUE )
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

obj5$gr_orig = obj5$gr
obj5$fn_orig = obj5$fn
obj5$fn <- function( vec ){
  Fn = obj5$gr_orig(vec)
  if( any(is.na(Fn ))) capture.output( matrix(Fn,ncol=1,dimnames=list(names(obj5$par),NULL)), file=paste0(DiagnosticDir,"Fn5.txt"), append = TRUE )
  return( Fn )
}
obj5$gr = function( vec ){
  Gr = obj5$gr_orig(vec)
  if( any(is.na(Gr))) capture.output( matrix(Gr,ncol=1,dimnames=list(names(obj5$par),NULL)), file=paste0(DiagnosticDir,"gr5.txt"), append = TRUE )
  return( Gr )
}

# First run
obj5$fn( obj5$par )
# Check for parameters that don't do anything
Which = which( obj5$gr( obj5$par )==0 )

# Run model
opt5 = nlminb(start=obj5$env$last.par.best[-c(obj5$env$random)], objective=obj5$fn, gradient=obj5$gr, control=list(eval.max=1e4, iter.max=1e4, trace=1, rel.tol=1e-14) )
opt5[["final_gradient"]] = obj5$gr( opt5$par )
opt5[["AIC"]] = 2*opt5$objective + 2*length(opt5$par)

ParHat <- obj5( opt5$par )
Report5 = obj5$report()
SD5 = sdreport( obj5, bias.correct=TRUE )
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

obj6$gr_orig = obj6$gr
obj6$fn_orig = obj6$fn
obj6$fn <- function( vec ){
  Fn = obj6$gr_orig(vec)
  if( any(is.na(Fn ))) capture.output( matrix(Fn,ncol=1,dimnames=list(names(obj6$par),NULL)), file=paste0(DiagnosticDir,"Fn6.txt"), append = TRUE )
  return( Fn )
}
obj6$gr = function( vec ){
  Gr = obj6$gr_orig(vec)
  if( any(is.na(Gr))) capture.output( matrix(Gr,ncol=1,dimnames=list(names(obj6$par),NULL)), file=paste0(DiagnosticDir,"gr6.txt"), append = TRUE )
  return( Gr )
}

# First run
obj6$fn( obj6$par )
# Check for parameters that don't do anything
Which = which( obj6$gr( obj6$par )==0 )

# Run model
opt6 = nlminb(start=obj6$env$last.par.best[-c(obj6$env$random)], objective=obj6$fn, gradient=obj6$gr, control=list(eval.max=1e4, iter.max=1e4, trace=1, rel.tol=1e-14) )
opt6[["final_gradient"]] = obj6$gr( opt6$par )
opt6[["AIC"]] = 2*opt6$objective + 2*length(opt6$par)

Report6 = obj6$report()
SD6 = sdreport( obj6, bias.correct=FALSE )

save(obj6, Report6, SD6, file = "Output/Best_Model_Output.RData")

#--------------------------------------------------

#----------------- Spatial + Temporal ------------------
# Turn off random effects in v1f (0 means exclude a component, except for ObsModel)
Options_vec = c("SpatialTF"=1, "TemporalTF"=1, "SpatiotemporalTF"=0, "DetectabilityTF"=1, "ObsModel"=1)

# Make inputs
Inputs <- makeInput(family = family, c_ip = c_ip, options = Options_vec, X = X_ij, t_i = t_i, version = Version)

# Make object
dyn.load( dynlib(paste0("Code/", Version )))
obj7 <- MakeADFun(data=Inputs$Data, parameters=Inputs$Params, random=Inputs$Random, map=Inputs$Map, hessian=FALSE, inner.control=list(maxit=1000) )
Report = obj7$report()

obj7$gr_orig = obj7$gr
obj7$fn_orig = obj7$fn
obj7$fn <- function( vec ){
  Fn = obj7$gr_orig(vec)
  if( any(is.na(Fn ))) capture.output( matrix(Fn,ncol=1,dimnames=list(names(obj7$par),NULL)), file=paste0(DiagnosticDir,"Fn7.txt") )
  return( Fn )
}
obj7$gr = function( vec ){
  Gr = obj7$gr_orig(vec)
  if( any(is.na(Gr))) capture.output( matrix(Gr,ncol=1,dimnames=list(names(obj7$par),NULL)), file=paste0(DiagnosticDir,"gr7.txt") )
  return( Gr )
}

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
#--------------------------------------------------


#--------------- AIC -------------
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


