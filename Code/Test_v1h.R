# Testing_Version1h

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
SD1 = sdreport( obj1, bias.correct=FALSE )
#--------------------------------------------------

#----------------- Observation-Detection Only no overdisersion ------------------
# Turn off random effects in v1f (0 means exclude a component, except for ObsModel)
Options_vec = c("SpatialTF"=0, "TemporalTF"=0, "SpatiotemporalTF"=1, "DetectabilityTF"=0, "ObsModel"=1, "OverdispersedTF"=0)

# Make inputs
Inputs <- makeInput(family = family, c_ip = c_ip, options = Options_vec, X = X_ij, t_i = t_i, version = Version)

# Make object
dyn.load( dynlib(paste0("Code/", Version )))
obj1b <- MakeADFun(data=Inputs$Data, parameters=Inputs$Params, random=Inputs$Random, map=Inputs$Map, hessian=FALSE, inner.control=list(maxit=1000) )
Report1 = obj1b$report()

# First run
obj1b$fn( obj1b$par )
# Check for parameters that don't do anything
Which = which( obj1b$gr( obj1b$par )==0 )

# Run model
opt1b = nlminb(start=obj1b$env$last.par.best[-c(obj1b$env$random)], objective=obj1b$fn, gradient=obj1b$gr, control=list(eval.max=1e4, iter.max=1e4, trace=1, rel.tol=1e-14) )
opt1b[["final_gradient"]] = obj1b$gr( opt1b$par )
opt1b[["AIC"]] = 2*opt1b$objective + 2*length(opt1b$par)

Report1b = obj1b$report()
SD1b = sdreport( obj1b, bias.correct=FALSE )
#--------------------------------------------------

#----------------- Temporal Only ------------------
# Turn off random effects in v1f (0 means exclude a component, except for ObsModel)
Options_vec = c("SpatialTF"=0, "TemporalTF"=1, "SpatiotemporalTF"=0, "DetectabilityTF"=1, "ObsModel"=1, "OverdispersedTF"=1)

# Make inputs
Inputs <- makeInput(family = family, c_ip = c_ip, options = Options_vec, X = X_ij, t_i = t_i, version = Version)

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
SD2 = sdreport( obj2, bias.correct=TRUE )
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
SD3 = sdreport( obj3, bias.correct=TRUE )
#--------------------------------------------------


#----------------- Spatial Only no overdispersion ------------------
# Turn off random effects in v1f (0 means exclude a component, except for ObsModel)
Options_vec = c("SpatialTF"=1, "TemporalTF"=0, "SpatiotemporalTF"=0, "DetectabilityTF"=1, "ObsModel"=1, "OverdispersedTF"=0)

# Make inputs
Inputs <- makeInput(family = family, c_ip = c_ip, options = Options_vec, X = X_ij, t_i = t_i, version = Version)

# Make object
dyn.load( dynlib(paste0("Code/", Version )))
obj3b <- MakeADFun(data=Inputs$Data, parameters=Inputs$Params, random=Inputs$Random, map=Inputs$Map, hessian=FALSE, inner.control=list(maxit=1000) )
Report = obj3b$report()

# First run
obj3b$fn( obj3b$par )
# Check for parameters that don't do anything
Which = which( obj3b$gr( obj3b$par )==0 )

# Run model
opt3b = nlminb(start=obj3b$env$last.par.best[-c(obj3b$env$random)], objective=obj3b$fn, gradient=obj3b$gr, control=list(eval.max=1e4, iter.max=1e4, trace=1, rel.tol=1e-14) )
opt3b[["final_gradient"]] = obj3b$gr( opt3b$par )
opt3b[["AIC"]] = 2*opt3b$objective + 2*length(opt3b$par)

Report3b = obj3b$report()
SD3b = sdreport( obj3b, bias.correct=TRUE )

#--------------------------------------------------

#----------------- Temporal + Spatial ------------------
# Turn off random effects in v1f (0 means exclude a component, except for ObsModel)
Options_vec = c("SpatialTF"=1, "TemporalTF"=1, "SpatiotemporalTF"=0, "DetectabilityTF"=1, "ObsModel"=1, "OverdispersedTF"=1)

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

# Run model
opt4 = nlminb(start=obj4$env$last.par.best[-c(obj4$env$random)], objective=obj4$fn, gradient=obj4$gr, control=list(eval.max=1e4, iter.max=1e4, trace=1, rel.tol=1e-14) )
opt4[["final_gradient"]] = obj4$gr( opt4$par )
opt4[["AIC"]] = 2*opt4$objective + 2*length(opt4$par)

Report4 = obj4$report()
SD4 = sdreport( obj4, bias.correct=TRUE )
#--------------------------------------------------