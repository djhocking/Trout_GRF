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
library(ggplot2)
source("Functions/Input_Functions.R")

#######################
# Load data
#######################
load("Data/Prepared_Data_W_Susquehanna.RData")

# remove year from X_ij now so it doesn't mess with testing of temporal and temporal-spatial mdoels
#X_ij <- X_ij[ , c("length_std", "effort_std")]
covs <- X_ij
X_ij <- as.matrix(dplyr::select(covs, length_std, forest_std, surfcoarse_std))

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
Options_vec = c("SpatialTF"=0, "TemporalTF"=0, "SpatiotemporalTF"=0, "DetectabilityTF"=1, "ObsModel"=1, "OverdispersedTF"=1)

# Make inputs
Inputs <- makeInput(family = family, c_ip = c_ip, Options_vec = Options_vec, X = X_ij, t_i = t_i, version = Version)

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
# 
opt1b <- bobyqa(par = obj1$env$last.par.best[-c(obj1$env$random)], fn = obj1$fn)
Report1b = obj1$report()
opt1b[["AIC"]] = 2*opt1b$fval + 2*length(opt1b$par)
SD1b <- sdreport(obj1, bias.correct=FALSE )

#--------------------------------------------------


#----------------- Temporal Only ------------------
# Turn off random effects in v1f (0 means exclude a component, except for ObsModel)
Options_vec_vec = c("SpatialTF"=0, "TemporalTF"=1, "SpatiotemporalTF"=0, "DetectabilityTF"=1, "ObsModel"=1, "OverdispersedTF"=1)

# Make inputs
Inputs <- makeInput(family = family, c_ip = c_ip, Options_vec = Options_vec_vec, X = X_ij, t_i = t_i, version = Version)

# Make object
dyn.load( dynlib(paste0("Code/", Version )))
obj2 <- MakeADFun(data=Inputs$Data, parameters=Inputs$Params, random=Inputs$Random, map=Inputs$Map, hessian=FALSE, inner.control=list(maxit=1000) )

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
opt2b <- bobyqa(par = obj2$env$last.par.best[-c(obj2$env$random)], fn = obj2$fn)
Report2b = obj2$report()
opt2b[["AIC"]] = 2*opt2b$fval + 2*length(opt2b$par)
SD2b <- sdreport(obj2, bias.correct=FALSE )
#--------------------------------------------------

#----------------- Spatial Only ------------------
# Turn off random effects in v1f (0 means exclude a component, except for ObsModel)
Options_vec = c("SpatialTF"=1, "TemporalTF"=0, "SpatiotemporalTF"=0, "DetectabilityTF"=1, "ObsModel"=1, "OverdispersedTF"=1)

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
#--------------------------------------------------


#----------------- Spatiotemporal Only ------------------
# Turn off random effects in v1f (0 means exclude a component, except for ObsModel)
Options_vec = c("SpatialTF"=0, "TemporalTF"=0, "SpatiotemporalTF"=1, "DetectabilityTF"=1, "ObsModel"=1, "OverdispersedTF"=1)

# Make inputs
Inputs <- makeInput(family = family, c_ip = c_ip, Options_vec = Options_vec, X = X_ij, t_i = t_i, version = Version)

# Make object
dyn.load( dynlib(paste0("Code/", Version )))
obj4 <- MakeADFun(data=Inputs$Data, parameters=Inputs$Params, random=Inputs$Random, map=Inputs$Map, hessian=FALSE, inner.control=list(maxit=1000) )
Report = obj4$report()

# First run
obj4$fn( obj4$par )
# Check for parameters that don't do anything
Which = which( obj4$gr( obj4$par )==0 )

# Run model (nlminb is slow)
opt4 = nlminb(start=obj4$env$last.par.best[-c(obj4$env$random)], objective=obj4$fn, gradient=obj4$gr, control=list(eval.max=1e4, iter.max=1e4, trace=4, rel.tol=1e-14) )
opt4[["final_gradient"]] = obj4$gr( opt4$par )
opt4[["AIC"]] = 2*opt4$objective + 2*length(opt4$par)
Report4 = obj4$report()
SD4 = sdreport( obj4, bias.correct=FALSE )

opt4b <- bobyqa(par = obj4$env$last.par.best[-c(obj4$env$random)], fn = obj4$fn)
Report4b = obj4$report()
opt4b[["AIC"]] = 2*opt4b$fval + 2*length(opt4b$par)
SD4b <- sdreport(obj4, bias.correct=FALSE )

opt4_lbfgsb <- optim(obj4$env$last.par.best[-c(obj4$env$random)], fn = obj4$fn, gr = obj4$gr, method = "L-BFGS-B")
Report4_lbfgsb = obj4$report()
opt4_lbfgsb[["AIC"]] = 2*opt4_lbfgsb$fval + 2*length(opt4_lbfgsb$par)
SD4_lbfgsb<- sdreport(obj4)
#--------------------------------------------------

#----------------- Temporal + Spatiotemporal ------------------
# Turn off random effects in v1f (0 means exclude a component, except for ObsModel)
Options_vec = c("SpatialTF"=0, "TemporalTF"=1, "SpatiotemporalTF"=1, "DetectabilityTF"=1, "ObsModel"=1, "OverdispersedTF"=1)

# Make inputs
Inputs <- makeInput(family = family, c_ip = c_ip, Options_vec = Options_vec, X = X_ij, t_i = t_i, version = Version)

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
opt5 = nlminb(start=obj5$env$last.par.best[-c(obj5$env$random)], objective=obj5$fn, gradient=obj5$gr, control=list(eval.max=1e4, iter.max=1e4, trace=1, rel.tol=1e-14) )
opt5[["final_gradient"]] = obj5$gr( opt5$par )
opt5[["AIC"]] = 2*opt5$objective + 2*length(opt5$par)
Report5 = obj5$report()
SD5 = sdreport( obj5, bias.correct=FALSE )

opt5b <- bobyqa(par = obj5$env$last.par.best[-c(obj5$env$random)], fn = obj5$fn)
Report5b = obj5$report()
opt5b[["AIC"]] = 2*opt5b$fval + 2*length(opt5b$par)
SD5b <- sdreport(obj5, bias.correct=FALSE )

#--------------------------------------------------


#----------------- Spatial Temporal Spatiotemporal ------------------
# Turn off random effects in v1f (0 means exclude a component, except for ObsModel)
Options_vec = c("SpatialTF"=1, "TemporalTF"=1, "SpatiotemporalTF"=1, "DetectabilityTF"=1, "ObsModel"=1, "OverdispersedTF"=1)

# Make inputs
Inputs <- makeInput(family = family, c_ip = c_ip, Options_vec = Options_vec, X = X_ij, t_i = t_i, version = Version)

# Make object
dyn.load( dynlib(paste0("Code/", Version )))
obj6 <- MakeADFun(data=Inputs$Data, parameters=Inputs$Params, random=Inputs$Random, map=Inputs$Map, hessian=FALSE, inner.control=list(maxit=1000) )
Report = obj6$report()

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

opt6b <- bobyqa(par = obj6$env$last.par.best[-c(obj6$env$random)], fn = obj6$fn)
Report6b = obj6b$report()
opt6b[["AIC"]] = 2*opt6b$fval + 2*length(opt6b$par)
SD6b <- sdreport(obj6, bias.correct=FALSE )

#--------------------------------------------------

#----------------- Spatial + Temporal ------------------
# Turn off random effects in v1f (0 means exclude a component, except for ObsModel)
Options_vec = c("SpatialTF"=1, "TemporalTF"=1, "SpatiotemporalTF"=0, "DetectabilityTF"=1, "ObsModel"=1, "OverdispersedTF" = 1)

# Make inputs
Inputs <- makeInput(family = family, c_ip = c_ip, Options_vec = Options_vec, X = X_ij, t_i = t_i, version = Version)

# Make object
dyn.load( dynlib(paste0("Code/", Version )))
obj7 <- MakeADFun(data=Inputs$Data, parameters=Inputs$Params, random=Inputs$Random, map=Inputs$Map, hessian=FALSE, inner.control=list(maxit=1000) )

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

opt7b <- bobyqa(par = obj7$env$last.par.best[-c(obj7$env$random)], fn = obj7$fn)
Report7b = obj7$report(obj7$env$last.par.best)
opt7b[["AIC"]] = 2*opt7b$fval + 2*length(opt7b$par)
SD7b <- sdreport(obj7, bias.correct=FALSE )

#--------------------------------------------------


#----------------- Spatial + Spatiotemporal ------------------
# Turn off random effects in v1f (0 means exclude a component, except for ObsModel)
Options_vec = c("SpatialTF"=1, "TemporalTF"=0, "SpatiotemporalTF"=1, "DetectabilityTF"=1, "ObsModel"=1, "OverdispersedTF"=1)

# Make inputs
Inputs <- makeInput(family = family, c_ip = c_ip, Options_vec = Options_vec, X = X_ij, t_i = t_i, version = Version)

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
opt8 = nlminb(start=obj8$env$last.par.best[-c(obj8$env$random)], objective=obj8$fn, gradient=obj8$gr, control=list(eval.max=1e4, iter.max=1e4, trace=1, rel.tol=1e-14) )
opt8[["final_gradient"]] = obj8$gr( opt8$par )
opt8[["AIC"]] = 2*opt8$objective + 2*length(opt8$par)
Report8 = obj8$report()
SD8 = sdreport( obj8, bias.correct=FALSE )

opt8b <- bobyqa(par = obj8$env$last.par.best[-c(obj8$env$random)], fn = obj8$fn)
Report8b = objb$report(obj8$env$last.par.best)
opt8b[["AIC"]] = 2*opt8b$fval + 2*length(opt8b$par)
SD8b <- sdreport(obj8, bias.correct=FALSE )

#--------------------------------------------------

# Convergence 
df_convergence <- data.frame(model = 1:8, 
                             message = c(opt1$message, opt2$message, opt3$message, opt4$message, opt5$message, opt6$message, opt7$message, opt8$message), 
                             final_gr = c(max(opt1$final_gradient), max(opt2$final_gradient), max(opt3$final_gradient), max(opt4$final_gradient), max(opt5$final_gradient), max(opt6$final_gradient), max(opt7$final_gradient), max(opt8$final_gradient))) %>%
  dplyr::mutate(problem_gr = ifelse(final_gr > 0.001, TRUE, FALSE))
df_convergence

c(opt1b$ierr, opt2b$ierr, opt3b$ierr, opt4b$ierr, opt5b$ierr, opt6b$ierr, opt7b$ierr, opt8b$ierr)


#------------------------------------------------

SD_table <- data.frame(Parameter = names(SD3$value), 
                       SD1 = SD1$sd, 
                       SD2 = SD2$sd, 
                       SD3 = SD3$sd, 
                       #SD4 = SD4r_lbfgsb$sd, 
                       #SD5 = SD5r_lbfgsb$sd, 
                       SD6 = SD6$sd, # fails ports and bobyqa
                       SD7 = SD7$sd,
                       SD8 = SD8$sd, # fails ports and bobyqa
                       stringsAsFactors = FALSE)
for(i in 1:ncol(as.matrix(X_ij))) {
  SD_table$Parameter[i] <- colnames(as.matrix(X_ij))[i]
}
format(SD_table, digits = 2, scientific = 5)

#-------------------- Redo reduced Spatiotemporal -----------------
# In sqrt(diag(object$cov.fixed)) : NaNs produced
# reduce covariates and try
X_ij <- as.matrix(dplyr::select(covs, length_std))

Options_vec = c("SpatialTF"=0, "TemporalTF"=0, "SpatiotemporalTF"=1, "DetectabilityTF"=1, "ObsModel"=1, "OverdispersedTF"=1)

# Make inputs
Inputs <- makeInput(family = family, c_ip = c_ip, Options_vec = Options_vec, X = X_ij, t_i = t_i, version = Version)

# Make object
dyn.load( dynlib(paste0("Code/", Version )))
obj4r <- MakeADFun(data=Inputs$Data, parameters=Inputs$Params, random=Inputs$Random, map=Inputs$Map, hessian=FALSE, inner.control=list(maxit=1000) )

# First run
obj4r$fn( obj4r$par )
# Check for parameters that don't do anything
Which = which( obj4r$gr( obj4r$par )==0 )

opt4r_lbfgsb <- optim(obj4r$env$last.par.best[-c(obj4r$env$random)], fn = obj4r$fn, gr = obj4r$gr, method = "L-BFGS-B")
Report4r_lbfgsb = obj4r$report()
opt4r_lbfgsb[["AIC"]] = 2*opt4r_lbfgsb$fval + 2*length(opt4r_lbfgsb$par)
SD4r_lbfgsb<- sdreport(obj4r)

opt4r_boybqa <- bobyqa(par = obj4r$env$last.par.best[-c(obj4r$env$random)], fn = obj4r$fn)
Report4r_bobyqa = obj4r$report()
opt4r_bobyqa[["AIC"]] = 2*opt4r_bobyqa$fval + 2*length(opt4r_bobyqa$par)
SD4r_bobyqa <- sdreport(obj4r, bias.correct=FALSE )

#-------------------- Redo reduced T + ST -----------------
# In sqrt(diag(object$cov.fixed)) : NaNs produced
# reduce covariates and try
X_ij <- as.matrix(dplyr::select(covs, length_std))

Options_vec = c("SpatialTF"=0, "TemporalTF"=1, "SpatiotemporalTF"=1, "DetectabilityTF"=1, "ObsModel"=1, "OverdispersedTF"=1)

# Make inputs
Inputs <- makeInput(family = family, c_ip = c_ip, Options_vec = Options_vec, X = X_ij, t_i = t_i, version = Version)

# Make object
dyn.load( dynlib(paste0("Code/", Version )))
obj5r <- MakeADFun(data=Inputs$Data, parameters=Inputs$Params, random=Inputs$Random, map=Inputs$Map, hessian=FALSE, inner.control=list(maxit=1000) )

# First run
obj5r$fn( obj5r$par )
# Check for parameters that don't do anything
Which = which( obj5r$gr( obj5r$par )==0 )

opt5r_lbfgsb <- optim(obj5r$env$last.par.best[-c(obj5r$env$random)], fn = obj5r$fn, gr = obj5r$gr, method = "L-BFGS-B", control=list(eval.max=1e4, iter.max=1e4, trace=1, rel.tol=1e-14))
Report5r_lbfgsb = obj5r$report()
opt5r_lbfgsb[["AIC"]] = 2*opt5r_lbfgsb$fval + 2*length(opt5r_lbfgsb$par)
SD5r_lbfgsb<- sdreport(obj5r)

#-------------------- Redo reduced S + T + ST -----------------
# In sqrt(diag(object$cov.fixed)) : NaNs produced
# reduce covariates and try
X_ij <- as.matrix(dplyr::select(covs, length_std))

Options_vec = c("SpatialTF"=1, "TemporalTF"=1, "SpatiotemporalTF"=1, "DetectabilityTF"=1, "ObsModel"=1, "OverdispersedTF"=1)

# Make inputs
Inputs <- makeInput(family = family, c_ip = c_ip, Options_vec = Options_vec, X = X_ij, t_i = t_i, version = Version)

# Make object
dyn.load( dynlib(paste0("Code/", Version )))
obj6r <- MakeADFun(data=Inputs$Data, parameters=Inputs$Params, random=Inputs$Random, map=Inputs$Map, hessian=FALSE, inner.control=list(maxit=1000) )

# First run
obj6r$fn( obj6r$par )
# Check for parameters that don't do anything
Which = which( obj6r$gr( obj6r$par )==0 )

opt6r_lbfgsb <- optim(obj6r$env$last.par.best[-c(obj6r$env$random)], fn = obj6r$fn, gr = obj6r$gr, method = "L-BFGS-B", control=list(eval.max=1e4, iter.max=1e4, trace=1, rel.tol=1e-14))
Report6r_lbfgsb = obj6r$report()
opt6r_lbfgsb[["AIC"]] = 2*opt6r_lbfgsb$fval + 2*length(opt6r_lbfgsb$par)
SD6r_lbfgsb<- sdreport(obj6r)

#----------------------

#-------------------- Redo reduced S + ST -----------------
# In sqrt(diag(object$cov.fixed)) : NaNs produced
# reduce covariates and try
X_ij <- as.matrix(dplyr::select(covs, length_std))

Options_vec = c("SpatialTF"=1, "TemporalTF"=1, "SpatiotemporalTF"=1, "DetectabilityTF"=1, "ObsModel"=1, "OverdispersedTF"=1)

# Make inputs
Inputs <- makeInput(family = family, c_ip = c_ip, Options_vec = Options_vec, X = X_ij, t_i = t_i, version = Version)

# Make object
dyn.load( dynlib(paste0("Code/", Version )))
obj8r <- MakeADFun(data=Inputs$Data, parameters=Inputs$Params, random=Inputs$Random, map=Inputs$Map, hessian=FALSE, inner.control=list(maxit=1000) )

# First run
obj8r$fn( obj8r$par )
# Check for parameters that don't do anything
Which = which( obj8r$gr( obj8r$par )==0 )

opt8r_lbfgsb <- optim(obj8r$env$last.par.best[-c(obj8r$env$random)], fn = obj8r$fn, gr = obj8r$gr, method = "L-BFGS-B", control=list(eval.max=1e4, iter.max=1e4, trace=1, rel.tol=1e-14))
Report8r_lbfgsb = obj8r$report()
opt8r_lbfgsb[["AIC"]] = 2*opt8r_lbfgsb$fval + 2*length(opt8r_lbfgsb$par)
SD8r_lbfgsb<- sdreport(obj8r)

#----------------------

SD4r_lbfgsb$sd # fail
SD5r_lbfgsb$sd # fail
SD6r_lbfgsb$sd # fail
SD8r_lbfgsb$sd # fail

#--------------- AIC -------------
Model <- c("Obs", 
           "Temporal", 
           "Spatial",#, 
          # "Spatiotemporal", 
          # "Temporal + ST", 
          # "S+T+ST", 
           "Spatial + Temporal" 
           #"Spatial + ST"
) #
M_num <- c(1,
           2,
           3, #,
           #4,
           #5,
           #6,
           7 #,
           #8
)
AIC <- c(opt1$AIC, 
         opt2$AIC, 
         opt3$AIC, #, 
         #opt4$AIC, 
         #opt5$AIC, 
         #opt6$AIC, 
         opt7$AIC #, 
         #opt8$AIC
         ) # 
aic_table <- data.frame(M_num, Model, AIC, stringsAsFactors = FALSE)
names(aic_table) <- c("M_num", "Model", "AIC")
aic_table <- dplyr::arrange(aic_table, AIC)
aic_table$delta_AIC <- 0
for(i in 2:nrow(aic_table)) {
  aic_table$delta_AIC[i] <- aic_table$AIC[i] - aic_table$AIC[1]
}
aic_table

# Convergence 
df_convergence <- data.frame(model = 1:8, 
                             message = c(opt1$message, opt2$message, opt3$message, opt4$message, opt5$message, opt6$message, opt7$message, opt8$message), 
                             final_gr = c(max(opt1$final_gradient), max(opt2$final_gradient), max(opt3$final_gradient), max(opt4$final_gradient), max(opt5$final_gradient), max(opt6$final_gradient), max(opt7$final_gradient), max(opt8$final_gradient))) %>%
  dplyr::mutate(problem_gr = ifelse(final_gr > 0.001, TRUE, FALSE))
df_convergence

######### Conclusions #########

# models with spatiotemporal component definitely don't work. data is a bit sparse. few sites have good time series. Good spatial data, okay temporal but poor timeseries at many sites
# Convergence sketchy for most remaining models so redid all with bobyqa to check

c(opt1b$ierr, opt2b$ierr, opt3b$ierr, opt4b$ierr, opt5b$ierr, opt6b$ierr, opt7b$ierr, opt8b$ierr) # all converge (but SD doesn't work for models with spatiotemporal)

Model <- c("Obs", 
           "Temporal", 
           "Spatial",#, 
           # "Spatiotemporal", 
           # "Temporal + ST", 
           # "S+T+ST", 
           "Spatial + Temporal" 
           #"Spatial + ST"
) #
M_num <- c(1,
           2,
           3, #,
           #4,
           #5,
           #6,
           7 #,
           #8
)
AIC <- c(opt1b$AIC, 
         opt2b$AIC, 
         opt3b$AIC, #, 
         #opt4b$AIC, 
         #opt5b$AIC, 
         #opt6b$AIC, 
         opt7b$AIC #, 
         #opt8b$AIC
) # 
aic_table <- data.frame(M_num, Model, AIC, stringsAsFactors = FALSE)
names(aic_table) <- c("M_num", "Model", "AIC")
aic_table <- dplyr::arrange(aic_table, AIC)
aic_table$delta_AIC <- 0
for(i in 2:nrow(aic_table)) {
  aic_table$delta_AIC[i] <- aic_table$AIC[i] - aic_table$AIC[1]
}
aic_table

# compare coefficient estimates
LCI <- SD7b$value - (1.96 * SD7b$sd) # lower CI rough estimate for best model
UCI <- SD7b$value + (1.96 * SD7b$sd)

coef_table <- data.frame(Parameter = names(SD7b$value), Estimate = SD7b$value, SD = SD7b$sd, LCI, UCI, stringsAsFactors = FALSE)
for(i in 1:ncol(as.matrix(X_ij))) {
  coef_table$Parameter[i] <- colnames(as.matrix(X_ij))[i]
}
format(coef_table, digits = 2, scientific = 5)

theta_low <- exp(coef_table[which(coef_table$Parameter == "log_theta"), ]$LCI)
theta_high <- exp(coef_table[which(coef_table$Parameter == "log_theta"), ]$UCI)

dist <- seq(from = 0, to = 10, length.out = 101)
cor_N <- exp(-1*SD7b$value["theta"]*dist)
cor_low <- exp(-1*theta_low*dist)
cor_high <- exp(-1*theta_high*dist)
df_theta <- data.frame(dist, cor_N)
ggplot(df_theta, aes(dist, cor_N)) + geom_line() + theme_bw() + xlab("Distance (km)") + ylab("Mean expected correlation") + geom_line(aes(dist, cor_high), lty = 2) + geom_line(aes(dist, cor_low), lty = 2) # + geom_point()  

# Plot predictions
df_observed <- df %>%
  dplyr::filter(!is.na(pass_1))

lambda_dt <- data.frame(Report7b$lambda_dt)
names(lambda_dt) <- min(t_i):max(t_i)
lambda_dt$child_b <- Inputs$Data$child_b
lambda_dt <- left_join(dplyr::select(df_observed, child_b, child_name, parent_b, NodeLat, NodeLon, featureid), lambda_dt, by = "child_b")
#lambda_dt$child_b <- as.character(lambda_dt$child_b)

foo <- lambda_dt %>%
  dplyr::select(-child_name, -parent_b, -NodeLat, -NodeLon, -featureid) %>%
  tidyr::gather(key = "year", value = lambda, -child_b, convert = T)
foo$lambda <- as.numeric(foo$lambda)

child_list <- unique(df_observed$child_b)
bar <- dplyr::filter(foo, child_b %in% child_list[1:20])

ggplot(bar, aes(year, lambda, group = child_b, colour = child_b)) + geom_line() + geom_point()


# plot observed vs. predicted counts
chat_ip <- Report7b$chat_ip
bar <- data.frame(chat_ip, row = 1:nrow(chat_ip))
bar <- gather(bar, key = row, value = chat, convert = TRUE)
sna <- data.frame(c_ip)
fu <- sna %>% gather(pass, count)
df_counts <- data.frame(fu, bar)
df_counts <- dplyr::filter(df_counts, complete.cases(df_counts))
ggplot(df_counts, aes(count, chat)) + geom_point() + geom_abline(aes(0,1), colour = "blue")

rmse <- function(error, na.rm = T) {
  sqrt(mean(error^2, na.rm = T))
}
rmse(df_counts$)



N_dt <- data.frame(Report3b$N_dt)
names(N_dt) <- min(t_i):max(t_i)
N_dt$child_b <- Inputs$Data$child_b
#lambda_dt$child_b <- as.character(lambda_dt$child_b)

N_ip <- data.frame(N = Report3b$N_ip[ , 1])
N_ip$child_b <- df$child_b
N_ip$year <- df$year
N_ip <- N_ip %>%
  dplyr::arrange(child_b, year)
tail(N_ip, 50)
df_observed <- df %>%
  dplyr::filter(!is.na(pass_1))
N_ip <- left_join(dplyr::select(df_observed, child_name, child_b, parent_b, NodeLat, NodeLon, featureid, pass_1), N_ip, by = "child_b")

child_list <- unique(df_observed$child_b)
bar <- dplyr::filter(N_ip, child_b %in% child_list[1:20])

ggplot(bar, aes(year, N, group = child_b, colour = child_b)) + geom_line() + geom_point() #+ geom_line(aes(year, pass_1))

# Spatial model
# plot predictions
lambda_dt <- data.frame(Report3b$lambda_dt)
names(lambda_dt) <- min(t_i):max(t_i)
lambda_dt$child_b <- Inputs$Data$child_b
lambda_dt <- left_join(lambda_dt, dplyr::select(df, child_name, parent_b, NodeLat, NodeLon, featureid), by = "child_b")
lambda_dt <- left_join(dplyr::select(df_observed, child_name, child_b, parent_b, NodeLat, NodeLon, featureid, pass_1), lambda_dt, by = "child_b")

foo <- lambda_dt %>%
  dplyr::select(-child_name, -parent_b, -NodeLat, -NodeLon, -featureid, -pass_1) %>%
  group_by(child_b) %>%
  tidyr::gather(key = "year", value = lambda, -child_b, convert = T)
foo <- left_join(foo, N_ip[ , c("child_b", "year", "N")]) 

foo <- foo %>%
  dplyr::rename(lambda_hat = lambda) %>%
  dplyr::mutate(lambda_stoc = lambda_hat*exp(rnorm(nrow(foo), 0, Report3b$sigmaIID)),
                N_hat = ifelse(is.na(N), lambda_stoc, N))

# plot time series for subset of sites
child_list <- unique(foo$child_b)
bar <- dplyr::filter(foo, child_b %in% child_list[1:10]) %>%
  dplyr::mutate(child_b = as.character(child_b))

ggplot(bar, aes(year, N_hat, group = child_b, colour = child_b)) + geom_line() + geom_point() + theme_bw() #+ geom_line(aes(year, pass_1))
ggplot(bar, aes(year, N, group = child_b, colour = child_b)) + geom_line() + geom_point() + theme_bw() 


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


