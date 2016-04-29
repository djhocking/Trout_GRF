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
source("Functions/runOUGMRF.R")

#######################
# Load data
#######################
load("Data/Prepared_Data_W_Susquehanna_Combined.RData")

# df_yoy$date <- as.Date(df_yoy$date)
# df_yoy$month <- month(df_yoy$date)
# c_ip_yoy[which(df_yoy$month < 6), c("pass_1", "pass_2", "pass_3")] <- NA

# remove year from X_ij now so it doesn't mess with testing of temporal and temporal-spatial mdoels
#X_ij <- X_ij[ , c("length_std", "effort_std")]
covs <- X_ij
X_ij <- as.matrix(dplyr::select(covs, forest, surfcoarse, temp_mean_fall_1, temp_mean_winter, temp_mean_spring, prcp_mean_fall_1, prcp_mean_winter, prcp_mean_spring))

offset <-  as.numeric(df_yoy$length_sample)

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

# initial site-years to get SD for lambda (cycle through just for best model)
start <- 1
end <- 2
Calc_lambda_ip <- rep(NA, length.out = nrow(c_ip))
Calc_lambda_ip[start:end] <- 1
Calc_lambda_ip[is.na(Calc_lambda_ip)] <- 0

#----------------- Observation-Detection Only ------------------
# Turn off random effects in v1f (0 means exclude a component, except for ObsModel)
Options_vec = c("SpatialTF"=0, "TemporalTF"=0, "SpatiotemporalTF"=0, "DetectabilityTF"=1, "ObsModel"=1, "OverdispersedTF"=1, "abundTF"=0)

# Make inputs
Inputs <- makeInput(family = family, df = df_yoy, c_ip = c_ip_yoy, options = Options_vec, X = X_ij, t_i = t_i, version = Version, offset_i = offset, CalcSD_lambda_ip = Calc_lambda_ip)

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
Options_vec = c("SpatialTF"=0, "TemporalTF"=1, "SpatiotemporalTF"=0, "DetectabilityTF"=1, "ObsModel"=1, "OverdispersedTF"=1, "abundTF"=0)

# Make inputs
Inputs <- makeInput(family = family, df = df_yoy, c_ip = c_ip_yoy, options = Options_vec, X = X_ij, t_i = t_i, version = Version, offset_i = offset, CalcSD_lambda_ip = Calc_lambda_ip)

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
Options_vec = c("SpatialTF"=1, "TemporalTF"=0, "SpatiotemporalTF"=0, "DetectabilityTF"=1, "ObsModel"=1, "OverdispersedTF"=1, "abundTF"=0)

# Make inputs
Inputs <- makeInput(family = family, df = df_yoy, c_ip = c_ip_yoy, options = Options_vec, X = X_ij, t_i = t_i, version = Version, offset_i = offset, CalcSD_lambda_ip = Calc_lambda_ip)

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
Options_vec = c("SpatialTF"=0, "TemporalTF"=0, "SpatiotemporalTF"=1, "DetectabilityTF"=1, "ObsModel"=1, "OverdispersedTF"=1, "abundTF"=0)

# Make inputs
Inputs <- makeInput(family = family, df = df_yoy, c_ip = c_ip_yoy, options = Options_vec, X = X_ij, t_i = t_i, version = Version, offset_i = offset, CalcSD_lambda_ip = Calc_lambda_ip)

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
Options_vec = c("SpatialTF"=0, "TemporalTF"=1, "SpatiotemporalTF"=1, "DetectabilityTF"=1, "ObsModel"=1, "OverdispersedTF"=1, "abundTF"=0)

# Make inputs
Inputs <- makeInput(family = family, df = df_yoy, c_ip = c_ip_yoy, options = Options_vec, X = X_ij, t_i = t_i, version = Version, offset_i = offset, CalcSD_lambda_ip = Calc_lambda_ip)

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
Options_vec = c("SpatialTF"=1, "TemporalTF"=1, "SpatiotemporalTF"=1, "DetectabilityTF"=1, "ObsModel"=1, "OverdispersedTF"=1, "abundTF"=0)

# Make inputs
Inputs <- makeInput(family = family, df = df_yoy, c_ip = c_ip_yoy, options = Options_vec, X = X_ij, t_i = t_i, version = Version, offset_i = offset, CalcSD_lambda_ip = Calc_lambda_ip)

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
Inputs <- makeInput(family = family, df = df_yoy, c_ip = c_ip_yoy, options = Options_vec, X = X_ij, t_i = t_i, version = Version, offset_i = offset, CalcSD_lambda_ip = Calc_lambda_ip)

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
Options_vec = c("SpatialTF"=1, "TemporalTF"=0, "SpatiotemporalTF"=1, "DetectabilityTF"=1, "ObsModel"=1, "OverdispersedTF"=1, "abundTF"=0)

# Make inputs
Inputs <- makeInput(family = family, df = df_yoy, c_ip = c_ip_yoy, options = Options_vec, X = X_ij, t_i = t_i, version = Version, offset_i = offset, CalcSD_lambda_ip = Calc_lambda_ip)

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
                       SD3 = SD3b$sd,
                       SD4 = SD4b$sd, 
                       SD5 = SD5b$sd, 
                       SD6 = SD6$sd, # fails
                       SD7 = SD7$sd,
                       SD8 = SD8$sd, #
                       stringsAsFactors = FALSE)
for(i in 1:ncol(as.matrix(X_ij))) {
  SD_table$Parameter[i] <- colnames(as.matrix(X_ij))[i]
}
format(SD_table, digits = 2, scientific = 5)
 
dplyr::arrange(data.frame(opt1b$AIC, opt2b$AIC, opt3b$AIC, opt4b$AIC, opt5b$AIC, opt6b$AIC, opt7b$AIC, opt8b$AIC))


############## most models fail - try reduced models ############


#----------------- Spatial Only ------------------
X_ij <- as.matrix(dplyr::select(covs, forest, temp_mean_winter, temp_mean_spring, prcp_mean_winter, prcp_mean_spring))

# Turn off random effects in v1f (0 means exclude a component, except for ObsModel)
Options_vec = c("SpatialTF"=1, "TemporalTF"=0, "SpatiotemporalTF"=0, "DetectabilityTF"=1, "ObsModel"=1, "OverdispersedTF"=1, "abundTF"=0)

# Make inputs
Inputs <- makeInput(family = family, df = df_yoy, c_ip = c_ip_yoy, options = Options_vec, X = X_ij, t_i = t_i, version = Version, offset_i = offset, CalcSD_lambda_ip = Calc_lambda_ip)

Mod3r <- runOUGMRF(Inputs)

Mod3r$SD$sd # fails

# Further Reduce
X_ij <- as.matrix(dplyr::select(covs, forest, temp_mean_spring, prcp_mean_winter))

Options_vec = c("SpatialTF"=1, "TemporalTF"=0, "SpatiotemporalTF"=0, "DetectabilityTF"=1, "ObsModel"=1, "OverdispersedTF"=1, "abundTF"=0)

# Make inputs
Inputs <- makeInput(family = family, df = df_yoy, c_ip = c_ip_yoy, options = Options_vec, X = X_ij, t_i = t_i, version = Version, offset_i = offset, CalcSD_lambda_ip = Calc_lambda_ip)

Mod3r2 <- runOUGMRF(Inputs)
Mod3r2$SD$sd # success


#----------------- Spatiotemporal Only ------------------
X_ij <- as.matrix(dplyr::select(covs, forest, temp_mean_winter, temp_mean_spring, prcp_mean_winter, prcp_mean_spring))

# Turn off random effects in v1f (0 means exclude a component, except for ObsModel)
Options_vec = c("SpatialTF"=0, "TemporalTF"=0, "SpatiotemporalTF"=1, "DetectabilityTF"=1, "ObsModel"=1, "OverdispersedTF"=1, "abundTF"=0)

# Make inputs
Inputs <- makeInput(family = family, df = df_yoy, c_ip = c_ip_yoy, options = Options_vec, X = X_ij, t_i = t_i, version = Version, offset_i = offset, CalcSD_lambda_ip = Calc_lambda_ip)

Mod4r <- runOUGMRF(Inputs)

str(Mod4r)
Mod4r$opt$convergence
Mod4r$SD$sd

# Further Reduce
X_ij <- as.matrix(dplyr::select(covs, forest, temp_mean_spring, prcp_mean_winter))

# Turn off random effects in v1f (0 means exclude a component, except for ObsModel)
Options_vec = c("SpatialTF"=0, "TemporalTF"=0, "SpatiotemporalTF"=1, "DetectabilityTF"=1, "ObsModel"=1, "OverdispersedTF"=1, "abundTF"=0)

# Make inputs
Inputs <- makeInput(family = family, df = df_yoy, c_ip = c_ip_yoy, options = Options_vec, X = X_ij, t_i = t_i, version = Version, offset_i = offset, CalcSD_lambda_ip = Calc_lambda_ip)

Mod4r2 <- runOUGMRF(Inputs)
str(Mod4r2)
Mod4r2$opt$convergence
Mod4r2$SD$sd

# final reduction
X_ij <- as.matrix(dplyr::select(covs, forest))

# Turn off random effects in v1f (0 means exclude a component, except for ObsModel)
Options_vec = c("SpatialTF"=0, "TemporalTF"=0, "SpatiotemporalTF"=1, "DetectabilityTF"=1, "ObsModel"=1, "OverdispersedTF"=1, "abundTF"=0)

# Make inputs
Inputs <- makeInput(family = family, df = df_yoy, c_ip = c_ip_yoy, options = Options_vec, X = X_ij, t_i = t_i, version = Version, offset_i = offset, CalcSD_lambda_ip = Calc_lambda_ip)

Mod4r3 <- runOUGMRF(Inputs)
str(Mod4r3) 
Mod4r3$opt$convergence
Mod4r3$SD$sd

#----------------- Temporal + ST ------------------
X_ij <- as.matrix(dplyr::select(covs, forest, temp_mean_winter, temp_mean_spring, prcp_mean_winter, prcp_mean_spring))

# Turn off random effects in v1f (0 means exclude a component, except for ObsModel)
Options_vec = c("SpatialTF"=0, "TemporalTF"=1, "SpatiotemporalTF"=1, "DetectabilityTF"=1, "ObsModel"=1, "OverdispersedTF" = 1)

# Make inputs
Inputs <- makeInput(family = family, df = df_yoy, c_ip = c_ip_yoy, options = Options_vec, X = X_ij, t_i = t_i, version = Version, offset_i = offset, CalcSD_lambda_ip = Calc_lambda_ip)

Mod5r <- runOUGMRF(Inputs)
Mod5r$SD$sd # Fail

# Further Reduce
X_ij <- as.matrix(dplyr::select(covs, forest, temp_mean_spring, prcp_mean_winter))

# Make inputs
Inputs <- makeInput(family = family, df = df_yoy, c_ip = c_ip_yoy, options = Options_vec, X = X_ij, t_i = t_i, version = Version, offset_i = offset, CalcSD_lambda_ip = Calc_lambda_ip)

Mod5r2 <- runOUGMRF(Inputs)
str(Mod5r2)
Mod5r2$opt$convergence
Mod5r2$SD$sd # Fails

# Final Reduction
X_ij <- as.matrix(dplyr::select(covs, forest))

# Make inputs
Inputs <- makeInput(family = family, df = df_yoy, c_ip = c_ip_yoy, options = Options_vec, X = X_ij, t_i = t_i, version = Version, offset_i = offset, CalcSD_lambda_ip = Calc_lambda_ip)

Mod5r3 <- runOUGMRF(Inputs)
str(Mod5r3)
Mod5r3$opt$convergence
Mod5r3$SD$sd # Fails


#----------------- Spatial Temporal Spatiotemporal Reduced ------------------
# Turn off random effects in v1f (0 means exclude a component, except for ObsModel)
Options_vec = c("SpatialTF"=1, "TemporalTF"=1, "SpatiotemporalTF"=1, "DetectabilityTF"=1, "ObsModel"=1, "OverdispersedTF"=1, "abundTF"=0)

# Final Reduction
X_ij <- as.matrix(dplyr::select(covs, forest))

# Make inputs
Inputs <- makeInput(family = family, df = df_yoy, c_ip = c_ip_yoy, options = Options_vec, X = X_ij, t_i = t_i, version = Version, offset_i = offset, CalcSD_lambda_ip = Calc_lambda_ip)

Mod6r3 <- runOUGMRF(Inputs)
str(Mod6r3)
Mod6r3$opt$convergence
Mod6r3$SD$sd # Fails

#----------------- Spatial + Temporal ------------------
X_ij <- as.matrix(dplyr::select(covs, forest, temp_mean_winter, temp_mean_spring, prcp_mean_winter, prcp_mean_spring))

# Turn off random effects in v1f (0 means exclude a component, except for ObsModel)
Options_vec = c("SpatialTF"=1, "TemporalTF"=1, "SpatiotemporalTF"=0, "DetectabilityTF"=1, "ObsModel"=1, "OverdispersedTF" = 1)

# Make inputs
Inputs <- makeInput(family = family, df = df_yoy, c_ip = c_ip_yoy, options = Options_vec, X = X_ij, t_i = t_i, version = Version, offset_i = offset, CalcSD_lambda_ip = Calc_lambda_ip)

Mod7r <- runOUGMRF(Inputs)
Mod7r$SD$sd # success

# Further Reduce
X_ij <- as.matrix(dplyr::select(covs, forest, temp_mean_spring, prcp_mean_winter))

# Make inputs
Inputs <- makeInput(family = family, df = df_yoy, c_ip = c_ip_yoy, options = Options_vec, X = X_ij, t_i = t_i, version = Version, offset_i = offset, CalcSD_lambda_ip = Calc_lambda_ip)

Mod7r2 <- runOUGMRF(Inputs)
str(Mod7r2)
Mod7r2$opt$convergence
Mod7r2$SD$sd # success

# final reduction - no time-varying covs
X_ij <- as.matrix(dplyr::select(covs, forest))

# Make inputs
Inputs <- makeInput(family = family, df = df_yoy, c_ip = c_ip_yoy, options = Options_vec, X = X_ij, t_i = t_i, version = Version, offset_i = offset, CalcSD_lambda_ip = Calc_lambda_ip)

Mod7r3 <- runOUGMRF(Inputs)
str(Mod7r3)
Mod7r3$opt$convergence
Mod7r3$SD$sd 

#----------------- Spatial + Spatiotemporal ------------------
X_ij <- as.matrix(dplyr::select(covs, forest, temp_mean_winter, temp_mean_spring, prcp_mean_winter, prcp_mean_spring))

Options_vec = c("SpatialTF"=1, "TemporalTF"=0, "SpatiotemporalTF"=1, "DetectabilityTF"=1, "ObsModel"=1, "OverdispersedTF"=1, "abundTF"=0)

# Make inputs
Inputs <- makeInput(family = family, df = df_yoy, c_ip = c_ip_yoy, options = Options_vec, X = X_ij, t_i = t_i, version = Version, offset_i = offset, CalcSD_lambda_ip = Calc_lambda_ip)

Mod8r <- runOUGMRF(Inputs)
str(Mod8r) # fails

# Further Reduce
X_ij <- as.matrix(dplyr::select(covs, forest, temp_mean_spring, prcp_mean_winter))

# Make inputs
Inputs <- makeInput(family = family, df = df_yoy, c_ip = c_ip_yoy, options = Options_vec, X = X_ij, t_i = t_i, version = Version, offset_i = offset, CalcSD_lambda_ip = Calc_lambda_ip)

Mod8r2 <- runOUGMRF(Inputs)
str(Mod8r2)
Mod8r2$opt$convergence
Mod8r2$SD$sd # success

# final reduction - no time-varying covs
X_ij <- as.matrix(dplyr::select(covs, forest))

# Make inputs
Inputs <- makeInput(family = family, df = df_yoy, c_ip = c_ip_yoy, options = Options_vec, X = X_ij, t_i = t_i, version = Version, offset_i = offset, CalcSD_lambda_ip = Calc_lambda_ip)

Mod8r3 <- runOUGMRF(Inputs)
str(Mod8r3)
Mod8r3$opt$convergence
Mod8r3$SD$sd 

#--------------------------------------------


# AIC
Model <- c("Obs", 
           "Temporal", 
           "Spatial",
           "Spatial reduced",
           "Spatial reduced 2",
           "Spatiotemporal",
           "Spatiotemporal reduced",
           # "Spatiotemporal reduced 2",
           "Spatiotemporal reduced 3",
            "Temporal + ST", 
           # "S+T+ST", 
           "Spatial + Temporal", 
           "Spatial + Temporal reduced", 
           "Spatial + ST reduced 3",
           "Spatial + ST reduced 3"
) #
M_num <- c("1",
           "2",
           "3",
           "3r",
           "3r2",
           "4",
           "4r",
          # "4r2",
           "4r3",
           "5",
           #6,
           "7",
           "7r",
           "8r2",
           "8r3"
)
AIC <- c(opt1b$AIC, 
             opt2b$AIC, 
             opt3b$AIC,
             Mod3r$opt$AIC, #, 
             Mod3r2$opt$AIC,
         opt4b$AIC,
             Mod4r$opt$AIC,
            # Mod4r2$opt$AIC,
             Mod4r3$opt$AIC,
         opt5b$AIC,
             opt7b$AIC,
             Mod7r$opt$AIC, 
             Mod8r2$opt$AIC,
             Mod8r3$opt$AIC
) # 
aic_table <- data.frame(M_num, Model, AIC, stringsAsFactors = FALSE)
names(aic_table) <- c("M_num", "Model", "AIC")
aic_table <- dplyr::arrange(aic_table, AIC)
aic_table$delta_AIC <- 0
for(i in 2:nrow(aic_table)) {
  aic_table$delta_AIC[i] <- aic_table$AIC[i] - aic_table$AIC[1]
}
aic_table

# no additional spatial autocorrelation - maybe some temporal

# compare models with reduced covariates for best models and best model with no overdispersion

#----------- reduced Obs models w/no spatial ----------------
X_ij <- as.matrix(dplyr::select(covs, forest, temp_mean_winter, temp_mean_spring, prcp_mean_winter, prcp_mean_spring))

Options_vec = c("SpatialTF"=0, "TemporalTF"=0, "SpatiotemporalTF"=0, "DetectabilityTF"=1, "ObsModel"=1, "OverdispersedTF"=1, "abundTF"=0)

# Make inputs
Inputs <- makeInput(family = family, df = df_yoy, c_ip = c_ip_yoy, options = Options_vec, X = X_ij, t_i = t_i, version = Version, offset_i = offset, CalcSD_lambda_ip = Calc_lambda_ip)

Mod1r <- runOUGMRF(Inputs)
Mod1r$SD$sd

# Further Reduce
X_ij <- as.matrix(dplyr::select(covs, forest, temp_mean_spring, prcp_mean_winter))

# Make inputs
Inputs <- makeInput(family = family, df = df_yoy, c_ip = c_ip_yoy, options = Options_vec, X = X_ij, t_i = t_i, version = Version, offset_i = offset, CalcSD_lambda_ip = Calc_lambda_ip)

Mod1r2 <- runOUGMRF(Inputs)
Mod1r2$SD$sd

# final reduction - no time-varying covs
X_ij <- as.matrix(dplyr::select(covs, forest))

# Make inputs
Inputs <- makeInput(family = family, df = df_yoy, c_ip = c_ip_yoy, options = Options_vec, X = X_ij, t_i = t_i, version = Version, offset_i = offset, CalcSD_lambda_ip = Calc_lambda_ip)

Mod1r3 <- runOUGMRF(Inputs)
Mod1r2$SD$sd


#----------- reduced temporal AR1 models w/no spatial ----------------
X_ij <- as.matrix(dplyr::select(covs, forest, temp_mean_winter, temp_mean_spring, prcp_mean_winter, prcp_mean_spring))

Options_vec = c("SpatialTF"=0, "TemporalTF"=1, "SpatiotemporalTF"=0, "DetectabilityTF"=1, "ObsModel"=1, "OverdispersedTF"=1, "abundTF"=0)

# Make inputs
Inputs <- makeInput(family = family, df = df_yoy, c_ip = c_ip_yoy, options = Options_vec, X = X_ij, t_i = t_i, version = Version, offset_i = offset, CalcSD_lambda_ip = Calc_lambda_ip)

Mod2r <- runOUGMRF(Inputs)
Mod2r$SD$sd

# Further Reduce
X_ij <- as.matrix(dplyr::select(covs, forest, temp_mean_spring, prcp_mean_winter))

# Make inputs
Inputs <- makeInput(family = family, df = df_yoy, c_ip = c_ip_yoy, options = Options_vec, X = X_ij, t_i = t_i, version = Version, offset_i = offset, CalcSD_lambda_ip = Calc_lambda_ip)

Mod2r2 <- runOUGMRF(Inputs)
Mod2r2$SD$sd

# final reduction - no time-varying covs
X_ij <- as.matrix(dplyr::select(covs, forest))

# Make inputs
Inputs <- makeInput(family = family, df = df_yoy, c_ip = c_ip_yoy, options = Options_vec, X = X_ij, t_i = t_i, version = Version, offset_i = offset, CalcSD_lambda_ip = Calc_lambda_ip)

Mod2r3 <- runOUGMRF(Inputs)
Mod2r3$SD$sd

#---------------------------

# AIC
Model <- c("Obs",
           "Obs reduced",
           "Obs reduced 2",
           "Obs reduced 3",
           "Temporal",
           "Temporal reduced",
           "Temporal reduced 2",
           "Temporal reduced 3",
           "Spatiotemporal",
           "Spatiotemporal reduced",
           #"Spatiotemporal reduced 2",
           "Spatiotemporal reduced 3",
           #"Spatial",#, 
           # "Spatiotemporal", 
           # "Temporal + ST", 
           # "S+T+ST", 
           "Spatial + Temporal",
           "Spatial + Temporal reduced",
           "Spatial + Temporal reduced 2",
           "Spatial + Temporal reduced 3"
           # "Spatial + ST"
) #
M_num <- c("1",
           "1r",
           "1r2",
           "1r3",
           "2",
           "2r",
           "2r2",
           "2r3",
           "4",
           "4r",
          # "4r2",
           "4r3",
           "7",
           "7r",
           "7r2",
           "7r3"
)
#k <- c(Mod1r$opt$k)
AIC <- c(opt1b$AIC, 
         Mod1r$opt$AIC, 
         Mod1r2$opt$AIC,
         Mod1r3$opt$AIC,
         opt2b$AIC,
         Mod2r$opt$AIC,
         Mod2r2$opt$AIC,
         Mod2r3$opt$AIC,
         opt4b$AIC,
         Mod4r$opt$AIC,
         #Mod4r2$opt$AIC,
         Mod4r3$opt$AIC,
         opt7b$AIC,
         Mod7r$opt$AIC,
         Mod7r2$opt$AIC,
         Mod7r3$opt$AIC
) # 
aic_table2 <- data.frame(M_num, AIC, stringsAsFactors = FALSE)
names(aic_table2) <- c("M_num", "AIC")
aic_table2 <- dplyr::arrange(aic_table2, AIC)
aic_table2$delta_AIC <- 0
for(i in 2:nrow(aic_table2)) {
  aic_table2$delta_AIC[i] <- aic_table2$AIC[i] - aic_table2$AIC[1]
}
aic_table2

# better to keep all fixed effects
X_ij <- as.matrix(dplyr::select(covs, forest, surfcoarse, temp_mean_fall_1, temp_mean_winter, temp_mean_spring, prcp_mean_fall_1, prcp_mean_winter, prcp_mean_spring))
Parameters = colnames(as.matrix(X_ij))

X_r2 <- as.matrix(dplyr::select(covs, forest, temp_mean_spring, prcp_mean_winter))
Parameters_r2 <- colnames(as.matrix(X_r2))


summary(SD7b, "fixed", p.value = TRUE)
SD7b$par.fixed

Mod5 <- list(Report = Report5b, opt = opt5b, SD = SD5b)

save.image("Output/W_Susquehanna_YOY.RData")
save(Mod5, SD_table, aic_table, aic_table2, covs, X_ij, Parameters, df = df_yoy, family, file = "Output/W_Susquehanna_YOY_Summary.RData")

saveRDS(list(df_yoy = df_yoy,
        family = family,
        Report5=Mod5$Report, 
        opt5=Mod5$opt, 
        SD5=Mod5$SD, 
        SD_table=SD_table, 
        Parameters = Parameters,
        covs = covs,
        aic_table=aic_table, 
        aic_table2 = aic_table2),
        file = "Output/W_Susquehanna_YOY_Summary_RDS.RData")






