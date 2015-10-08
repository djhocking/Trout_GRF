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

opt1_nelder_mead <- optim(obj1$env$last.par.best[-c(obj1$env$random)], fn = obj1$fn, gr = obj1$gr)
opt1_BFGS <- optim(obj1$env$last.par.best[-c(obj1$env$random)], fn = obj1$fn, gr = obj1$gr, method = "BFGS")

cbind(NLMINB = c(opt1$par, opt1$convergence), nelder_mead = c(opt1_nelder_mead$par, opt1_nelder_mead$convergence), BFGS = c(opt1_BFGS$par, opt1_BFGS$convergence))

SD1_BFGS <- sdreport(obj1)
data.frame(SD1$sd, SD1_BFGS$sd)
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

# # Diagnose problems of convergence and SD estimation
# DiagnosticDir <- "Diagnostics/"
# # create code directory if doesn't exist
# if (!file.exists(DiagnosticDir)) {
#   dir.create(DiagnosticDir)
# }
# 
# obj4$gr_orig = obj4$gr
# obj4$fn_orig = obj4$fn
# obj4$fn <- function( vec ){
#   Fn = obj4$fn_orig(vec)
#   if( any(is.na(Fn ))) capture.output( matrix(Fn,ncol=1,dimnames=list(names(obj4$par),NULL)), file=paste0(DiagnosticDir,"Fn4.txt"), append = TRUE )
#   return( Fn )
# }
# obj4$gr = function( vec ){
#   Gr = obj4$gr_orig(vec)
#   if( any(is.na(Gr))) capture.output( matrix(Gr,ncol=1,dimnames=list(names(obj4$par),NULL)), file=paste0(DiagnosticDir,"gr4.txt"), append = TRUE )
#   return( Gr )
# }


# First run
obj4$fn( obj4$par )
# Check for parameters that don't do anything
Which = which( obj4$gr( obj4$par )==0 )

# Run model (nlminb is slow)
time1 <- Sys.time()
opt4 = nlminb(start=obj4$env$last.par.best[-c(obj4$env$random)], objective=obj4$fn, gradient=obj4$gr, control=list(eval.max=1e4, iter.max=1e4, trace=1, rel.tol=1e-14) )
opt4[["final_gradient"]] = obj4$gr( opt4$par )
opt4[["AIC"]] = 2*opt4$objective + 2*length(opt4$par)

Report4 = obj4$report()
SD4 = sdreport( obj4, bias.correct=FALSE )
time2 <- Sys.time()
time.nlminb <- time2 - time1

# BFGS - fast
time1 <- Sys.time()
opt4_BFGS <- optim(obj4$env$last.par.best[-c(obj4$env$random)], fn = obj4$fn, gr = obj4$gr, method = "BFGS")
SD4_BFGS <- sdreport(obj4) # fails
time2 <- Sys.time()
time.bfgs <- time2 - time1

# NM = slow!!!!!!!
time1 <- Sys.time()
opt4_NM <- optim(obj4$env$last.par.best[-c(obj4$env$random)], fn = obj4$fn, gr = obj4$gr, method = "Nelder-Mead")
SD4_NM <- sdreport(obj4)
time2 <- Sys.time()
time.nm <- time2 - time1

# CG
# time1 <- Sys.time()
# opt4_CG <- optim(obj4$env$last.par.best[-c(obj4$env$random)], fn = obj4$fn, gr = obj4$gr, method = "CG")
# SD4_CG <- sdreport(obj4)
# time2 <- Sys.time()
# time.cg <- time2 - time1

# author of CG, John Nash, says not good and should Rcgmin instead
library(Rcgmin)
time1 <- Sys.time()
opt4_CG <- Rcgmin(obj4$env$last.par.best[-c(obj4$env$random)], fn = obj4$fn, gr = obj4$gr)
SD4_CG <- sdreport(obj4)
time2 <- Sys.time()
time.cg <- time2 - time1

# LBFGSB - not sure if better than BFGS if don't set bounds: aparently yes
time1 <- Sys.time()
opt4_LBFGSB <- optim(obj4$env$last.par.best[-c(obj4$env$random)], fn = obj4$fn, gr = obj4$gr, method = "L-BFGS-B")
SD4_LBFGSB <- sdreport(obj4)
time2 <- Sys.time()
time.lbfgsb <- time2 - time1

# bobyqa
time1 <- Sys.time()
opt4_bobyqa <- bobyqa(obj4$env$last.par.best[-c(obj4$env$random)], fn = obj4$fn)
SD4_bobyqa <- sdreport(obj4)
time2 <- Sys.time()
time.bobyqa <- time2 - time1

# SANN will be super slow - only run as a LAST resort
# time1 <- Sys.time()
# opt4_SANN <- optim(obj4$env$last.par.best[-c(obj4$env$random)], fn = obj4$fn, gr = obj4$gr, method = "SANN")
# SD4_SANN<- sdreport(obj4)
# time2 <- Sys.time()
# time.sann <- time2 - time1

cbind(NLMINB = c(opt4$par, opt4$convergence), 
      BFGS = c(opt4_BFGS$par, opt4_BFGS$convergence),
      CG = c(opt4_CG$par, opt4_CG$convergence),
      LBFGSB = c(opt4_LBFGSB$par, opt4_LBFGSB$convergence),
      NM = c(opt4_NM$par, opt4_NM$convergence),
      BOBYQA = c(opt4_bobyqa$par, opt4_bobyqa$ierr))
data.frame(SD4$sd, 
           SD4_NM$sd,
           #SD4_BFGS$sd,
           SD4_CG$sd,
           SD4_bobyqa$sd,
           SD4_LBFGSB$sd)
data.frame(optimizer = c("nlmninb",
             "nm",
             "bfgs",
             "cg",
             "bobyqa",
             "lbfgsb"),
           time = c(time.nlminb,
           time.nm*60,
           time.bfgs,
           time.cg,
           time.bobyqa,
           time.lbfgsb), 
           SD_success = c(FALSE,
                       TRUE,
                       FALSE,
                       TRUE,
                       TRUE,
                       TRUE),
           convergence = c(opt4$convergence,
                           opt4_NM$convergence,
                           opt4_BFGS$convergence,
                           opt4_CG$convergence,
                           opt4_bobyqa$ierr,
                           opt4_LBFGSB$convergence))

# BOBYQA and LBFGSB are the best options here

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

# Diagnose problems of convergence and SD estimation
# DiagnosticDir <- "Diagnostics/"
# # create code directory if doesn't exist
# if (!file.exists(DiagnosticDir)) {
#   dir.create(DiagnosticDir)
# }

# obj5$gr_orig = obj5$gr
# obj5$fn_orig = obj5$fn
# obj5$fn <- function( vec ){
#   Fn = obj5$fn_orig(vec)
#   capture.output( as.matrix((t(c(Fn = Fn, obj5$env$last.par[-c(obj5$env$random)])))), file = paste0(DiagnosticDir, "Fn5.txt"), append = TRUE) # 
#  # if( any(is.na(Fn ))) capture.output( matrix(Fn,ncol=1,dimnames=list(names(obj5$par),NULL)), file=paste0(DiagnosticDir,"Fn5.txt"), append = FALSE )
#   return( Fn )
# }
# obj5$gr = function( vec ){
#   Gr = obj5$gr_orig(vec)
#  # if( any(is.na(Gr))) capture.output( matrix(Gr,ncol=1,dimnames=list(names(obj5$par),NULL)), file=paste0(DiagnosticDir,"gr5.txt"), append = FALSE )
#   return( Gr )
# }

# First run
obj5$fn( obj5$par )
#fn_test < obj5$fn(obj5$par)

# Check for parameters that don't do anything
Which = which( obj5$gr( obj5$par )==0 )

# Run model
opt5 = nlminb(start=obj5$env$last.par.best[-c(obj5$env$random)], objective=obj5$fn, gradient=obj5$gr, control=list(eval.max=1e4, iter.max=1e4, trace=1, rel.tol=1e-14) )
opt5[["final_gradient"]] = obj5$gr( opt5$par )
opt5[["AIC"]] = 2*opt5$objective + 2*length(opt5$par)

ParHat <- obj5( opt5$par )
Report5 = obj5$report()
SD5 = sdreport( obj5, bias.correct=F )

# LBFGSB
time1.lbfgsb.5 <- Sys.time()
opt5_LBFGSB <- optim(obj5$env$last.par.best[-c(obj5$env$random)], fn = obj5$fn, gr = obj5$gr, method = "L-BFGS-B")
SD5_LBFGSB <- sdreport(obj5)
time2.lbfgsb.5 <- Sys.time()
time.lbfgsb.5 <- time2.lbfgsb.5 - time1.lbfgsb.5

Report5.lbfgsb = obj5$report()
SD5.lbfgsb = sdreport( obj5, bias.correct=F ) # 

# bobyqa
time1_bobyqa_5 <- Sys.time()
opt5_bobyqa <- bobyqa(obj5$env$last.par.best[-c(obj5$env$random)], fn = obj5$fn)
SD5_bobyqa <- sdreport(obj5)
time2_bobyqa_5 <- Sys.time()
time_bobyqa_5 <- time2_bobyqa_5 - time1_bobyqa_5

Report5.bobyqa = obj5$report()
SD5.bobyqa = sdreport( obj5, bias.correct=F ) # 
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

# Diagnose problems of convergence and SD estimation
DiagnosticDir <- "Diagnostics/"
# create code directory if doesn't exist
if (!file.exists(DiagnosticDir)) {
  dir.create(DiagnosticDir)
}

# obj6$gr_orig = obj6$gr
# obj6$fn_orig = obj6$fn
# obj6$fn <- function( vec ){
#   Fn = obj6$fn_orig(vec)
#   if( any(is.na(Fn ))) capture.output( matrix(Fn,ncol=1,dimnames=list(names(obj6$par),NULL)), file=paste0(DiagnosticDir,"Fn6.txt"), append = TRUE )
#   return( Fn )
# }
# obj6$gr = function( vec ){
#   Gr = obj6$gr_orig(vec)
#   if( any(is.na(Gr))) capture.output( matrix(Gr,ncol=1,dimnames=list(names(obj6$par),NULL)), file=paste0(DiagnosticDir,"gr6.txt"), append = TRUE )
#   return( Gr )
# }

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

# bobyqa
time1_bobyqa_6 <- Sys.time()
opt6_bobyqa <- bobyqa(obj6$env$last.par.best[-c(obj6$env$random)], fn = obj6$fn)
SD6_bobyqa <- sdreport(obj6)
time2_bobyqa_6 <- Sys.time()
time_bobyqa_6 <- time2_bobyqa_6 - time1_bobyqa_6

Report6_bobyqa = obj6$report()

data.frame(names(SD6$value), SD6$sd, SD6_bobyqa$sd)
opt6_bobyqa[["final_gradient"]] = obj6$gr( opt6_bobyqa$par )
opt6_bobyqa$ierr # converges but SD fails
opt6_bobyqa[["AIC"]] = 2*opt6_bobyqa$fval + 2*length(opt6_bobyqa$par)

# model can't differentiate between spatial and spatiotemporal

# LBFGSB
lower6g <- c(log(0), log(1e-4), log(0), log(0), 0, log(1e-4), log(1e-4), -100, log(0), -10, -10, -10, -10, log(0.01))
upper6g <- c(log(800), log(100), log(800), log(100), 100, log(100), log(100), 100, log(100), 10, 10, 10, 10, log(10))

time1.lbfgsb.6 <- Sys.time()
opt6_LBFGSB <- optim(obj6$env$last.par.best[-c(obj6$env$random)], fn = obj6$fn, gr = obj6$gr, method = "L-BFGS-B", lower = lower6g, upper = upper6g)
SD6_LBFGSB <- sdreport(obj6)
time2.lbfgsb.6 <- Sys.time()
time.lbfgsb.6 <- time2.lbfgsb.6 - time1.lbfgsb.6

Report6.lbfgsb = obj6$report()
SD6.lbfgsb = sdreport( obj6, bias.correct=F )

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

# Diagnose problems of convergence and SD estimation
# DiagnosticDir <- "Diagnostics/"
# # create code directory if doesn't exist
# if (!file.exists(DiagnosticDir)) {
#   dir.create(DiagnosticDir)
# }
# 
# obj7$gr_orig = obj7$gr
# obj7$fn_orig = obj7$fn
# obj7$fn <- function( vec ){
#   Fn = obj7$fn_orig(vec)
#   if( any(is.na(Fn ))) capture.output( matrix(Fn,ncol=1,dimnames=list(names(obj7$par),NULL)), file=paste0(DiagnosticDir,"Fn.txt"), append = TRUE )
#   return( Fn )
# }
# obj7$gr = function( vec ){
#   Gr = obj7$gr_orig(vec)
#   if( any(is.na(Gr))) capture.output( matrix(Gr,ncol=1,dimnames=list(names(obj7$par),NULL)), file=paste0(DiagnosticDir,"gr.txt"), append = TRUE )
#   return( Gr )
# }

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

# bobyqa
time1_bobyqa_7 <- Sys.time()
opt7_bobyqa <- bobyqa(obj7$env$last.par.best[-c(obj7$env$random)], fn = obj7$fn)
SD7_bobyqa <- sdreport(obj7)
time2_bobyqa_7 <- Sys.time()
time_bobyqa_7 <- time2_bobyqa_7 - time1_bobyqa_7

Report7_bobyqa = obj7$report()

cbind(SD7$sd, SD7_bobyqa$sd)
opt7_bobyqa[["final_gradient"]] = obj7$gr( opt7_bobyqa$par )
opt7_bobyqa$ierr # converges
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
# ports
opt8 = nlminb(start=obj8$env$last.par.best[-c(obj8$env$random)], objective=obj8$fn, gradient=obj8$gr, control=list(eval.max=1e4, iter.max=1e4, trace=1, rel.tol=1e-14) )
opt8[["final_gradient"]] = obj8$gr( opt8$par )
opt8[["AIC"]] = 2*opt8$objective + 2*length(opt8$par)

ParHat <- obj8( opt8$par )
Report8 = obj8$report()
SD8 = sdreport( obj8, bias.correct=F )

# LBFGSB
time1.lbfgsb.8 <- Sys.time()
opt8_LBFGSB <- optim(obj8$env$last.par.best[-c(obj8$env$random)], fn = obj8$fn, gr = obj8$gr, method = "L-BFGS-B")
SD8_LBFGSB <- sdreport(obj8)
time2.lbfgsb.8 <- Sys.time()
time.lbfgsb.8 <- time2.lbfgsb.8 - time1.lbfgsb.8

Report8.lbfgsb = obj8$report()
SD8.lbfgsb = sdreport( obj8, bias.correct=F ) # 

# bobyqa
time1_bobyqa_8 <- Sys.time()
opt8_bobyqa <- bobyqa(obj8$env$last.par.best[-c(obj8$env$random)], fn = obj8$fn)
SD8_bobyqa <- sdreport(obj8)
time2_bobyqa_8 <- Sys.time()
time_bobyqa_8 <- time2_bobyqa_8 - time1_bobyqa_8

Report8.bobyqa = obj8$report()
SD8.bobyqa = sdreport( obj8, bias.correct=F ) # 
#--------------------------------------------------

################# Try without IID variation in lambda ##########
Version = "OU_GMRF_v1h"

# Compile
if(FALSE) {
  dyn.unload(dynlib(paste0("Code/", Version)))
  file.remove( paste0("Code/", Version,c(".o",".dll")) )
}
compile( paste0("Code/", Version,".cpp") )

#-------------- S + T + ST ----------------
# Turn off random effects in v1f (0 means exclude a component, except for ObsModel)
Options_vec = c("SpatialTF"=1, "TemporalTF"=1, "SpatiotemporalTF"=1, "DetectabilityTF"=1, "ObsModel"=1, "OverdispersedTF"=0)

# Make inputs
Inputs <- makeInput(family = family, c_ip = c_ip, options = Options_vec, X = X_ij, t_i = t_i, version = Version)

# Make object
dyn.load( dynlib(paste0("Code/", Version )))
obj6h <- MakeADFun(data=Inputs$Data, parameters=Inputs$Params, random=Inputs$Random, map=Inputs$Map, hessian=FALSE, inner.control=list(maxit=1000) )
Report = obj6h$report()

# First run
obj6h$fn( obj6h$par )
# Check for parameters that don't do anything
Which = which( obj6h$gr( obj6h$par )==0 )

# Run model
# bobyqa
library(optimx)
lower6h <- c(log(0), log(1e-4), log(0), log(0), 0, log(1e-4), log(1e-4), -100, log(0), -10, -10, -10, -10, log(0.01))
upper6h <- c(log(800), log(100), log(800), log(100), 100, log(100), log(100), 100, log(100), 10, 10, 10, 10, log(10))

time1_bobyqa_6h <- Sys.time()
opt6h_bobyqa <- bobyqa(obj6h$env$last.par.best[-c(obj6h$env$random)], fn = obj6h$fn) # optimx(obj6h$par, fn = obj6h$fn, gradient=obj6h$gr, method = "bobyqa", lower = lower6h, upper = upper6h)
SD6h_bobyqa <- sdreport(obj6h)
time2_bobyqa_6h <- Sys.time()
time_bobyqa_6h <- time2_bobyqa_6h - time1_bobyqa_6h

Report6h_bobyqa = obj6h$report()

data.frame(names(SD6h_bobyqa$value), SD6h_bobyqa$sd)
opt6h_bobyqa[["final_gradient"]] = obj6h$gr( opt6h_bobyqa$par )
opt6h_bobyqa$ierr # 
opt6h_bobyqa[["AIC"]] = 2*opt6h_bobyqa$fval + 2*length(opt6h_bobyqa$par)

# LBFGSB
time1.lbfgsb.6h <- Sys.time()
opt6h_LBFGSB <- optimx(obj6h$par, fn = obj6h$fn, gradient=obj6h$gr, method = "L-BFGS-B", lower = lower6h, upper = upper6h)
SD6h_LBFGSB <- sdreport(obj6h) # fails with NaN
time2.lbfgsb.6h <- Sys.time()
time.lbfgsb.6h <- time2.lbfgsb.6h - time1.lbfgsb.6h # quick = 7 min

Report6h.lbfgsb = obj6h$report()


#-------------- S + T no overdispersion ----------------
# Turn off random effects in v1f (0 means exclude a component, except for ObsModel)
Options_vec = c("SpatialTF"=1, "TemporalTF"=1, "SpatiotemporalTF"=0, "DetectabilityTF"=1, "ObsModel"=1, "OverdispersedTF"=0)

# Make inputs
Inputs <- makeInput(family = family, c_ip = c_ip, options = Options_vec, X = X_ij, t_i = t_i, version = Version)

# Make object
dyn.load( dynlib(paste0("Code/", Version )))
obj8 <- MakeADFun(data=Inputs$Data, parameters=Inputs$Params, random=Inputs$Random, map=Inputs$Map, hessian=FALSE, inner.control=list(maxit=1000) )
Report = obj8$report()

# First run
obj8$fn( obj8$par )
# Check for parameters that don't do anything
Which = which( obj8$gr( obj8$par )==0 )

# Run model
# bobyqa
library(optimx)
#lower8 <- c(log(0), log(1e-4), log(0), log(0), 0, log(1e-4), log(1e-4), -100, log(0), -10, -10, -10, -10, log(0.01))
#upper8 <- c(log(800), log(100), log(800), log(100), 100, log(100), log(100), 100, log(100), 10, 10, 10, 10, log(10))

time1_bobyqa_8 <- Sys.time()
opt8_bobyqa <- bobyqa(obj8$env$last.par.best[-c(obj8$env$random)], fn = obj8$fn) # optimx(obj8$par, fn = obj8$fn, gradient=obj8$gr, method = "bobyqa", lower = lower8, upper = upper8)
SD8_bobyqa <- sdreport(obj8) # fails
time2_bobyqa_8 <- Sys.time()
time_bobyqa_8 <- time2_bobyqa_8 - time1_bobyqa_8

Report8_bobyqa = obj8$report()

data.frame(names(SD8_bobyqa$value), SD8_bobyqa$sd)
opt8_bobyqa[["final_gradient"]] = obj8$gr( opt8_bobyqa$par )
opt8_bobyqa$ierr # 
opt8_bobyqa[["AIC"]] = 2*opt8_bobyqa$fval + 2*length(opt8_bobyqa$par)

# LBFGSB
time1.lbfgsb.8 <- Sys.time()
opt8_LBFGSB <- optimx(obj8$par, fn = obj8$fn, gradient=obj8$gr, method = "L-BFGS-B")
SD8_LBFGSB <- sdreport(obj8) # fails with NaN
time2.lbfgsb.8 <- Sys.time()
time.lbfgsb.8 <- time2.lbfgsb.8 - time1.lbfgsb.8 # quick = 7 min

Report8.lbfgsb = obj8$report()

# nm
time1.nm.8 <- Sys.time()
opt8_nm <- optimx(obj8$par, fn = obj8$fn, gradient=obj8$gr, method = "Nelder-Mead")
SD8_nm <- sdreport(obj8) # fails with NaN
time2.nm.8 <- Sys.time()
time.nm.8 <- time2.nm.8 - time1.nm.8 # 

Report8.nm = obj8$report()

#-------------- ST ----------------
# Turn off random effects in v1f (0 means exclude a component, except for ObsModel)
Options_vec = c("SpatialTF"=0, "TemporalTF"=0, "SpatiotemporalTF"=1, "DetectabilityTF"=1, "ObsModel"=1, "OverdispersedTF"=0)

# Make inputs
Inputs <- makeInput(family = family, c_ip = c_ip, options = Options_vec, X = X_ij, t_i = t_i, version = Version)

# Make object
dyn.load( dynlib(paste0("Code/", Version )))
obj9 <- MakeADFun(data=Inputs$Data, parameters=Inputs$Params, random=Inputs$Random, map=Inputs$Map, hessian=FALSE, inner.control=list(maxit=1000) )
Report = obj9$report()

# First run
obj9$fn( obj9$par )
# Check for parameters that don't do anything
Which = which( obj9$gr( obj9$par )==0 )

# Run model
# bobyqa
library(optimx)
lower9 <- c(log(0), log(1e-4), log(0), log(0), 0, log(1e-4), log(1e-4), -100, log(0), -10, -10, -10, -10, log(0.01))
upper9 <- c(log(800), log(100), log(800), log(100), 100, log(100), log(100), 100, log(100), 10, 10, 10, 10, log(10))

time1_bobyqa_9 <- Sys.time()
opt9_bobyqa <- bobyqa(obj9$env$last.par.best[-c(obj9$env$random)], fn = obj9$fn) # optimx(obj9$par, fn = obj9$fn, gradient=obj9$gr, method = "bobyqa", lower = lower9, upper = upper9)
SD9_bobyqa <- sdreport(obj9) # fails
time2_bobyqa_9 <- Sys.time()
time_bobyqa_9 <- time2_bobyqa_9 - time1_bobyqa_9

Report9_bobyqa = obj9$report()

data.frame(names(SD9_bobyqa$value), SD9_bobyqa$sd)
opt9_bobyqa[["final_gradient"]] = obj9$gr( opt9_bobyqa$par )
opt9_bobyqa$ierr # 
opt9_bobyqa[["AIC"]] = 2*opt9_bobyqa$fval + 2*length(opt9_bobyqa$par)

#-------------- Obs ----------------
# Turn off random effects in v1f (0 means exclude a component, except for ObsModel)
Options_vec = c("SpatialTF"=0, "TemporalTF"=0, "SpatiotemporalTF"=0, "DetectabilityTF"=1, "ObsModel"=1, "OverdispersedTF"=0)

# Make inputs
Inputs <- makeInput(family = family, c_ip = c_ip, options = Options_vec, X = X_ij, t_i = t_i, version = Version)

# Make object
dyn.load( dynlib(paste0("Code/", Version )))
obj10 <- MakeADFun(data=Inputs$Data, parameters=Inputs$Params, random=Inputs$Random, map=Inputs$Map, hessian=FALSE, inner.control=list(maxit=1000) )
Report = obj10$report()

# First run
obj10$fn( obj10$par )
# Check for parameters that don't do anything
Which = which( obj10$gr( obj10$par )==0 )

# Run model
# bobyqa
library(optimx)
lower10 <- c(log(0), log(1e-4), log(0), log(0), 0, log(1e-4), log(1e-4), -100, log(0), -10, -10, -10, -10, log(0.01))
upper10 <- c(log(800), log(100), log(800), log(100), 100, log(100), log(100), 100, log(100), 10, 10, 10, 10, log(10))

time1_bobyqa_10 <- Sys.time()
opt10_bobyqa <- bobyqa(obj10$env$last.par.best[-c(obj10$env$random)], fn = obj10$fn) # optimx(obj10$par, fn = obj10$fn, gradient=obj10$gr, method = "bobyqa", lower = lower10, upper = upper10)
SD10_bobyqa <- sdreport(obj10) # fails
time2_bobyqa_10 <- Sys.time()
time_bobyqa_10 <- time2_bobyqa_10 - time1_bobyqa_10

Report10_bobyqa = obj10$report()

data.frame(names(SD10_bobyqa$value), SD10_bobyqa$sd)
opt10_bobyqa[["final_gradient"]] = obj10$gr( opt10_bobyqa$par )
opt10_bobyqa$ierr # 
opt10_bobyqa[["AIC"]] = 2*opt10_bobyqa$fval + 2*length(opt10_bobyqa$par)

###############################



# Convergence 
opt1$final_gradient
opt2$final_gradient
opt3$final_gradient
opt4$final_gradient
opt5$final_gradient
opt6$final_gradient
opt7$final_gradient

c(opt1$message, opt2$message, opt3$message, opt4$message, opt5$message, opt6$message, opt7$message))

#--------------- AIC -------------
Model <- c("Obs", "Temporal", "Spatial", "Spatiotemporal", "Temporal + ST", "S+T+ST", "Spatial + Temporal") #
M_num <- 1:length(Model)
AIC <- c(opt1$AIC, opt2$AIC, opt3$AIC, opt4$AIC, opt5$AIC, opt6$AIC, opt7$AIC) # 
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
  
coef_table <- data.frame(Parameter = names(SD7$value), Estimate = SD7$value, BOBYQA = SD7_bobyqa$value, SD = SD7$sd, LCI, UCI, stringsAsFactors = FALSE)
for(i in 1:ncol(as.matrix(X_ij))) {
  coef_table$Parameter[i] <- colnames(as.matrix(X_ij))[i]
}
format(coef_table, digits = 2, scientific = 5)

SD_table <- data.frame(Parameter = names(SD3$value), 
                       SD1 = SD1$sd, 
                       SD2 = SD2$sd, 
                       SD3 = SD3$sd, 
                       SD4 = SD4_bobyqa$sd, 
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


