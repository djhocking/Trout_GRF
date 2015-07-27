

#setwd("C:/Users/James.Thorson/Desktop/UW Hideaway/Collaborations/2015 -- river network GMRF/exploratory data/")
#######################
# Load libraries
#######################
library(TMB)

#######################
# Load data
#######################
load("Data/Prepared_Data.RData")

family <- df

#######################
# Fit in TMB
#######################

Version = "Code/OU_GMRF_v1c"
# v1a -- Original version
# v1b -- added covariates matrix X_ij
# v1c- adds linear predictors to SD output and multinomial count process
#setwd( TmbFile )
if(FALSE){
  dyn.unload(dynlib(paste0("Code/", Version)))
  file.remove( paste0("Code/", Version,c(".o",".dll")) )
}
compile( paste0(Version,".cpp") )

# Make inputs

######### temp make covariates just intercept ##########
#X_ij <- cbind(rep(x = 1, length.out = dim(c_ip)[1])) # didn't help
##############

if(Version=="Code/OU_GMRF_v1a") Data = list( "n_i"=length(c_i), "n_b"=nrow(family), "c_i"=c_i, "d_i"=family[,'child_b']-1, "parent_b"=family[,'parent_b']-1, "child_b"=family[,'child_b']-1, "dist_b"=family[,'dist_b'])
if(Version=="Code/OU_GMRF_v1b") Data = list( "n_i"=length(c_i), "n_b"=nrow(family), "c_i"=c_i, "d_i"=family[,'child_b']-1, "X_ij"=X_ij, "parent_b"=family[,'parent_b']-1, "child_b"=family[,'child_b']-1, "dist_b"=family[,'dist_b'])
if(Version=="Code/OU_GMRF_v1c") Data = list( "n_i"=dim(c_ip)[1], "n_b"=nrow(family), "c_ip"=as.matrix(c_ip), "d_i"=family[,'child_b']-1, "X_ij"=X_ij, "parent_b"=family[,'parent_b']-1, "child_b"=family[,'child_b']-1, "dist_b"=family[,'dist_b'])

if(Version=="Code/OU_GMRF_v1a") Params = list( "log_theta"=log(1), "log_SD"=log(1), "log_mean"=log(1), "Epsiloninput_d"=rnorm(Data$n_b))
if(Version=="Code/OU_GMRF_v1b") Params = list( "log_theta"=log(1), "log_SD"=log(1), "log_mean"=log(1), "gamma_j"=rep(0,ncol(Data$X_ij)), "Epsiloninput_d"=rnorm(Data$n_b))
if(Version=="Code/OU_GMRF_v1c") Params = list( "log_theta"=log(1), "log_SD"=log(1), "log_mean"=log(1), "gamma_j"=rep(0,ncol(Data$X_ij)), "log_detectrate"=log(0.2), "log_extradetectrate_i"=log(rep(1,Data$n_i)), "Epsiloninput_d"=rnorm(Data$n_b))
Random = c( "Epsiloninput_d", "log_extradetectrate_i" )

Map = NULL

# Make object
dyn.load( dynlib(paste0(Version) ))
obj <- MakeADFun(data=Data, parameters=Params, random=Random, map=Map, hessian=FALSE, inner.control=list(maxit=1000) )

# First run
obj$fn( obj$par )
# Check for parameters that don't do anything
which( obj$gr( obj$par )==0 )

# Run model
opt = nlminb(start=obj$env$last.par.best[-c(obj$env$random)], objective=obj$fn, gradient=obj$gr, control=list(eval.max=1e4, iter.max=1e4, trace=1, rel.tol=1e-14) )
opt[["final_gradient"]] = obj$gr( opt$par )

# Get standard errors
Report = obj$report()
Sdreport = sdreport( obj )

if(!exists(file.path("Output"))) dir.create(path = "Output")
capture.output( Sdreport, file="Output/SD.txt")

# compare predicted vs. observed on original scale
c_est <- exp(Report[["Epsiloninput_d"]]+Report[["log_mean"]]+Report[["eta_i"]])
plot(c_i, c_est)
abline( a=0, b=1, lty="dotted")

