

#setwd("C:/Users/James.Thorson/Desktop/UW Hideaway/Collaborations/2015 -- river network GMRF/exploratory data/")
#setwd("C:/Users/James.Thorson/Desktop/Project_git/Trout_GRF/")
#######################
# Load libraries
#######################
library(TMB)
library(dplyr)

#######################
# Load data
#######################
load("Data/Prepared_Data.RData")

# remove year from X_ij now so it doesn't mess with testing of temporal and temporal-spatial mdoels
X_ij <- X_ij[ , c("(Intercept)", "length_std", "width_std")]

# df = dataframe with all data including sites with multiple passes (multiple instances of each child)
# family = dataframe with unique child rows. Other columns are parents of each child, lat, lon, and other data associated with each child node.
# C_ip matrix of counts at site-year i on electrofish survey pass p
# X_ij matrix of covariates (j) for each site-year i as a design matrix including the intercept (column of 1s)
# t_i vector of length i indicating the survey year for each site-year visit
# df_stds dataframe with three columns of the parameter name, means, stds used for z-score standardization of the continuously-distributed independent variables in X_ij

#######################
# Fit in TMB
#######################

Version = "OU_GMRF_v1f"
# v1a -- Original version
# v1b -- added covariates matrix X_ij
# v1c- adds linear predictors to SD output and multinomial count process (HAS A BUG!)
# v1d- adds makes random variation in detection probability random, and fixed bug in detectprob calculation
# v1e- adds Temporal variation as AR1 process
# v1f- adds Spatiotemporal variation as kronecker product of O-U and AR1 process, plus interface to turn off components
#setwd( TmbFile )

# Turn off random effects in v1f (0 means exclude a component, except for ObsModel)
Options_vec = c("SpatialTF"=1, "TemporalTF"=1, "SpatiotemporalTF"=1, "DetectabilityTF"=1, "ObsModel"=1)

# YearSet
YearSet = min(t_i):max(t_i)

# Compile
if(FALSE) {
  dyn.unload(dynlib(paste0("Code/", Version)))
  file.remove( paste0("Code/", Version,c(".o",".dll")) )
}
compile( paste0("Code/", Version,".cpp") )

# Make inputs
rmatrix = function( nrow=1, ncol=1, mean=0, sd=1, ... ){
  Return = matrix( rnorm(nrow*ncol,mean=mean,sd=sd), nrow=nrow, ncol=ncol, ...)
  return( Return )
}

# convert 3-pass counts to abundance using Carle & Strub 1978 methods
if(Version=="OU_GMRF_v1b") {
  #source("http://www.rforge.net/FSA/InstallFSA.R")
  if(!"FSA" %in% installed.packages()) {
    if (!require('devtools')) install.packages('devtools')
    require('devtools')
    devtools::install_github('droglenc/FSA')
  }
  library(FSA)
  c_obs <- dplyr::filter(c_ip, !is.na(pass_1))
  N_cs <- as.data.frame(matrix(NA, nrow(c_obs), 8))
  for(i in 1:dim(c_obs)[1]){
    foo <- FSA::removal(catch = as.vector(c_obs[i, ]), method = c("CarleStrub"), just.ests = TRUE)
    N_cs[i, ] <- foo
  } # end carle-strub for loop
  names(N_cs) <- names(foo)
  data.frame(c_obs, c_sum = rowSums(c_obs), N_cs)
  foo <- c_ip %>%
    dplyr::select(pass_1) %>%
    dplyr::filter(is.na(pass_1))
  c_i <- as.vector(dplyr::bind_rows(foo, dplyr::select(N_cs, No))$No)
} # end if version b statement

if(Version=="OU_GMRF_v1a") Data = list( "n_i"=length(c_i), "n_b"=nrow(family), "c_i"=c_i, "d_i"=family[,'child_b']-1, "parent_b"=family[,'parent_b']-1, "child_b"=child_b-1, "dist_b"=family[,'dist_b'])
if(Version=="OU_GMRF_v1b") Data = list( "n_i"=length(c_i), "n_b"=nrow(family), "c_i"=c_i, "d_i"=family[,'child_b']-1, "X_ij"=X_ij, "parent_b"=family[,'parent_b']-1, "child_b"=child_b-1, "dist_b"=family[,'dist_b'])
if(Version%in%c("OU_GMRF_v1c","OU_GMRF_v1d")) Data = list( "n_i"=dim(c_ip)[1], "n_b"=nrow(family), "c_ip"=as.matrix(c_ip), "d_i"=df[,'child_b']-1, "X_ij"=X_ij, "parent_b"=family[ ,'parent_b']-1, "child_b"=family[ ,'child_b']-1, "dist_b"=family[,'dist_b'])
if(Version%in%c("OU_GMRF_v1e")) Data = list( "n_i"=dim(c_ip)[1], "n_b"=nrow(family), "n_t"=length(YearSet), "c_ip"=as.matrix(c_ip), "d_i"=df[,'child_b']-1, "X_ij"=X_ij, "t_i"=t_i-min(t_i), "parent_b"=family[ ,'parent_b']-1, "child_b"=family[ ,'child_b']-1, "dist_b"=family[,'dist_b'])
if(Version%in%c("OU_GMRF_v1f")) Data = list( "Options_vec"=Options_vec, "n_i"=dim(c_ip)[1], "n_b"=nrow(family), "n_t"=length(YearSet), "c_ip"=as.matrix(c_ip), "d_i"=df[,'child_b']-1, "X_ij"=X_ij, "t_i"=t_i-min(t_i), "parent_b"=family[ ,'parent_b']-1, "child_b"=family[ ,'child_b']-1, "dist_b"=family[,'dist_b'])

if(Version=="OU_GMRF_v1a") Params = list( "log_theta"=log(1), "log_SD"=log(1), "log_mean"=log(1), "Epsiloninput_d"=rnorm(Data$n_b))
if(Version=="OU_GMRF_v1b") Params = list( "log_theta"=log(1), "log_SD"=log(1), "log_mean"=log(1), "gamma_j"=rep(0,ncol(Data$X_ij)), "Epsiloninput_d"=rnorm(Data$n_b))
if(Version=="OU_GMRF_v1c") Params = list( "log_theta"=log(1), "log_SD"=log(1), "log_mean"=log(1), "gamma_j"=rep(0,ncol(Data$X_ij)), "log_detectrate"=log(0.2), "log_extradetectrate_i"=log(rep(1,Data$n_i)), "Epsiloninput_d"=rnorm(Data$n_b))
if(Version=="OU_GMRF_v1d") Params = list( "log_theta"=log(1), "log_SD"=log(1), "log_mean"=log(1), "log_extradetectionSD"=log(1), "gamma_j"=rep(0,ncol(Data$X_ij)), "log_detectrate"=log(0.2), "log_extradetectrate_i"=log(rep(1,Data$n_i)), "Epsiloninput_d"=rnorm(Data$n_b))
if(Version=="OU_GMRF_v1e") Params = list( "log_theta"=log(1), "log_SD"=log(1), "log_mean"=log(1), "log_extradetectionSD"=log(1), "rhot"=0, "log_sigmat"=log(1), "gamma_j"=rep(0,ncol(Data$X_ij)), "log_detectrate"=log(0.2), "log_extradetectrate_i"=log(rep(1,Data$n_i)), "Epsiloninput_d"=rnorm(Data$n_b), "Deltainput_t"=rnorm(Data$n_t))
if(Version=="OU_GMRF_v1f") Params = list( "log_theta"=log(1), "log_SD"=log(1), "log_theta_sp"=log(1), "log_SD_sp"=log(1), "rho_sp"=0, "log_mean"=log(1), "log_extradetectionSD"=log(1), "rhot"=0, "log_sigmat"=log(1), "gamma_j"=rep(0,ncol(Data$X_ij)), "log_detectrate"=log(0.2), "log_extradetectrate_i"=rnorm(Data$n_i,sd=0.01), "Epsiloninput_d"=rnorm(Data$n_b,sd=0.01), "Deltainput_t"=rnorm(Data$n_t,sd=0.01), "Nu_dt"=rmatrix(Data$n_b,Data$n_t,sd=0.01))

if(Version%in%c("OU_GMRF_v1a","OU_GMRF_v1b")) Random = c( "Epsiloninput_d" )
if(Version%in%c("OU_GMRF_v1c","OU_GMRF_v1d")) Random = c( "Epsiloninput_d", "log_extradetectrate_i" )
if(Version%in%c("OU_GMRF_v1e")) Random = c( "Epsiloninput_d", "log_extradetectrate_i", "Deltainput_t" )
if(Version%in%c("OU_GMRF_v1f")) Random = c( "Epsiloninput_d", "log_extradetectrate_i", "Deltainput_t", "Nu_dt" )

# Turn off random effects if desired
Map = list()
if( Version%in%c("OU_GMRF_v1f","OU_GMRF_v1e","OU_GMRF_v1d") & Options_vec[["SpatialTF"]]==FALSE ){
  Map[["log_theta"]] = factor(NA)
  Map[["log_SD"]] = factor(NA)
  Params[["Epsiloninput_d"]] = rep(0,length("Epsiloninput_d"))
  Map[["Epsiloninput_d"]] = factor( rep(NA,length("Epsiloninput_d")) )
}
if( Version%in%c("OU_GMRF_v1f","OU_GMRF_v1e","OU_GMRF_v1d") & Options_vec[["TemporalTF"]]==FALSE ){
  Map[["rhot"]] = factor(NA)
  Map[["log_sigmat"]] = factor(NA)
  Params[["Deltainput_t"]] = rep(0,length(Params[["Deltainput_t"]]))
  Map[["Deltainput_t"]] = factor( rep(NA,length(Params[["Deltainput_t"]])) )
}
if( Version%in%c("OU_GMRF_v1f","OU_GMRF_v1e","OU_GMRF_v1d") & Options_vec[["SpatiotemporalTF"]]==FALSE ){
  Map[["log_theta_sp"]] = factor(NA)
  Map[["log_SD_sp"]] = factor(NA)
  Map[["rho_sp"]] = factor(NA)
  Params[["Nu_dt"]] = array(0,dim(Params[["Nu_dt"]]))
  Map[["Nu_dt"]] = factor( array(NA,dim(Params[["Nu_dt"]])) )
}
if( Version%in%c("OU_GMRF_v1f","OU_GMRF_v1e","OU_GMRF_v1d") & Options_vec[["DetectabilityTF"]]==FALSE ){
  Map[["ExtraDetectionSD"]] = factor(NA)
  Params[["log_extradetectrate_i"]] = rep(0,Data$n_i)
}

# Make object
dyn.load( dynlib(paste0("Code/", Version) ))
obj <- MakeADFun(data=Data, parameters=Params, random=Random, map=Map, hessian=FALSE, inner.control=list(maxit=1000) )

# First run
obj$fn( obj$par )
# Check for parameters that don't do anything
Which = which( obj$gr( obj$par )==0 )

# Run model
opt = nlminb(start=obj$env$last.par.best[-c(obj$env$random)], objective=obj$fn, gradient=obj$gr, control=list(eval.max=1e4, iter.max=1e4, trace=1, rel.tol=1e-14) )
opt[["final_gradient"]] = obj$gr( opt$par )

# Get standard errors
Report = obj$report()
Sdreport = sdreport( obj )
SD = sdreport( obj, bias.correct=TRUE )
SD$unbiased$value

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


