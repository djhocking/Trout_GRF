

#setwd("C:/Users/James.Thorson/Desktop/UW Hideaway/Collaborations/2015 -- river network GMRF/exploratory data/")
#######################
# Load libraries
#######################
library(TMB)
library(dplyr)

#######################
# Load data
#######################
load("Data/Prepared_Data.RData")

# df = dataframe with all data including sites with multiple passes (multiple instances of each child)
# family = dataframe with unique child rows. Other columns are parents of each child, lat, lon, and other data associated with each child node.
# C_ip matrix of counts at site-year i on electrofish survey pass p
# X_ij matrix of covariates (j) for each site-year i as a design matrix including the intercept (column of 1s)
# t_i vector of length i indicating the survey year for each site-year visit
# df_stds dataframe with three columns of the parameter name, means, stds used for z-score standardization of the continuously-distributed independent variables in X_ij

#######################
# Fit in TMB
#######################

Version = "OU_GMRF_v1d"
# v1a -- Original version
# v1b -- added covariates matrix X_ij
# v1c- adds linear predictors to SD output and multinomial count process (HAS A BUG!)
# v1d- adds makes random variation in detection probability random, and fixed bug in detectprob calculation
#setwd( TmbFile )

# Other settings
ExtraDetectionSD = TRUE

# Compile
if(FALSE) {
  dyn.unload(dynlib(paste0("Code/", Version)))
  file.remove( paste0("Code/", Version,c(".o",".dll")) )
}
compile( paste0("Code/", Version,".cpp") )

# Make inputs

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

if(Version=="OU_GMRF_v1a") Params = list( "log_theta"=log(1), "log_SD"=log(1), "log_mean"=log(1), "Epsiloninput_d"=rnorm(Data$n_b))
if(Version=="OU_GMRF_v1b") Params = list( "log_theta"=log(1), "log_SD"=log(1), "log_mean"=log(1), "gamma_j"=rep(0,ncol(Data$X_ij)), "Epsiloninput_d"=rnorm(Data$n_b))
if(Version=="OU_GMRF_v1c") Params = list( "log_theta"=log(1), "log_SD"=log(1), "log_mean"=log(1), "gamma_j"=rep(0,ncol(Data$X_ij)), "log_detectrate"=log(0.2), "log_extradetectrate_i"=log(rep(1,Data$n_i)), "Epsiloninput_d"=rnorm(Data$n_b))
if(Version=="OU_GMRF_v1d") Params = list( "log_theta"=log(1), "log_SD"=log(1), "log_mean"=log(1), "log_extradetectionSD"=log(1), "gamma_j"=rep(0,ncol(Data$X_ij)), "log_detectrate"=log(0.2), "log_extradetectrate_i"=log(rep(1,Data$n_i)), "Epsiloninput_d"=rnorm(Data$n_b))

if(Version%in%c("OU_GMRF_v1a","OU_GMRF_v1b")) Random = c( "Epsiloninput_d" )
if(Version%in%c("OU_GMRF_v1c","OU_GMRF_v1d")) Random = c( "Epsiloninput_d", "log_extradetectrate_i" )

Map = list()
if(Version=="OU_GMRF_v1c"){
  Map[["log_extradetectrate_i"]] = factor( rep(NA,Data$n_i) )
}
if( Version%in%c("OU_GMRF_v1d") & ExtraDetectionSD==FALSE ){
  Map[["ExtraDetectionSD"]] = factor(NA)
}

# Make object
dyn.load( dynlib(paste0("Code/", Version) ))
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

df_N <- data.frame(child_name = df$child_name, c_sum = rowSums(df[ , c("pass_1", "pass_2", "pass_3")]), N_i, N_unbias, N_sd, p = Report$detectprob_ip, pass_1 = df$pass_1, pass_2 = df$pass_2, pass_3 = df$pass_3)

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
  dplyr::select(c_sum, N_i, N_unbias, N_sd, pass_1, pass_2, pass_3, p_miss_total) %>%
  dplyr::filter(!is.na(c_sum))

# save output for Kyle to map
df_kyle <- df_N %>%
  dplyr::select(child_name, N_100, rho_b) # if rho doesn't change filter to unique child

write.csv(df_kyle, file = "Output/N_correlation.csv")


