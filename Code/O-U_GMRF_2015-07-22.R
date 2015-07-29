

# let the .rproj file set the working directory then do subdirectories relative to that for code portability

# setwd("C:/Users/James.Thorson/Desktop/UW Hideaway/Collaborations/2015 -- river network GMRF/exploratory triple-pass depletion/")

#######################
# Load data
#######################

load( "Data/sample.RData")
colnames(family)[1] = "child_name"
family = cbind( family, "child_b"=1:nrow(family) )

#######################
# Simulate GMRF following O-U process
#######################

# parameters
theta = 1
SD = 1
log_mean = 2

# object
condSD_b = x_b = rep(NA, nrow(family))

# seed at top of network
WhichRoot = which( is.na(family[,'parent_b']) )
condSD_b[WhichRoot] = sqrt( SD^2 / 2*theta )
x_b[WhichRoot] = rnorm(1, mean=0, sd=condSD_b[WhichRoot])

# Loop through network
while( TRUE ){
  for(i in 1:nrow(family)){
    if( is.na(x_b[i]) ){
      SimulatedNodes = which(!is.na(x_b))
      Match = match( family[i,'parent_b'], SimulatedNodes ) # Which
      if(length(Match)==1){
        condSD_b[i] = sqrt( SD^2/(2*theta) * (1-exp(-2*theta*family[i,'dist_b'])) )
        x_b[i] = x_b[SimulatedNodes[Match]] + rnorm(1, mean=0, sd=condSD_b[i])
      }
    }
  }
  # Stopping condition
  if( all(!is.na(x_b)) ) break()
}

# Simulate poisson count process (for latent density at each site)
c_i = rpois( length(x_b), lambda=exp(x_b+log_mean))
X_ij = cbind( rnorm(length(c_i)) )

# Simulate multinomial count process (for count of each depletion pass)
detectprob = 0.2
passprob = (1-detectprob)^(0:3) * detectprob^(c(1,1,1,0))
c_ip = NULL
for(i in 1:length(c_i)) c_ip = rbind(c_ip, rmultinom(n=1, size=c_i[i], prob=passprob)[,1] )
c_ip = c_ip[,1:3]

#######################
# Fit in TMB
#######################

library( TMB )

Version = "Code/OU_GMRF_v1c"
# v1a -- Original version
# v1b -- added covariates matrix X_ij
# v1c -- added triple-pass depletion with unobservable variation in area-swept/catchability among sites
#setwd( TmbFile )
if(FALSE){
  dyn.unload(dynlib(Version))
  file.remove( paste0(Version,c(".o",".dll")) )
}
compile( paste0(Version,".cpp") )

# Make inputs
if(Version=="Code/OU_GMRF_v1a") Data = list( "n_i"=length(c_i), "n_b"=nrow(family), "c_i"=c_i, "d_i"=family[,'child_b']-1, "parent_b"=family[,'parent_b']-1, "child_b"=family[,'child_b']-1, "dist_b"=family[,'dist_b'])
if(Version=="Code/OU_GMRF_v1b") Data = list( "n_i"=length(c_i), "n_b"=nrow(family), "c_i"=c_i, "d_i"=family[,'child_b']-1, "X_ij"=X_ij, "parent_b"=family[,'parent_b']-1, "child_b"=family[,'child_b']-1, "dist_b"=family[,'dist_b'])
if(Version=="Code/OU_GMRF_v1c") Data = list( "n_i"=length(c_i), "n_b"=nrow(family), "c_ip"=c_ip, "d_i"=family[,'child_b']-1, "X_ij"=X_ij, "parent_b"=family[,'parent_b']-1, "child_b"=family[,'child_b']-1, "dist_b"=family[,'dist_b'])
if(Version=="Code/OU_GMRF_v1a") Params = list( "log_theta"=log(1), "log_SD"=log(1), "log_mean"=log(1), "Epsiloninput_d"=rnorm(Data$n_b))
if(Version=="Code/OU_GMRF_v1b") Params = list( "log_theta"=log(1), "log_SD"=log(1), "log_mean"=log(1), "gamma_j"=rep(0,ncol(Data$X_ij)), "Epsiloninput_d"=rnorm(Data$n_b))
if(Version=="Code/OU_GMRF_v1c") Params = list( "log_theta"=log(1), "log_SD"=log(1), "log_mean"=log(1), "gamma_j"=rep(0,ncol(Data$X_ij)), "log_detectrate"=log(0.2), "log_extradetectrate_i"=log(rep(1,Data$n_i)), "Epsiloninput_d"=rnorm(Data$n_b))
Random = c( "Epsiloninput_d", "log_extradetectrate_i" )

# Turn off site-level heterogeneity in catchability (a.k.a. unobserved area-swept)
Map = list()
#Map[["log_extradetectrate_i"]] = factor( c(NA,2:Data$n_i) )
Map[["log_extradetectrate_i"]] = factor( rep(NA,Data$n_i) )

# Make object
dyn.load( dynlib(Version) )
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

# Compare predicted vs. observed
plot( x=Report[["Epsiloninput_d"]]+Report[["log_mean"]], y=x_b+log_mean )
abline( a=0, b=1, lty="dotted")

# estimated catchability
detectprob_hat = 1 - exp(-exp(opt$par['log_detectrate']))
colSums(Report$detectprob_ip) / sum(Report$detectprob_ip)
colSums( c_ip ) / sum(c_ip)

Report$detectprob_ip
