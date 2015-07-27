
#######################
#
# SUMMARY
# 1. Multinomial in R works
# 2. Poisson in R doesn't work (probably because it doesn't converge)
# 3. Poisson in TMB works
#
#######################

# Data inputs
Data_ij = matrix( c(3,61,232, 80,137,400, 216,131,301, 268,76,203, 197,50,188, 150,24,164, 91,10,183), byrow=TRUE, ncol=3, nrow=7)
colSums(Data) / sum(Data)

###############
# Multinomial regression in R
###############
Multinom_fn = function( Par, Return="NLL" ){
  Prob_j = exp( c(0,Par) )
  Prob_j = Prob_j / sum(Prob_j)
  NLL = 0
  for(i in 1:nrow(Data_ij)) NLL = NLL + -1*dmultinom( Data_ij[i,], prob=Prob_j, log=TRUE )
  if(Return=="NLL") Return = NLL
  if(Return=="Par_hat") Return = list("Prob_j"=Prob_j)
  return(Return)
}
Multinom_fn( rep(0,ncol(Data_ij)-1) )
opt1 = optim( par=rep(0,ncol(Data_ij)-1), fn=Multinom_fn, hessian=TRUE )
Multinom_fn( opt1$par, Return="Par_hat")

###############
# Poisson regression in R
###############
Pois_fn = function( Par, Return="NLL" ){
  Eta = Par[1]
  Alpha_j = c(0, Par[2:(ncol(Data_ij))])
  Prob_j = exp(Alpha_j) / sum(exp(Alpha_j))
  Beta_i = c(Par[-c(1:(ncol(Data_ij)))])
  NLL = 0
  Lambda = exp(Eta + outer(Beta_i,rep(1,ncol(Data_ij))) + outer(rep(1,nrow(Data_ij)),Alpha_j) )
  NLL = -1*sum( dpois(x=Data_ij, lambda=Lambda, log=TRUE) )
  if(Return=="NLL") Return = NLL
  if(Return=="Par_hat") Return = list("Prob_j"=Prob_j, "Eta"=Eta, "Alpha_j"=Alpha_j, "Beta_i"=Beta_i, "Lambda"=Lambda)
  return(Return)
}
Pois_fn( rep(0,ncol(Data_ij)+nrow(Data_ij)-0) )
opt2 = optim( par=rep(0,ncol(Data_ij)+nrow(Data_ij)-0), fn=Pois_fn, hessian=TRUE )
( ParHat = Pois_fn( opt2$par, Return="Par_hat") )
all( eigen(opt2$hessian)$values > (-1e-6) )
colSums( ParHat$Lambda )
rowSums( ParHat$Lambda )

###############
# Poisson in TMB
###############
library( TMB )
setwd( "C:/Users/James.Thorson/Desktop/UW Hideaway/Software/R/Thorson/Multinomial-Poisson transformation/")
compile( "poisson_v1.cpp" )

# Make object
dyn.load( dynlib("poisson_v1"))
Data = list("Data_ij"=Data_ij) #"n_i"=nrow(Data_ij), "n_j"=ncol(Data_ij))
Param = list("Alpha_j"=rep(0,ncol(Data_ij)), "Beta_i"=rep(0,nrow(Data_ij)) )
Map = list( "Beta_i"=factor(c(NA,rep(1:(length(Param[["Beta_i"]])-1)))) )
Random = c("Beta_i")
obj = MakeADFun( data=Data, parameters=Param, map=Map, random=Random)

# Optimize
opt = nlminb( start=obj$par, objective=obj$fn, gradient=obj$gr)
Report = obj$report()
Report$Prob_j
SD = sdreport( obj )
summary(SD)
# True
colSums(Data_ij) / sum(Data_ij)
