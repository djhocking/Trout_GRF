runOUGMRF <- function(inputs) {
  # Make object
  dyn.load( dynlib(paste0("Code/", Version )))
  obj1 <- MakeADFun(data=inputs$Data, parameters=inputs$Params, random=inputs$Random, map=inputs$Map, hessian=FALSE, inner.control=list(maxit=1000) )
  
  # First run
  obj1$fn( obj1$par )
  
  # set initial so return and booleans don't fail
  opt1 <- NULL
  Report1 <- NULL
  SD1 <- NULL
  
  try({
    # Run model
    opt1 = nlminb(start=obj1$env$last.par.best[-c(obj1$env$random)], objective=obj1$fn, gradient=obj1$gr, control=list(eval.max=50, iter.max=50, trace=1, rel.tol=1e-14) )
    opt1[["Param"]] = names( opt1$par )
    opt1[["final_gradient"]] = obj1$gr( opt1$par )
    opt1[["AIC"]] = 2*opt1$objective + 2*length(opt1$par)
    opt1[["k"]] = length(opt1$par)
    Report1 = obj1$report()
    Report1[["Optimizer"]] <- "nlminb"
    SD1 = sdreport( obj1, bias.correct=FALSE )
    
    # Use BOBYQA optimization if PORTS fails
    #if(is.null(SD1)) {
    if(opt1$convergence == 1 | max(abs(opt1$final_gradient)) > 0.001 | is.null(SD1)) {
      opt1 <- bobyqa(par = obj1$env$last.par.best[-c(obj1$env$random)], fn = obj1$fn)
      opt1[["convergence"]] <- opt1$ierr
      Report1 = obj1$report()
      Report1[["Optimizer"]] <- "bobyqa"
      opt1[["Param"]] = names( opt1$par )
      opt1[["k"]] = length(opt1$par)
      opt1[["AIC"]] = 2*opt1$fval + 2*length(opt1$par)
      SD1 <- sdreport(obj1, bias.correct=TRUE )
    #}
     }
    }) # end try - consider trycatch and having error return NULL or NA or something
    return(list(opt = opt1, Report = Report1, SD = SD1, obj=obj1))
  #}) # end try - consider trycatch and having error return NULL or NA or something
}

