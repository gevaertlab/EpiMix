#' The "Cov_Adjustment" function
#'
#' Internal. Function to necessary Fitted Observation Values (among other required variables) for Covariant-Adjusted Data.
#' @param dep CpG site for which the ideal # of Components is to be calculated.
#' @param fixed Variable within the data set to be considered the "fixed variable" in the application of Covariant-adjustment.
#' @param data matrix with methylation data for all samples (genes in columns, samples in rows).
#' @param k Number of Components to run the BIC calculation for.
#' @param weight weights of the data. Vector or colname of data. Default is NULL
#' @param family underlying type density function as a character ("gaussian" or "poisson").
#' @param pop.at.risk population at risk: An offset that could be used to determine a mixture model for Poisson
#'  data from unequally large populations at risk. Vector or colname of data. Default is NULL.
#' @param var.lnOR variances of the data: These variances might be given when working with meta analyses. 
#' Vector or colname of data. Default is NULL.
#' @param maxiter parameter to control the maximal number of iterations in the VEM and EM loops. Default is 5000.
#' @param acc convergence criterion. VEM and EM loops stop when deltaLL<acc. Default is 10^(-7).
#' @param returnHomogenousModel boolean to indicate whether the homogeneous model (simple glm) should be returned too. Default is FALSE..
#' @param family Type of distribution to perform. Default is "Gaussian".
#' @return "CAMAN_glm_object" containing calculated data.
#' @keywords internal
#'
Cov_Adjustment <- function(dep,fixed,random="",data,k,weight=NULL, pop.at.risk=NULL, 
                           var.lnOR=NULL, family="gaussian", maxiter=50, 
                           acc=10^-7, returnHomogeneousModel = FALSE){
  
  ## Initialize variables for calculation
  cl <- match.call()
  nn <- nrow(data)	
  pop.at.risk <- rep(1, nn)
  var.lnOR <- rep(1, nn)
  
  
  
  # Usual homogenous model 
  
  #y0 <- as.vector(model.extract(model.frame(formula(m), data = data), response))
  
  form<-as.formula(paste(paste(dep,"~"),paste(fixed,collapse="+"))) #dependencies for linear model
  m<-glm(form,family=family,weights=weight,data=data,x=T,na.action=na.omit, offset=log(pop.at.risk)) #compute linear model
  
  ixFixed <- which(names(data) %in% fixed)
  ixRandom <- which(names(data) %in% random)
  ixDep <- which(names(data) == dep)
  
  y <- data[,ixDep]
  y0 <- y
  ixColSort <- c(ixDep, ixFixed, ixRandom )
  
  
  idxControl_dep <- 1; names(idxControl_dep) = dep;
  idxControl_fixed <- integer(0)
  idxControl_Random <- integer(0)
  
  if(length(ixFixed)>0) {idxControl_fixed <- 2:(2+length(ixFixed)-1); names(idxControl_fixed) = names(data)[ixFixed];}	
  if(length(ixRandom)>0) {idxControl_Random <- 2:(2+length(ixRandom)-1); names(idxControl_Random) = names(data)[ixRandom];}
  idxControl <- list(ixDep = idxControl_dep, ixFixed=idxControl_fixed, ixRandom=idxControl_Random)
  colnmes = names(data)[ixColSort]
  
  data <- as.data.frame(data[,ixColSort])  #rearrange data so that the dependent variable is in the first column
  #idxControl has the indices of the rearranged data!!
  names(data) = colnmes
  mix_data <- data.frame(data[,1], rep(1,nn), pop.at.risk, var.lnOR)
  #no weights
  
  if(is.null(weight) && family=="gaussian") {
    dfg1 <- m$df.residual
    wt0<- deviance(m)/dfg1
  }
  else wt0<-1./weight
  
  pre<-fitted(m)
  
  #compute Log-Likelihood of the homogeneous model 
  logl0=0
  if (family=="gaussian")
    logl0<-sum(log(dnorm(y0,pre,sqrt(wt0))))
  
  
  LL_homo <- logl0
  
  ## If Distribution has only one Component, return the following
  if(k==1) {
    return(list(m, LL=LL_homo)) #only one component
  }
  
  
  form <- paste(paste(dep,"~",sep=""),paste("Z+",collapse="+"))
  form <- paste(form,paste(fixed,collapse="+")) 
  ### form <- paste(fixed[fixed!="1"], collapse="+")
  if((random[1] !="") || (length(random)>1)){
    random <- random[random!="1"] #cut out the intercept 
    form_rand<-paste("Z/",random,sep="",collapse="+")
    form<-paste(form,form_rand,sep="+")
  }
  form<-as.formula(paste(form,"-1"))
  
  # starting values
  p <- rep(1/k,k) #components are equal distributed 
  
  
  ixCols_fixedEffects = which(colnames(data) %in% fixed);
  ixCols_randomEffects = which(colnames(data) %in% random);
  
  if ((length(fixed)==1)&&(fixed =="1")) numPara_fixed<-0  #no fixed effects
  else if (sum(sapply(data[,ixCols_fixedEffects], is.factor))== 0) numPara_fixed<-length(fixed) ##?? numPara_fixed<-length(fixed)-1 or sum(fixed!="1") 
  else{
    if (length(ixCols_fixedEffects)==1) ixFactorsFixed = which(lapply(data, is.factor)[[ixCols_fixedEffects]])
    else ixFactorsFixed = which(unlist(lapply(data, is.factor)[ixCols_fixedEffects]))
    #if clause is just needed to seperate between [[single_idx] and [several_idx]
    number_of_levels_Fixed = unlist(lapply(data, function(xx){length(levels(xx))}))[ixCols_fixedEffects[ixFactorsFixed]]  
    numPara_fixed <- (length(fixed)-length(ixFactorsFixed)) + sum(number_of_levels_Fixed -1)
  }
  
  if((length(random) == 1) && (random == "")) numPara_random <- 0 #no random effects
  
  numPara<-numPara_fixed+k+numPara_random  #total no. of parameters
  
  b<-rep(0,numPara) #parameter initialization
  b <- seq(min(y0), max(y0), length.out=k)
  ## b[1:k] -> estimated intercepts of the components (starting values)
  
  if (numPara>k) b[(k+1):numPara] = 0
  ## b[k+1:numPara] -> starting values for other parameters
  
  # obtain solution of EM-algorithm
  mem<-mix.perform_glm(form,data,k,p,y,b,var=wt0,family=family, 
                       maxiter=maxiter, acc=acc, expected = pop.at.risk)
  logem<-mem$logl
  p <- mem$p
  #wp <- p[iii]
  x <- mem$x
  xf <- mem$xf
  m1 <- mem$m1 
  steps <- mem$n_iter
  coef<-coef(mem$m1)
  
  
  residVar <- as.numeric(NA)
  if (family=="gaussian"){ 
    residVar <- mem$residVar
  }
  
  
  
  pPosteriori<-matrix(mem$pPosteriori,nrow=nn,ncol=k)
  
  resultObj <- new("CAMAN.glm.object", dat=mix_data, family=family, 
                   num.k=k, p=p, numPara=numPara, depVar=dep, 
                   form=form, glmModel=m1, residVar = residVar)
  
  posterior_matrix <- mix.densDistr_Epi(resultObj)
  resultObj@classification = as.numeric(apply(posterior_matrix, 1, which.max)) #as.numeric(apply(posterior_matrix, 1, which.max))
  # resultObj@prob <- posterior_matrix# posterior_matrix
  resultObj@fittedObs <- fitted(m1)[resultObj@classification + ((0:(nn-1))*k)]
  #	cat("fittedObs", resultObj@fittedObs)
  
  if (returnHomogeneousModel){
    return(list(mixModel = resultObj, homoModel = list(lm=m, LL=LL_homo, BIC=BIC_homo) ) )
  }
  return(resultObj)
}

#' The "mix.densDistr_Epi" function
#'
#' Internal. Function to take data calculated in the Cov_Adjustment function and return a matrix.
#' @param mix "CAMAN_glm_obj" created within the Cov_Adjustment function and used to calculated the fitted observations
#' @return matrix containing calculated data.
#' @keywords internal
#' 
mix.densDistr_Epi <- function(mix){
  # computes the probability for each observation (1..n - row of mix@dat) belonging to each component (1..k)
  # returns a matrix of dimension n x k
  dat <- mix@dat[,1]	
  res <- matrix(ncol=mix@num.k, nrow=length(dat))
  p <- mix@p
  obs.hat <- matrix(as.vector(fitted(mix@glmModel)), ncol=mix@num.k, byrow=TRUE)
  for (i in 1:mix@num.k){
    if (mix@family == "gaussian") {
      mu <- mix@coefMatrix[1:mix@num.k,1]
      mix.sd <- sqrt(mix@residVar)
      for (j in 1:mix@num.k) res[,j] <- apply(cbind(dat,obs.hat), 1, 
                                              function(x){p[j]*dnorm(x[1], x[j+1], mix.sd ) / sum(p*dnorm(x[1], x[-1],mix.sd ))})				
    }
  }
  return(res)
}