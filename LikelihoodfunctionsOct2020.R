
#####Main function: fits Heterogeneous Logit Model with data in long format, call by name

###########Explanatory variables: distinguishes between  
##########################        global variables (not category-specific) like gender,age
##########################        category-specific variables like price in transportation mode
##########################        if all variables are global specify any variable in nameshet, but use Indcats=0


GHMNL <- function(dat, namesglob,Indglob,namescats, Indcats,nameshet,nameresp, k,penalty){
  
  ############################################################################
  ###  Input
  
  ### dat:         data in long format (k rows for each observation ) 
  ### namesglob:   names of variables with global effects (subject-specific)
  ### Indglob >0:  global are included, otherwise ignored
  ### namescats:   names of category-specific variables (global parameters) 
  #                if no category-specific variables are in the data set, 
  #                choose any variable and set Indcats =0
  ### Indcats >0:  category-specific variables are included, otherwise ignored
  ### nameshet:    names of variables in heterogeneity term
  #                if nameshet <- c()  no heterogeneity term in the model
  ### nameresp:    name  of response variable
  ### k:           number of response category
  ### penalty:     number >= 0, if penalty=0 pure maximum likelihood
  #                             if penalty>0 ridge penalized maximum likelihood 
  
  ###  Output
  
  ### $parglob:    estimates global parameters
  ### $parglobstd: standard errors global parameter estimates
  ### $parcats:    estimates category-specific parameters
  ### $parcatsstd: standard errors category-specific parameter estimates
  ### $parhet:     estimates heterogeneity parameters
  ### $parhetstd:  standard errors heterogeneity parameter estimates
  ### $Loglik:     log-likelihood of fitted model

  ####################################################################
  
  
  hetcomp <- length(nameshet)
  #if (Indglob > 0){
  datxlength <- select(dat, namesglob)
  datxlength <- as.matrix(datxlength)
  p <- dim(datxlength)[2]
  #}
  if (Indglob <= 0){p <-0}
  
  ###
  
  pdatw <- 0
  if (Indcats > 0){ 
    datw <- select(dat, namescats)
    datw <- as.matrix(datw)
    pdatw <- dim(datw)[2]}
  
  #new
  dathlength<- select(dat, nameshet)
  dathlength <- as.matrix(dathlength) 
  #hetcomp <- dim(dathlength)[2]
  
  datresplong <- select(dat, nameresp)
  datresplong <- as.matrix(datresplong)
  
  
  #######
  #Fit<-  fitHMLMCatspec(datxlength,p,datw, pdatw, datresplong ,k,dathlength,hetcomp,penalty)
  Fit<-fitHMLMCatspecNames(datxlength,p,datw, pdatw, datresplong ,k,dathlength,hetcomp,pen,
                           namesglob,namescats,nameshet)
  
  
  #newList <- list("p" = p,"parcats"=pdatw,"parhet"= hetcomp)
  
  #return(newList)
  return(Fit)
} 
#########################################################################


######################################################################### 
#########################################################################
#### Fits heterogeneous choice model, input are data matrices, calls loglikelihhod

fitHMLMCatspecNames <- function(datxlength,p,datw, pdatw, Rblow,k,dathlength,hetcomp,pen, nameglob,namecatsp,namehet){
  
  ###jwith derivative
  ### datx       (n*k) x p matrix  global variables
  ### p          number of variables (without intercept)  
  ### datw       (n*k) x pdatw matrix design category specific variables
  ### pdatw      number of category specific variables
  ### Rblow      (n*k) x 1 matrix, 0-1 dummies for response
  ### k          number of response categories  
  ###dathlength  (n*k) x hetcomp    design matrix heterogeneity
  ###hetcomp     number of columns in dathlength (number of heterogeneity variables) 
  ###pen         penalty
  
  #output  $parglob, $parglobstd, $parcats, $parcatsstd  $parhet, $parhetstd
  
  # if (pdatw > 0) {datxblown <- blowupcatspec(datxlength,datw, k)} 
  # if (pdatw <= 0) {  datxblown <- blowupcatspec(datxlength,datw,k)  }
  
  ##new
  #datxblown <- blowupcatspec2(datxlength,datw,pdatw,k)
  datxblown <- blowupcatspec3(datxlength,p,datw,pdatw,k)
  
  ##################
  betalengthnew <- (p+1)*(k-1)+pdatw+hetcomp
  betain <- rep(0.01,betalengthnew)
  
  ###new
  #numglob <- dim(datxlength)[2]
  numglob <- p
  ###
  
  Hetz <- dathlength
  n <- dim(datxlength)[1]/k
  #mtrial <- loglikMLogitHetPenCatS(betain,df,Rblow, Hetz, n,k,numglob, 0.0001)
  
  ###jetzt mit M version matrix computation
  #esthPen <- optim(betain,fn= loglikMLogitHetPenCatSM, pen=pen , datxblown=datxblown, Hetz=Hetz,
  #                 Rblow=Rblow, n=n, k=k, numglob =numglob,gr = NULL,  
  #                 method = "BFGS", control = list(), hessian = TRUE)
  
  #####new
  if (hetcomp > 0) {esthPen <- optim(betain,fn= loglikMLogitHetPenCatSM, gr =loglikMLogitHetPenCatSMgrad,
                                     pen=pen , datxblown=datxblown, Hetz=Hetz,
                                     Rblow=Rblow, n=n, k=k, numglob =numglob, 
                                     method = "BFGS", control = list(), hessian = TRUE)}
  
  if (hetcomp <= 0){esthPen <- optim(betain,fn= loglikMLogitHetPenCatSMwithouth, gr =loglikMLogitHetPenCatSMgradwithouth,
                                     pen=pen , datxblown=datxblown, Hetz=Hetz,
                                     Rblow=Rblow, n=n, k=k, numglob =numglob, 
                                     method = "BFGS", control = list(), hessian = TRUE)}
  
  invh <-inv(esthPen$hessian)
  
  std1h<- diag(invh)
  stdhPen<-sqrt(std1h)
  
  ndum <- (k-1)*(p+1) 
  esthPen$par1 <- matrix(0,1,ndum)
  for(i in 1:ndum){esthPen$par1[1,i] <- esthPen$par[i]}
  matrixpar <- matrix(esthPen$par1,k-1,p+1,byrow=TRUE)
  
  if (pdatw > 0){ esthPen$par11 <- matrix(0,1,pdatw)
  for(i in 1:pdatw){esthPen$par11[1,i] <- esthPen$par[ndum+i]}}
  
  
  
  if (hetcomp > 0){
    esthPen$par2 <- matrix(0,1,hetcomp)
    for(i in 1:hetcomp){esthPen$par2[1,i] <- esthPen$par[ndum+pdatw+i]}
  }
  ##
  #matrixpar <- matrix(esthPen$par1,k-1,(p+1)+pdatw,byrow=TRUE)
  
  esthPen$par3 <- matrix(0,1,ndum)
  for(i in 1:ndum){esthPen$par3[1,i] <- stdhPen[i]}
  matrixparerr <- matrix(esthPen$par3,k-1,(p+1),byrow=TRUE)
  
  
  if (pdatw > 0) {esthPen$par31 <- matrix(0,1,pdatw)
  for(i in 1:pdatw){esthPen$par31[1,i] <- stdhPen[ndum+i]}}
  
  
  
  if (hetcomp > 0){esthPen$par4 <- matrix(0,1,hetcomp)
  for(i in 1:hetcomp){esthPen$par4[1,i] <- stdhPen[ndum+pdatw+i]}
  }
   
  
  if (p  > 0) {nameglob1 <- c("intercept",nameglob)}
  if (p  <= 0) {nameglob1 <- c("intercept")}
  colnames(matrixpar)<- nameglob1
  
  if (pdatw > 0) {colnames(esthPen$par11)<- namecatsp}
  if (hetcomp > 0){colnames(esthPen$par2)<- namehet}
  
  ####with Loglik
  #newList <- list("parglob" = matrixpar,"parcats"=esthPen$par11, "parhet"= esthPen$par2,"parglobstd" = matrixparerr,
  #                "parcatsstd"=esthPen$par31, "parhetstd"= esthPen$par4)
  
  newList <- list("parglob" = matrixpar,"parcats"=esthPen$par11, "parhet"= esthPen$par2,"parglobstd" = matrixparerr,
                  "parcatsstd"=esthPen$par31, "parhetstd"= esthPen$par4,"Loglik"=esthPen$value)
  
  return(newList)
}
##############end function
#################################################################################



############################################################################### 
###############################################################################
##########################Loglikelihood 
####with category-specific variables
##in datxblown already all variables, first global, then category_specific
# numglob: number of global variables


loglikMLogitHetPenCatSM <- function(betain,datxblown,Rblow, Hetz, n,k,numglob, pen){
   
  xdim<-dim(datxblown)
  pp <- dim(datxblown)[2]
  
  #lenghtx <- pbas1*(k-1)
  ####alternativ
  lenghtx <-pp- (numglob+1)
  
  lenghtz <- dim(Hetz)[2]
  
  #####
  betax <- matrix(0,1,lenghtx)
  betaz <- matrix(0,1, lenghtz)
  for(i in 1:lenghtx){betax[i]=betain[i]}
  for(i in 1:lenghtz){betaz[i] <- betain[lenghtx+i]}
  
  ###erste kategorie auf 0
  repd <- rep(0,numglob+1)
  betab <- c(repd,betax)
  betam <- matrix(betab)
  betaz <- matrix(betaz)
  
  ###neu
  
  #exp(Hetz%*%betaz)
  
  pred <-(datxblown%*%betam)*exp(Hetz%*%betaz)
  exppred <-  exp(pred)
  
  summ <-(diag(n)%x%matrix(1,1,k))%*%exppred
  sumb <- summ%x%matrix(1,k,1) 
  prob <- exppred*(1/sumb)
  logprob <- log(prob)
  
  l <- colSums(logprob*Rblow)- pen* (betain%*%betain)
  
  ####
  
  l <- -l
  return(l)
}   
##############end function
#################################################################################



############################################################################### 
###############################################################################
##########################Loglikelihood without heterogeneity
####with category-specific variables
##in datxblown already all variables, first global, then category_specific
# numglob: number of global variables

loglikMLogitHetPenCatSMwithouth <- function(betain,datxblown,Rblow, Hetz, n,k,numglob, pen){
  # loglikMLogitHetPenCatSM <- function(betain,datxblown,Rblow, Hetz, n,k,numglob, pen){
  
  xdim<-dim(datxblown)
  pp <- dim(datxblown)[2]
  
  #lenghtx <- pbas1*(k-1)
  ####alternativ
  lenghtx <-pp- (numglob+1)
  
  #new
  #lenghtz <- dim(Hetz)[2]
  
  #####
  betax <- matrix(0,1,lenghtx)
  #betaz <- matrix(0,1, lenghtz)
  for(i in 1:lenghtx){betax[i]=betain[i]}
  #for(i in 1:lenghtz){betaz[i] <- betain[lenghtx+i]}
  
  ###erste kategorie auf 0
  repd <- rep(0,numglob+1)
  betab <- c(repd,betax)
  betam <- matrix(betab)
  #betaz <- matrix(betaz)
  
  ###neu
  
  #exp(Hetz%*%betaz)
  #new
  #pred <-(datxblown%*%betam)*exp(Hetz%*%betaz)
  pred <-(datxblown%*%betam)
  
  exppred <-  exp(pred)
  
  summ <-(diag(n)%x%matrix(1,1,k))%*%exppred
  sumb <- summ%x%matrix(1,k,1) 
  prob <- exppred*(1/sumb)
  logprob <- log(prob)
  
  
  
  l <- colSums(logprob*Rblow)- pen* (betain%*%betain)
  
  ####
  
  
  l <- -l
  return(l)
}   
##############end function
#################################################################################



############################################################################### 
###############################################################################
##########################Gradient of Loglikelihood 
####with category-specific variables
##in datxblown already all variables, first global, then category_specific
# numglob: number of global variables


loglikMLogitHetPenCatSMgrad <- function(betain,datxblown,Rblow, Hetz, n,k,numglob, pen){
  
  xdim<-dim(datxblown)
  pp <- dim(datxblown)[2]
  
  #lenghtx <- pbas1*(k-1)
  ####alternativ
  lenghtx <-pp- (numglob+1)
  
  lenghtz <- dim(Hetz)[2]
  
  #####
  betax <- matrix(0,1,lenghtx)
  betaz <- matrix(0,1, lenghtz)
  for(i in 1:lenghtx){betax[i]=betain[i]}
  for(i in 1:lenghtz){betaz[i] <- betain[lenghtx+i]}
  
  ###erste kategorie auf 0
  repd <- rep(0,numglob+1)
  betab <- c(repd,betax)
  betam <- matrix(betab)
  betaz <- matrix(betaz)
  
  
  
  #exp(Hetz%*%betaz)
  
  pred <-(datxblown%*%betam)*exp(Hetz%*%betaz)
  exppred <-  exp(pred)
  
  summ <-(diag(n)%x%matrix(1,1,k))%*%exppred
  sumb <- summ%x%matrix(1,k,1) 
  prob <- exppred*(1/sumb)
  #logprob <- log(prob)
  
  
  m <- matrix(0,nrow=k-1,ncol=pp) 
  eta <- matrix(0,nrow=k,ncol=1)
  resp <- matrix(0,nrow=k-1,ncol=1)
  probv <- matrix(0,nrow=k-1,ncol=1)
  
  #i <- 2
  Grad <- matrix(0,lenghtx+lenghtz, n)
  for(i in 1:n){ hetm <- exp(Hetz[(i-1)*k +1,]%*%betaz) 
  
  for(j in 1:k-1){m[j,]<- datxblown[(i-1)*k +j+1,]
  resp[j]<- Rblow[(i-1)*k +j+1]
  probv[j] <- prob[(i-1)*k +j+1]
  }
  
  prod <- (t(m)%*%(resp-probv))*(as.numeric(hetm))
  
  nn <- numglob+2
  grad1 <- prod[nn:pp,1]
  grad11 <-  matrix(grad1)
  
  
  grad1 <-  grad11-(2*pen*matrix(betax))
  #######
  prod2 <- betax%*%grad11
  
  prod22 <- (matrix(Hetz[(i-1)*k +1,], lenghtz,1))
  
  v <- as.numeric(prod2)
  grad2 <- prod22*v-(2*pen*(betaz))
  
  grad <- c(grad1,grad2) 
  
  grad <- matrix(grad)
  
  Grad[,i] <- grad
  #
  }
  
  
  Grad <- -rowSums(Grad)
  
  return(Grad)
}   
##############end function
#################################################################################



############################################################################### 
###############################################################################
##########################Gradient of Loglikelihood without heterogeneity
####with category-specific variables
##in datxblown already all variables, first global, then category_specific
# numglob: number of global variables


loglikMLogitHetPenCatSMgradwithouth <- function(betain,datxblown,Rblow, Hetz, n,k,numglob, pen){
  
  xdim<-dim(datxblown)
  pp <- dim(datxblown)[2]
  
  #lenghtx <- pbas1*(k-1)
  ####alternativ
  lenghtx <-pp- (numglob+1)
  
  #lenghtz <- dim(Hetz)[2]
  
  #####
  betax <- matrix(0,1,lenghtx)
  #betaz <- matrix(0,1, lenghtz)
  for(i in 1:lenghtx){betax[i]=betain[i]}
  #for(i in 1:lenghtz){betaz[i] <- betain[lenghtx+i]}
  
  ###erste kategorie auf 0
  repd <- rep(0,numglob+1)
  betab <- c(repd,betax)
  betam <- matrix(betab)
  #betaz <- matrix(betaz)
  
  
  
  #exp(Hetz%*%betaz)
  
  pred <-(datxblown%*%betam)
  exppred <-  exp(pred)
  
  summ <-(diag(n)%x%matrix(1,1,k))%*%exppred
  sumb <- summ%x%matrix(1,k,1) 
  prob <- exppred*(1/sumb)
  #logprob <- log(prob)
  
  
  m <- matrix(0,nrow=k-1,ncol=pp) 
  eta <- matrix(0,nrow=k,ncol=1)
  resp <- matrix(0,nrow=k-1,ncol=1)
  probv <- matrix(0,nrow=k-1,ncol=1)
  
  #i <- 2
  Grad <- matrix(0,lenghtx, n)
  for(i in 1:n){ hetm <- 1 
  
  for(j in 1:k-1){m[j,]<- datxblown[(i-1)*k +j+1,]
  resp[j]<- Rblow[(i-1)*k +j+1]
  probv[j] <- prob[(i-1)*k +j+1]
  }
  
  prod <- (t(m)%*%(resp-probv))*(as.numeric(hetm))
  
  nn <- numglob+2
  grad1 <- prod[nn:pp,1]
  grad11 <-  matrix(grad1)
  
  
  grad1 <-  grad11-(2*pen*matrix(betax))
  #######
  #prod2 <- betax%*%grad11
  
  #prod22 <- (matrix(Hetz[(i-1)*k +1,], lenghtz,1))
  
  #v <- as.numeric(prod2)
  #grad2 <- prod22*v-(2*pen*(betaz))
  
  grad <- grad1 
  
  grad <- matrix(grad)
  
  Grad[,i] <- grad
  #
  }
  
  
  Grad <- -rowSums(Grad)
  
  return(Grad)
}   
##############end function
#################################################################################



############################################################################### 
###############################################################################
##################Split coding for ordinal predictors

fitSplitOrd <- function(datvec, cat){
  #transforms dat in k-1 dummy variables - split coding
  # works for long or not long format
  
  ##dat : data matrix nX1 variable with categories 1,...,k
  ##cat   : number categories
  
  n1 <- dim(datvec)[1]  
  
  dum <-matrix(0,n1,cat-1) 
  for(i in 1:n1){
    for(j in 2:cat){if (datvec[i,1] >= j) dum[i,j-1] <-1}
  }    
  return(dum)
}
##############end function
#################################################################################



############################################################################### 
###############################################################################
##########################Generates long format

blowupcatspec3 <- function(datxlength,pdatx,datw,pdatw, k){
  xdim<-dim(datxlength)   
  row <- xdim[1]
  
  ##changed
  #c <- xdim[2]
  c <- pdatx
  ####
  r <- row/k
  
  ##############global
  
  datx <- matrix(1,r,1)
  
  if (c > 0) {
    vec <- seq(1, row, by=k) 
    vev <- matrix(vec,r,1)
    one <- rep(1, r)
    datx <- matrix(0,r,c)
    for(i in 1:r){datx[i,]<-  datxlength[vev[i,1],]}
    
    datx<-cbind(one,datx)
  }
  
  
  d <-diag(k)
  
  
  M <- t(matrix(datx[1,]))
  kr <-kronecker(d,M)
  kr
  
  
  for(i in 2:r){M1 <- t(matrix(datx[i,]))
  kr1 <-kronecker(d,M1)
  kr  <- rbind(kr,kr1) 
  }
  
  datxblown <- kr
  datblown <-datxblown
  #######################category-specific
  if (pdatw > 0) { 
    datwdes <- matrix(0,dim(datw)[1],dim(datw)[2])
    for(i in 1:r)
    {for(j in 1:k) {datwdes[(i-1)*k+j,]<-  datw[(i-1)*k+j,]-datw[(i-1)*k+1,]
    }}
    datblown <- cbind(datxblown,datwdes)
  }
  
  return(datblown)
}
##############end function
#################################################################################



############################################################################### 
###############################################################################
##########################Gradient of Loglikelihood without heterogeneity
###########################################################possibly not needed
blowupcatspec2 <- function(datxlength,datw,pdatw, k){
  xdim<-dim(datxlength)   
  row <- xdim[1]
  c <- xdim[2]
  r <- row/k
  
  ##############global
  vec <- seq(1, row, by=k) 
  vev <- matrix(vec,r,1)
  datx <- matrix(0,r,c)
  for(i in 1:r){datx[i,]<-  datxlength[vev[i,1],]}
  
  one <- rep(1, r)
  datx<-cbind(one,datx)
  
  
  
  d <-diag(k)
  
  
  M <- t(matrix(datx[1,]))
  kr <-kronecker(d,M)
  kr
  
  
  for(i in 2:r){M1 <- t(matrix(datx[i,]))
  kr1 <-kronecker(d,M1)
  kr  <- rbind(kr,kr1) 
  }
  
  datxblown <- kr
  datblown <-datxblown
  #######################category-specific
  if (pdatw > 0) { 
    datwdes <- matrix(0,dim(datw)[1],dim(datw)[2])
    for(i in 1:r)
    {for(j in 1:k) {datwdes[(i-1)*k+j,]<-  datw[(i-1)*k+j,]-datw[(i-1)*k+1,]
    }}
    datblown <- cbind(datxblown,datwdes)
  }
  
  return(datblown)
}

##############end function
#################################################################################



############################################################################### 
###############################################################################
##########################Log-likelihood with two ridge penalty terms
####pen1 for location term parameters
####pen2 for heterogeneity term
###########################################

loglikMLogitHetPenCatSMdoublePen <- function(betain,datxblown,Rblow, Hetz, n,k,numglob, pen1,pen2){
#nur für globale Variablen!  
  
  xdim<-dim(datxblown)
  pp <- dim(datxblown)[2]
  
  #lenghtx <- pbas1*(k-1)
  ####alternativ
  lenghtx <-pp- (numglob+1)
  
  lenghtz <- dim(Hetz)[2]
  
  #####
  betax <- matrix(0,1,lenghtx)
  betaz <- matrix(0,1, lenghtz)
  for(i in 1:lenghtx){betax[i]=betain[i]}
  for(i in 1:lenghtz){betaz[i] <- betain[lenghtx+i]}
  
  ###erste kategorie auf 0
  repd <- rep(0,numglob+1)
  betab <- c(repd,betax)
  betam <- matrix(betab)
  betaz <- matrix(betaz)
  
  ###neu
  
  #exp(Hetz%*%betaz)
  
  pred <-(datxblown%*%betam)*exp(Hetz%*%betaz)
  exppred <-  exp(pred)
  
  summ <-(diag(n)%x%matrix(1,1,k))%*%exppred
  sumb <- summ%x%matrix(1,k,1) 
  prob <- exppred*(1/sumb)
  logprob <- log(prob)
  
  ##penalty vector
  
  idintercept <- seq(from = 1, to = k-1, by = numglob)
  #idhet <- seq((k-1)*(numglob+1)+1,lenghtx,1)
  
  betain1 <- betain
  for(i in 1:length(idintercept)){
    betain1[idintercept[i]]<-0
  } 
 
  
  
  l <- colSums(logprob*Rblow)- pen1* (betain1%*%betain1)-pen2* (t(betaz) %*% betaz)
  
  ####
  
  
  l <- -l
  return(l)
}   
##############end function
#################################################################################



############################################################################### 
###############################################################################
##########################Gradient of Log-likelihood with two ridge penalty terms
####pen1 for location term parameters
####pen2 for heterogeneity term
###########################################


loglikMLogitHetPenCatSMgraddoublePen <- function(betain,datxblown,Rblow, Hetz, n,k,numglob, pen1,pen2){
  
  xdim<-dim(datxblown)
  pp <- dim(datxblown)[2]
  
  #lenghtx <- pbas1*(k-1)
  ####alternativ
  lenghtx <-pp- (numglob+1)
  
  lenghtz <- dim(Hetz)[2]
  
  #####
  betax <- matrix(0,1,lenghtx)
  betaz <- matrix(0,1, lenghtz)
  for(i in 1:lenghtx){betax[i]=betain[i]}
  for(i in 1:lenghtz){betaz[i] <- betain[lenghtx+i]}
  
  ###erste kategorie auf 0
  repd <- rep(0,numglob+1)
  betab <- c(repd,betax)
  betam <- matrix(betab)
  betaz <- matrix(betaz)
  
  
  
  #exp(Hetz%*%betaz)
  
  pred <-(datxblown%*%betam)*exp(Hetz%*%betaz)
  exppred <-  exp(pred)
  
  summ <-(diag(n)%x%matrix(1,1,k))%*%exppred
  sumb <- summ%x%matrix(1,k,1) 
  prob <- exppred*(1/sumb)
  #logprob <- log(prob)
  
  
  m <- matrix(0,nrow=k-1,ncol=pp) 
  eta <- matrix(0,nrow=k,ncol=1)
  resp <- matrix(0,nrow=k-1,ncol=1)
  probv <- matrix(0,nrow=k-1,ncol=1)
  
  #i <- 2
  Grad <- matrix(0,lenghtx+lenghtz, n)
  for(i in 1:n){ hetm <- exp(Hetz[(i-1)*k +1,]%*%betaz) 
  
  for(j in 1:k-1){m[j,]<- datxblown[(i-1)*k +j+1,]
  resp[j]<- Rblow[(i-1)*k +j+1]
  probv[j] <- prob[(i-1)*k +j+1]
  }
  
  prod <- (t(m)%*%(resp-probv))*(as.numeric(hetm))
  
  nn <- numglob+2
  grad1 <- prod[nn:pp,1]
  grad11 <-  matrix(grad1)
  
  
  grad1 <-  grad11-(2*pen1*matrix(betax))
  #######
  prod2 <- betax%*%grad11
  
  prod22 <- (matrix(Hetz[(i-1)*k +1,], lenghtz,1))
  
  v <- as.numeric(prod2)
  grad2 <- prod22*v-(2*pen2*(betaz))
  
  grad <- c(grad1,grad2) 
  
  grad <- matrix(grad)
  
  Grad[,i] <- grad
  #
  }
  
  
  Grad <- -rowSums(Grad)
  
  return(Grad)
}   
##############end function
#################################################################################



############################################################################### 
###############################################################################
###########################################
############function blow up  

##### data without intercept blown up with intercept included
blowup <- function(datx,k){
  xdim<-dim(datx)  
  r <- xdim[1]
  c <- xdim[2]
  
  one <- rep(1, r)
  datx<-cbind(one,datx)
  
  
  
  d <-diag(k)
  
  
  M <- t(matrix(datx[1,]))
  kr <-kronecker(d,M)
  kr
  
  ####hier war n statt r
  for(i in 2:r){M1 <- t(matrix(datx[i,]))
  kr1 <-kronecker(d,M1)
  kr  <- rbind(kr,kr1) 
  }
  
  datxblown <- kr
  return(datxblown)
}
##############end function
#################################################################################



############################################################################### 
###############################################################################
###########################################
############function blow up  

############# blowup for heterogeneity component

blowuphet <- function(dathet,k){
  xdim<-dim(dathet)  
  r <- xdim[1]
  c <- xdim[2]
  
  
  d <-rep(1,k)
  d <- (matrix(d))
  
  
  M <- t(matrix(dathet[1,]))
  kr <-kronecker(d,M)
  kr
  
  
  for(i in 2:r){M1 <- t(matrix(dathet[i,]))
  kr1 <-kronecker(d,M1)
  kr  <- rbind(kr,kr1) 
  }
  
  datxblown <- kr
  return(datxblown)
}
##############end function
#################################################################################



############################################################################### 
###############################################################################
###########################################
############function blow up  response
################respblow


blowresp <- function(respmult,k){
  #xdim<-dim(dathet)  
  xdim<-dim(respmult)
  r <- xdim[1]
  c <- xdim[2]
  
  
  
  R1 <- matrix(0,k,1)
  res <- respmult[1]
  R1[res] <-1
  
  
  for(i in 2:r){
    Radd <- matrix(0,k,1)
    res <- respmult[i]
    Radd[res] <-1 
    R1 <- rbind(R1,Radd)
  }
  
  return(R1)
}

###################################

