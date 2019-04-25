
##### Simulate bivariate Lomax and geometric random variables

# n is number of observations
rblomaxgeo<- function(n,alpha,beta,p){
  N<-rgeom(n,p)+1
  X<-rgamma(n,shape =N,rate = beta )
  X<-X/rgamma(1,shape =1/alpha,rate = 1/alpha )
  return(data.frame(X,N))
}

# data is bivariate vector  (X,N) vector representing observations from bivariate Lomax and geometric model
dblomaxgeo<- function(data,alpha,beta,p){
  N<-data[,2]
  X<-data[,1]
  M<-((1-p)^(N-1))*(p*beta^(alpha*N))*(X^(alpha*N-1))*exp(-beta*X)[3]
  return(M)
}

# q is bivariate vector  (X,N) vector quantiles from bivariate Lomax and geometric model
pblomaxgeo<- function(q,alpha,beta,p, lower.tail=TRUE,log.p=FALSE){
  N<-q[,2]
  X<-q[,1]
  M<-NULL
  t=1
  for (i in 1:length(N)) {
    k=seq(1,N[i])
    S0<-0
    for (j in k) {
      S0<-S0 + p*((1-p)^(j-1))*pracma::gammainc(beta*X[i],j*alpha)[3] 
    }
    M[t]<- S0
    t=t+1
  }
  if (lower.tail==FALSE & log.p==TRUE){
    M<-log(1-M)
    return(M)
  }
  else if (lower.tail==TRUE & log.p==TRUE){
    M<-log(M)
    return(M) 
  }
  else if (lower.tail==TRUE & log.p==FALSE){
    return(M) 
  }
  else {
    return(1-M)
  }
}


## Eample
Data.df<-rblomaxgeo(100,alpha = 1,beta = 100,p=0.4)
prob<-pblomaxgeo(Data.df,alpha = 1,beta = 100,p=0.4)
prob


