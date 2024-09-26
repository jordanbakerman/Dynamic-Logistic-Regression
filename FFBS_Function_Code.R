
#############################################
# FFBS Function - Unknown Variance Components
#############################################

FFBS = function(x=NULL,zt=NULL,W=NULL,omega=NULL){
  
  # x     = matrix of the time t vector of covariates
  # zt    = the pseudo response using the latent variable augmentation approach
  # W     = the state vector covariance matrix
  # omega = polya-gamma latent variable 

  #Initialize
  n = length(zt)
  p = ncol(x)
  m = rep(0,p)
  C = diag(p)*1
  #W = diag(p)*1
  G = diag(p)
  V = 1/omega
  pred = NULL
  m.list = vector("list",length=n)
  C.list = vector("list",length=n)
  R.list = vector("list",length=n)
  a.list = vector("list",length=n)
  
  #Forward Filter
  for(i in 1:n){
    
    #Set Ft
    Ft = as.vector(x[i,])
    
    #One step ahead predictive distribution of \theta_t|y_{1:t-1}
    a = G%*%m
    R = G%*%C%*%t(G) + W
    
    #One step ahead predictive distribution of y_t|y_{1:t-1}
    f = as.numeric(Ft%*%a)
    Q = as.numeric(t(Ft)%*%R%*%Ft) + V[i]
    
    #The filtering distibution of \theta_t|y_{1:t}
    m = a + R%*%Ft*(zt[i]-f)/Q
    C = R - R%*%Ft%*%t(Ft)%*%R/Q
    
    #Collect m,C,R,a
    m.list[[i]] = m
    C.list[[i]] = C
    R.list[[i]] = R
    a.list[[i]] = a
    
  }
  
  #Backward Sampling
  h = m.list[[n]]
  H = C.list[[n]]
  theta.t1 = as.vector(t(chol(H))%*%rnorm(p) + h)
  a.t1 = as.vector(a.list[[n]])
  back.samp = rbind(theta.t1,NULL)
  
  for(i in (n-1):1){
    
    m = m.list[[i]]
    C = C.list[[i]]
    R = chol2inv(chol(R.list[[i+1]]))
    a.t1 = as.vector(a.list[[i+1]])
    
    h = m + C%*%G%*%R%*%(theta.t1-a.t1)
    H = C - C%*%G%*%R%*%G%*%C
    theta.t1 = as.vector(t(chol(H))%*%rnorm(p) + h)
    back.samp = rbind(theta.t1,back.samp)
    
  }
  
  return(back.samp)
  
}

