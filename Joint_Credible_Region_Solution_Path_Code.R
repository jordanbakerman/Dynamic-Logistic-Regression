
##########################################
# Joint Credible Region Solution Path Code
##########################################

mypath = function(x=x,y=y,beta_mn=NULL,beta_cov=NULL,beta_sd=NULL,max.steps=NULL,joint=TRUE){
  
  # x         = matrix of the time t vector of covariates
  # y         = vector of the binary response
  # beta_mn   = posterior mean of the state vector
  # beta_cov  = posterior covariance matrix of the state vector
  # beta_sd   = posterior standard deviation each state vector parameter
  # max.steps = maximum steps for the lars algorithm
  # joint     = if true creates joint credible regions, if false creates marginal credible regions

  p<-ncol(x)
  n<-nrow(x)
  
  if (is.null(max.steps)) {max.steps = 8 * min(n, p)}
  
  if (joint)
  {
    prec_est<-solve(beta_cov)
    
    x1 <- chol(prec_est, pivot=TRUE)
    pivot <- attr(x1, "pivot")
    x1 <- x1[, order(pivot)]
    
    y1 <- x1%*%beta_mn
    w <- beta_mn^2 
    xs<-scale(x1,center=FALSE,scale=1/w)
    adapt.cred.set<-lars(xs,y1,type="lasso",normalize=FALSE,intercept=FALSE,max.steps=max.steps,use.Gram=F)
    
    adapt.ind<-1*(coef(adapt.cred.set)!=0)
    test.df<-apply(adapt.ind,1,sum)
    
    coefs.joint<-matrix(0,nrow=length(test.df),ncol=p+1)
    SSE.joint<-rep(0,length(test.df))	
    
    for (dloc in 1:length(test.df))
    {
      if (test.df[dloc]==0){fits<-lm(y~1)}
      else
      {
        X_curr<-x[,adapt.ind[dloc,]!=0]
        fits<-lm(y~X_curr)
      }
      coefs.joint.temp<-rep(0,p)
      coefs.joint.temp[adapt.ind[dloc,]!=0]<-fits$coef[-1]
      coefs.joint[dloc,]<-c(fits$coef[1],coefs.joint.temp)
      SSE.joint[dloc]<-t(fits$res)%*%fits$res
    }
    coefs.joint=coefs.joint[test.df<=min((n-1),p),]
    SSE.joint=SSE.joint[test.df<=min((n-1),p)]
  }
  else
  {
    coefs.joint=matrix(c(lm(y~1)$coef,rep(0,p)),nrow=1)
    SSE.joint=var(y)*(n-1)
  }
  
  z.stats<-beta_mn/beta_sd
  ordered.z.s<-rank(-abs(z.stats),ties.method="random")
  
  coefs.marg<-matrix(0,nrow=min((n-1),p),ncol=p+1)
  SSE.marg<-rep(0,min((n-1),p))
  
  for (dloc in 1:min((n-1),p))
  {
    X_curr<-x[,(ordered.z.s<=dloc)]
    fits<-lm(y~X_curr)
    coefs.marg.temp<-rep(0,p)
    coefs.marg.temp[(ordered.z.s<=dloc)]<-fits$coef[-1]
    coefs.marg[dloc,]<-c(fits$coef[1],coefs.marg.temp)
    SSE.marg[dloc]<-t(fits$res)%*%fits$res
  }
  coefs.marg<-rbind(coefs.joint[1,],coefs.marg)
  SSE.marg<-c(SSE.joint[1],SSE.marg)
  
  if (!joint) 
  {
    coefs.joint=NULL
    SSE.joint=NULL		
    order.joint=NULL
    df.joint=NULL
  }
  else
  {
    order.joint=unlist(adapt.cred.set$action)
    df.joint=(1+test.df)
  }
  
  out <- list(SSE.joint=SSE.joint,SSE.marg=SSE.marg,coefs.joint=coefs.joint,
              coefs.marg=coefs.marg,df.joint=df.joint,df.marg=(1:min(n,(p+1))),
              order.joint=order.joint,order.marg=sort(ordered.z.s,index.return=T)$ix)
  return(out)
}

