# Weighted regression of y on X with weights w
WeightedRegression <- function(X,y,w=array(1,length(y),1) ){
  beta<-solve(t(X*w)%*%X,t(X*w)%*%y)
  return(beta)
}

# Weighted mean of y with weights w
WeightedMean <- function(y,w=array(1,c(length(y),1)) ){
  mu<-sum(y%*%w)/sum(w)
  return(mu)
}

# Estimator of mixed 4th moments
EstimateK4 <-function(ell,z.1,z.2,crosstrait.intercept=1,ldsc.intercept=1,weights=1/pmax(1,ell),sig.threshold=.Machine$integer.max,n.1=1,n.2=1,intercept.12=0,nargout=3){
  
  # LDSC regression on each trait
  if(ldsc.intercept==0){
    intercept.1<-1/n.1
    intercept.2<-1/n.2
    temp<-WeightedRegression(ell,z.1^2-intercept.1,weights)
    h2g.1<-temp[1]
    temp<-WeightedRegression(ell,z.2^2-intercept.2,weights)
    h2g.2<-temp[1]
  }
  else{
    # Exclude significant SNPs when calculating LDSC intercept
    keep.snps.1<-z.1^2 <= sig.threshold*mean(z.1^2)
    temp<-WeightedRegression(cbind(ell[keep.snps.1],matrix(1,sum(keep.snps.1))),z.1[keep.snps.1]^2,weights[keep.snps.1]);
    intercept.1<-temp[2];
    
    # Include significant SNPs when computing heritability
    temp<-WeightedRegression(ell,z.1^2-intercept.1,weights);
    h2g.1<-temp[1];
    
    keep.snps.2<-z.2^2 <= sig.threshold*mean(z.2^2)
    temp<-WeightedRegression(cbind(ell[keep.snps.2],matrix(1,sum(keep.snps.2))),z.2[keep.snps.2]^2,weights[keep.snps.2]);
    intercept.2<-temp[2];
    
    temp<-WeightedRegression(ell,z.2^2-intercept.2,weights);
    h2g.2<-temp[1];  
  }  

  # Cross-trait LDSC regression
  if(crosstrait.intercept==0) {
    temp<-WeightedRegression(ell,z.1*z.2-intercept.12,weights)
    rho<-temp/sqrt(h2g.1*h2g.2)
  }
  else{
    keep.snps.12=(z.1^2<sig.threshold*mean(z.1^2))*(z.2^2<sig.threshold*mean(z.2^2))==1
    temp<-WeightedRegression(cbind(ell[keep.snps.12],matrix(1,sum(keep.snps.12))),z.1[keep.snps.12]*z.2[keep.snps.12],weights[keep.snps.12])
    intercept.12<-temp[2]
    temp<-WeightedRegression(ell,z.1*z.2-intercept.12,weights)
    rho<-temp[1]/sqrt(h2g.1*h2g.2)
  }
  
  # Normalize effect sizes 
  s.1<-sqrt(WeightedMean(z.1^2,weights)-intercept.1)
  s.2<-sqrt(WeightedMean(z.2^2,weights)-intercept.2)
  nz.1<-z.1/s.1
  nz.2<-z.2/s.2
  
  # Estimates of mixed 4th moments
  k41<-WeightedMean(nz.2*nz.1^3-3*nz.1*nz.2*(intercept.1/s.1^2)-3*(nz.1^2-intercept.1/s.1^2)*intercept.12/s.1/s.2,weights)
  k42<-WeightedMean(nz.1*nz.2^3-3*nz.1*nz.2*(intercept.2/s.2^2)-3*(nz.2^2-intercept.2/s.2^2)*intercept.12/s.1/s.2,weights)
  
  argout<-c(rho,k41,k42,intercept.12,s.1,s.2,intercept.1,intercept.2)
  return(argout[1:nargout])
  
}
