# RunLCV runs LCV on summary statistics for two traits.
#	Input vectors are sorted lists of LD scores and signed summary statistics. They
#		must be sorted by genomic position, as LCV uses a block-jackknife procedure
#		to compute standard errors; if consecutive SNPs are not approximately
#		contiguous, standard errors will be underestimated.
#
#   INPUT VARIBLES: ell, Mx1 vector of LD scores; 
#	z.1, Mx1 vector of estimated marginal per-normalized-genotype effects on trait 1
#   	(or Z scores; invariant to scaling);  
#	z.2, Mx2 vector of effects on trait2;  
#	crosstrait.intercept, 0 if cohorts are disjoint or overlap is known, 1 if cohorts 
#   	are possibly nondisjoint and necessary correction is unknown
#   ldsc.intercept, 0 if intercept should be fixed at 1 and 1 otherwise;
#   weights, Mx1 vector of regression weights; sig.threshold, threshold
#   	above which to discard chisq statistics for the purpose of estimating
#   	the LDSC intercept if they are above sig_threshold*mean(chisq);
#   no.blocks, number of jackknife blocks.
#   n1, sample size for trait 1, only needed if ldsc.intercept=1; 
#   n2, sample size for trait 2;  
#	intercept.12, covariance between sampling errors for Z1 and Z2, only needed if
#   	crosstrait_intercept=0.
#
#   OUTPUT VARIABLES: lcv.output, a data frame with columns:
#   "zscore", Z score for partial genetic causality; 
#   	to obtain 2-tailed p-value from zscore, compute:
#   	x<-pt(lcv.output[["zscore"]],no.blocks-2); pval<-2*min(x,1-x).
#   	Significantly positive zscore implies trait 1 partially genetically causal for trait 2.
#   "gcp.pm", posterior mean gcp (positive: trait 1 -> trait 2); 
#   "gcp.pse", posterior standard error for gcp; 
#   "rho.est", estimated genetic correlation; 
#   "rho.err", standard error of rho estimate;
#   "pval.fullycausal.1", p-value for null that gcp=-1; 
#   "pval.fullycausal.2", p-value for null that gcp=1; 
#   "h2.zscore.1", z score for trait 1 being significantly heritable;
#   	we recommend reporting results for h2.zscore.k > 7 (a very stringent threshold).
#   "h2.zscore.2", z score for trait 2 being significantly heritable.

RunLCV <- function(ell,z.1,z.2,no.blocks=100,crosstrait.intercept=1,ldsc.intercept=1,weights=1/pmax(1,ell),sig.threshold=.Machine$integer.max,n.1=1,n.2=1,intercept.12=0,nargout=5){
  grid<- (-100:100)/100;
  mm=length(z.1)
  
  source("MomentFunctions.R")
  
  # Jackknife estimates of moments
  size.blocks=floor(mm/no.blocks)
  jackknife=matrix(0,no.blocks,8)
  for(jk in 1:no.blocks){
    if(jk==1)
    {ind<-(size.blocks+1):mm}
    else if(jk==no.blocks)
    {ind <- 1:(size.blocks*(jk-1))}
    else
    {ind<-c(1:((jk-1)*size.blocks), (jk*size.blocks+1):mm)}
    jackknife[jk,] <- EstimateK4(ell[ind],z.1[ind],z.2[ind],crosstrait.intercept,ldsc.intercept,weights[ind],sig.threshold,n.1,n.2,intercept.12,8)
  }
  
  jackknife[,2:3]<-jackknife[,2:3]-3*jackknife[,1]
  
  estimate<-mean(jackknife)
  error<-sd(jackknife)*sqrt(no.blocks+1)
  
  # Likelihood of each gcp value
  gcp.likelihood=grid;gcp.likelihood[]=0
  for(kk in 1:length(grid)){
    xx<-grid[kk]
    fx<-abs(jackknife[,1])^(-xx)
    numer<-jackknife[,2]/fx-fx*jackknife[,3]
    denom=pmax(1/abs(jackknife[,1]),sqrt(jackknife[,2]^2/fx^2 + jackknife[,3]^2*fx^2 ))
    pct.diff<-numer/denom
    
    gcp.likelihood[kk]<-dt(mean(pct.diff)/sd(pct.diff)/sqrt(no.blocks+1),no.blocks-2)
    
    if(xx==-1){
      pval.fullycausal.2<-pt(mean(pct.diff)/sd(pct.diff)/sqrt(no.blocks+1),no.blocks-2)
    }
    if(xx==1){
      pval.fullycausal.1<-pt(-mean(pct.diff)/sd(pct.diff)/sqrt(no.blocks+1),no.blocks-2)
    }
    if(xx==0){
      zscore<- mean(pct.diff)/sd(pct.diff)/sqrt(no.blocks+1)
    }
    
  }
  
  rho.est<-mean(jackknife[,1])
  rho.err=sd(jackknife[,1])*sqrt(no.blocks+1)
  
  gcp.pm<-WeightedMean(grid,gcp.likelihood)
  gcp.pse<-sqrt(WeightedMean(grid^2,gcp.likelihood)-gcp.pm^2)
  
  h2.zscore.1<-mean(jackknife[,5])/sd(jackknife[,5])/sqrt(no.blocks+1)
  h2.zscore.2<-mean(jackknife[,6])/sd(jackknife[,6])/sqrt(no.blocks+1)
  
  lcv.output<-data.frame(zscore,gcp.pm,gcp.pse,rho.est,rho.err,
                         pval.fullycausal.1,pval.fullycausal.2,h2.zscore.1,h2.zscore.2)
    
  return(lcv.output)
  
  
  
  
}