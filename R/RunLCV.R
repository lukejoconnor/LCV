# RunLCV runs LCV on summary statistics for two traits.
#	Input vectors are sorted lists of LD scores and signed summary statistics. They
#		must be sorted by genomic position, as LCV uses a block-jackknife procedure
#		to compute standard errors; if consecutive SNPs are not approximately
#		contiguous, standard errors will be underestimated.
#
#   INPUT VARIBLES: 
#	ell, Mx1 vector of LD scores; 
#	z.1, Mx1 vector of estimated marginal per-normalized-genotype effects on trait 1
#   	(or Z scores; invariant to scaling);  
#	z.2, Mx2 vector of effects on trait2;  
#   no.blocks, number of jackknife blocks;
#	crosstrait.intercept, 0 if crosstrait LDSC intercept should be fixed and 1 otherwise;
#   ldsc.intercept, 0 if LDSC intercept should be fixed and 1 otherwise (1 is recommended);
#   weights, Mx1 vector of regression weights; 
#	sig.threshold, threshold above which to discard chisq statistics for the purpose of estimating
#   	the LDSC intercept; large-effect SNPs discarded above sig_threshold*mean(z.x^2);
#   n1, sample size for trait 1, only needed if ldsc.intercept=1; 
#   n2, sample size for trait 2, only needed if ldsc.intercept=1; 
#	intercept.12, covariance between sampling errors for Z1 and Z2, only needed if
#   	crosstrait_intercept=0. If the 2 GWAS are disjoint, this can be set to zero.
#
#   OUTPUT VARIABLES: 
#	lcv.output, a list with named entries:
#   "zscore", Z score for partial genetic causality. zscore>>0 implies gcp>0.
#	pval.gcpzero.2tailed, 2-tailed p-value for null that gcp=0.
#   "gcp.pm", posterior mean gcp (gcp=1: trait 1 -> trait 2; gcp=-1: trait 2-> trait 1); 
#   "gcp.pse", posterior standard error for gcp; 
#   "rho.est", estimated genetic correlation; 
#   "rho.err", standard error of rho estimate;
#   "pval.fullycausal" [2 entries], p-values for null that gcp=1 or that gcp=-1, respectively; 
#   "h2.zscore" [2 entries], z scores for trait 1 and trait 2 being heritable, respectively;
#   	we recommend reporting results for h2.zscore > 7 (a very stringent threshold).

RunLCV <- function(ell,z.1,z.2,no.blocks=100,crosstrait.intercept=1,ldsc.intercept=1,weights=1/pmax(1,ell),sig.threshold=.Machine$integer.max,n.1=1,n.2=1,intercept.12=0){
  mm=length(ell)
  if(length(z.1)!=mm||length(z.2)!=mm){
    stop('LD scores and summary statistics should have the same length')
  }
  
  source("MomentFunctions.R")
  grid<- (-100:100)/100;
  
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
    if(any(is.nan(jackknife))){
      stop('NaNs produced, probably due to negative heritability estimates. Check that summary statistics and LD scores are ordered correctly.')
    }
  }
  rho.est<-mean(jackknife[,1])
  rho.err=sd(jackknife[,1])*sqrt(no.blocks+1)
  flip=sign(rho.est)
  
  jackknife[,2:3]<-jackknife[,2:3]-3*jackknife[,1]
  
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
      pval.fullycausal.2<-pt(-flip*mean(pct.diff)/sd(pct.diff)/sqrt(no.blocks+1),no.blocks-2)
    }
    if(xx==1){
      pval.fullycausal.1<-pt(flip*mean(pct.diff)/sd(pct.diff)/sqrt(no.blocks+1),no.blocks-2)
    }
    if(xx==0){
      zscore<- flip*mean(pct.diff)/sd(pct.diff)/sqrt(no.blocks+1)
    }
    
  }
  
  pval.gcpzero.2tailed=pt(-abs(zscore),no.blocks-1)*2
  
  gcp.pm<-WeightedMean(grid,gcp.likelihood)
  gcp.pse<-sqrt(WeightedMean(grid^2,gcp.likelihood)-gcp.pm^2)
  
  h2.zscore.1<-mean(jackknife[,5])/sd(jackknife[,5])/sqrt(no.blocks+1)
  h2.zscore.2<-mean(jackknife[,6])/sd(jackknife[,6])/sqrt(no.blocks+1)
  
  if(h2.zscore.1<4 || h2.zscore.2<4){
    warning('Very noisy heritability estimates potentially leading to false positives')
  }
  else{
    if(h2.zscore.1<7 || h2.zscore.2<7){
      warning('Borderline noisy heritability estimates potentially leading to false positives')
    }
  }
  if(abs(rho.est/rho.err)<2){
    warning('No significantly nonzero genetic correlation, potentially leading to conservative p-values')
  }
  
  lcv.output<-list(zscore=zscore,pval.gcpzero.2tailed=pval.gcpzero.2tailed,gcp.pm=gcp.pm,gcp.pse=gcp.pse,rho.est=rho.est,rho.err=rho.err,
                         pval.fullycausal=c(pval.fullycausal.1,pval.fullycausal.2),h2.zscore=c(h2.zscore.1,h2.zscore.2))
  
  return(lcv.output)
  
  
  
  
}