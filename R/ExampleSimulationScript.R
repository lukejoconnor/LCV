# Example simulations script calling RunLCV.R

source("SimulateLCV.R")
source("RunLCV.R")
## Simulation parameters

M<-50000 # total number of SNPs

N.1<-20000 # sample size Y1 (note: this script doesn't actually generate individual-level data; set N.k as large as you want)
N.2<-50000 # sample size Y2

h2.1<-.3 # Y1 heritability
h2.2<-.3 # Y2 heritability

q.1<-1 # effect of L on Y1
q.2<-0.2 # effect of L on Y2

gcp.true<-log(q.2/q.1)/log(q.2*q.1) 
rho.true<-q.1*q.2 
sprintf("gcp=%.2f, rho=%.2f\n",gcp.true, rho.true)

p.pi<-.05 # Proportion of SNPs causal for L. p.pi*M should be an integer.
p.g1<-.05 # Proportion of SNPs causal for Y1 only
p.g2<-.2 # Proportion of SNPs causal for Y2 only

## LCV parameters
crosstrait.intercept<-0 # Whether to use intercept in cross-trait LDSC regression. Should be 0 if there is no LD (generally, recommend setting this to 1 in real-data analyses)

ldsc.intercept<-0 # Whether to use intercept in univariate LDSC regression. Should be 0 if there is no LD (generally, recommend setting this to 1 in real-data analyses)

# Significance threshold: throw out SNPs about this threshold times mean
#   chisq for computing LDSC intercept (all SNPs are included in other
#   parts of the analysis)
sig.threshold<-30

# Number of jackknife blocks
no.blocks<-100

# Value of crosstrait LDSC intercept. Only applicable if crosstrait.intercept==1.
# Should be zero for non-overlapping cohorts.
intercept.12<-0;

ell<-rep(1,M) # LD scores
weights <- rep(1,M) # regression weights

## Simulate from LCV model
alpha<-SimulateLCV(M,N.1,N.2,h2.1,h2.2,q.1,q.2,p.pi,p.g1,p.g2)

## Run LCV

LCV<- RunLCV(ell,alpha$a1,alpha$a2,no.blocks,crosstrait.intercept,ldsc.intercept,weights,sig.threshold,N.1,N.2,intercept.12)


sprintf("Estimated posterior gcp=%.2f(%.2f), log10(p)=%.1f; estimated rho=%.2f(%.2f)",LCV$gcp.pm, LCV$gcp.pse, log(LCV$pval.gcpzero.2tailed)/log(10), LCV$rho.est, LCV$rho.err)

