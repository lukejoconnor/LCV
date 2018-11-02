function     [ Z1,Z2, beta1,beta2, pi, gamma1, gamma2 ] = simulate_LCV( R,n1,n2,q1,q2,h2g1,h2g2,rho_e,pleiotrop,p_causalsnp,p_trait1_causal,p_trait2_causal,Rsqrt)
%SIMULATE_LCV simulates causal effect sizes and summary
%statistics under the LCV model
%   INPUT PARAMTERS: R, MxM LD matrix; n1, sample size for trait 1; n2, sample
%   size for trait 2; q1, effect of L on trait 1; q2, effect of L on trait 2; h2g1,
%   heritability of trait 1; h2g2, heritability of trait 2; rho_e, total 
%   phenotypic correlation between traits times percent sample overlap; 
%   pleiotrop, fraction of SNPs with uncorrelated pleiotropic effects on 
%   both traits; p_causalsnp, fraction of SNPs causal for L;
%   p_trait1_causal, fraction of SNPs affecting trait 1 only;
%   p_trait2_causal, fraction of SNPs affecting trait 2 only; Rsqrt
%   (optional), matrix square root of R.
%   OUTPUT PARAMETERS: Z1, estimated marginal per-normalized-genotype
%   effect sizes on trait 1; Z2, estimated marginal per-normalied-genotype 
%   effect sizes on trait 2; beta1, true causal effect sizes on trait 1; 
%   beta2, true causal effect sizes on trait 2; pi, true causal effect
%   sizes on L; gamma1, true causal effect sizes on trait 1 only; gamma2,
%   true causal effect sizes on trait 1 only.

m=size(R,1);

if any(diag(R)~=1)
    error('R should be a correlation matrix')
end
if any([q1,q2].^2>1)
    error('q1,q2 should be between -1 and 1')
end
if any([h2g1,h2g2,rho_e]>1) || any([h2g1,h2g2]<0) || rho_e<-1
    error('h2g should be between 0 and 1, and rho_e shouldbe between -1 and 1')
end

if pleiotrop+p_causalsnp+p_trait1_causal+p_trait2_causal>1
    error('Mixing weights must add up to at most 1')
end
if min([p_causalsnp,pleiotrop+p_trait1_causal,pleiotrop+p_trait2_causal])<=0
    error('Mixing weights must be positive')
end


samp=randsample(1:m,ceil(m*p_causalsnp));
temp=zeros(m,1);temp(samp)=1;
pi=randn(m,1).*temp;
pi=(pi-mean(pi));
pi=pi/sqrt(sum(pi.^2));

causalsnps_both=rand(m,1)<pleiotrop;
temp=rand(m,1);
causalsnps_trait1=temp<p_trait1_causal;
causalsnps_trait2=temp>1-p_trait2_causal;

causalsnps_trait1=(causalsnps_trait1+causalsnps_both)>0;
causalsnps_trait2=(causalsnps_trait2+causalsnps_both)>0;

gamma1=randn(m,1).*causalsnps_trait1;
gamma2=randn(m,1).*causalsnps_trait2;

gamma1(gamma1~=0)=gamma1(gamma1~=0)-mean(gamma1(gamma1~=0));
gamma1=gamma1*sqrt(1-q1^2)/sqrt(sum(gamma1.^2));
gamma2(gamma2~=0)=gamma2(gamma2~=0)-mean(gamma2(gamma2~=0));
gamma2=gamma2*sqrt(1-q2^2)/sqrt(sum(gamma2.^2));

if exist('Rsqrt')
    sumstat_noise=Rsqrt*mvnrnd([0 0],[[1/n1 rho_e/sqrt(n1*n2)];[rho_e/sqrt(n1*n2) 1/n2]],m);
else
    sumstat_noise=sqrtm(R)*mvnrnd([0 0],[[1/n1 rho_e/sqrt(n1*n2)];[rho_e/sqrt(n1*n2) 1/n2]],m);
end

Z1=R*(q1*pi+gamma1)*sqrt(h2g1)+sumstat_noise(:,1);
Z2=R*(q2*pi+gamma2)*sqrt(h2g2)+sumstat_noise(:,2);

beta1=(q1*pi+gamma1)*sqrt(h2g1);
beta2=(q2*pi+gamma2)*sqrt(h2g2);
end

