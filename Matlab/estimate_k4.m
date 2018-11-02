function [ rho,k41,k42,intercept12,s1,s2,intercept1,intercept2 ] = estimate_k4( ell,Z1,Z2,crosstrait_intercept,ldsc_intercept,weights,sig_threshold,n1,n2,noisecorr  )
%ESTIMATE_K4 estimates mixed 4th moments and other parameters from summary
%statistics
%   INPUT VARIBLES: ell, Mx1 vector of LD scores; Z1, Mx1 vector of
%   estimated marginal per-normalized-genotype effects on trait 1
%   (or Z scores; invariant to scaling); Z2, Mx2 vector of effects on trait
%   2; crosstrait_intercept, 0 if cohorts are disjoint, 1 if cohorts are
%   possibly nondisjoint and necessary correction is unknown, 2 if cohorts
%   are nondisjoint with known overlap and phenotypic covariance;
%   ldsc_intercept, 0 if intercept should be fixed and 1 otherwise;
%   weights, Mx1 vector of regression weights.
%   OPTIONAL INPUT VARIBLES: sig_threshold, threshold
%   above which to discard chisq statistics for the purpose of estimating
%   the LDSC intercept if they are above sig_threshold*mean_chisq; n1,
%   1/var(Z1), only needed if ldsc_intercept=1; n2, 1/var(Z2); noisecorr,
%   covariance between sampling errors for Z1 and Z2, only needed if
%   crosstrait_intercept=2 (if zero, equivalent to setting
%   crosstrait_intercept=0).
%   OUTPUT VARIABLES: rho, estimated genetic correlation; k41, estimate of
%   E(alpha1^3 alpha2); k42, estimate of E(alpha2^3 alpha1); intercept12,
%   crosstrait LDSC intercept; s1, normalization used for Z1, proportional
%   to sqrt(h2g); s2, normalization used for Z2; intercept1, LDSC intercept
%   for trait 1; intercept2, LDSC intercept for trait 2.

h2g_warning_tol=.02;
if ~exist('sig_threshold'); sig_threshold=inf; end% For excluding SNPs from LDSC

m=length(Z1);

% LDSC regression on each trait
if ~ldsc_intercept
    intercept1=1/n1;
    intercept2=1/n2;
    temp=weighted_regression([ell],Z1.^2-intercept1,weights);
    h2g1=temp(1);
    temp=weighted_regression([ell],Z2.^2-intercept2,weights);
    h2g2=temp(1);
else
    sig_snps1=Z1.^2>sig_threshold*mean(Z1.^2);
    temp=weighted_regression([ell(~sig_snps1), ones(sum(~sig_snps1),1)],Z1(~sig_snps1).^2,weights(~sig_snps1));
    intercept1=temp(end);
    temp=weighted_regression([ell],Z1.^2-intercept1,weights);
    h2g1=temp(1);
    sig_snps2=Z2.^2>sig_threshold*mean(Z2.^2);
    temp=weighted_regression([ell(~sig_snps2), ones(sum(~sig_snps2),1)],Z2(~sig_snps2).^2,weights(~sig_snps2));
    intercept2=temp(end);
    temp=weighted_regression([ell],Z2.^2-intercept2,weights);
    h2g2=temp(1);    
end

% Cross-trait LDSC regression
if crosstrait_intercept==0
    temp=weighted_regression([ell],(Z1.*Z2),weights);
    rho=temp/sqrt(h2g1*h2g2);
    intercept12=0;
elseif crosstrait_intercept==1
    sig_snps=1-(Z1.^2<sig_threshold*mean(Z1.^2)).*(Z2.^2<sig_threshold*mean(Z2.^2))==1;
    temp=weighted_regression([ell(~sig_snps), ones(sum(~sig_snps),1)],(Z1(~sig_snps).*Z2(~sig_snps)),weights(~sig_snps));
    intercept12=temp(end);
    temp=weighted_regression([ell],(Z1.*Z2-intercept12),weights);
    rho=temp(1)/sqrt(h2g1*h2g2);
elseif crosstrait_intercept==2
    intercept12=noisecorr;
    temp=weighted_regression([ell],(Z1.*Z2-intercept12),weights);
    rho=temp/sqrt(h2g1*h2g2);
end

% Normalize effect sizes 
s1=sqrt(weighted_mean(Z1.^2,weights)-intercept1);
s2=sqrt(weighted_mean(Z2.^2,weights)-intercept2);
nZ1=Z1/s1;
nZ2=Z2/s2;

% Estimates of mixed 4th moments
k41=weighted_mean(nZ2.*nZ1.^3-3*nZ1.*nZ2*(intercept1/s1^2)-3*(nZ1.^2-intercept1/s1^2)*intercept12/s1/s2,weights);%/sum(nZ1.*nZ2)^2;
k42=weighted_mean(nZ1.*nZ2.^3-3*nZ1.*nZ2*(intercept2/s2^2)-3*(nZ2.^2-intercept2/s2^2)*intercept12/s1/s2,weights);


end

