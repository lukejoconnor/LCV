% Example script to generate data under the LCV model and run LCV on the
% simulated data. These simulations generate simulations directly 
% (no individual-level genotypes).

%% Simulation parameters

n1=2*10^4;% Sample size for trait 1
n2=10^5;
M=5*10^4;% No. snps

R=speye(M);% LD matrix (no LD)
Rsqrt=R;% Matrix square root of LD matrix
ell=sum(R.^2)';% Vector of LD scores
weights=1./ell;% Regression weights

q1=1;% Effect of L on trait 1
q2=.1;
gcp_true=log(q2/q1)/log(q2*q1);% True value of gcp
rhog_true=q1*q2;% True value of rho_g
fprintf('gcp=%.2f, rhog=%.2f\n',gcp_true,rhog_true)

h2g1=.3;% Trait 1 heritability
h2g2=.3;

% Correlation between noise terms in summary stats for each trait. 
%   Should be 0 if crosstrait_intercept is 0.
rho_e=0;

p_pi=.05;% Proportion of snps causal for L (N/A if q1=q2=0)
p_g1=.05;% Proportion of SNPs causal for trait 1 only (N/A if q1=1)
p_g2=.2;
p_pleiotrop=0;% Proportion of SNPs with uncorrelated pleiotropic effects on both traits

% Whether to use intercept in cross-trait LDSC regression. Should be 0 if there is no LD.
crosstrait_intercept=0;
% Whether to use intercept in LDSC regression. Should be 0 if there is no LD.
ldsc_intercept=0;
% Significance threshold: throw out SNPs about this threshold times mean
%   chisq for computing LDSC intercept (all SNPs are included in other
%   parts of the analysis)
sig_threshold=30;
% Number of jackknife blocks
no_blocks=100;
% Correlation between noise terms in summary stats, if known (only
%   applicable if crosstrait_intercept=2)
cross_int=0;



%% Simulate data
[ Z1,Z2, beta1,beta2, pi, gamma1, gamma2 ] = ...
    simulate_LCV( R,n1,n2,q1,q2,h2g1,h2g2,rho_e,p_pleiotrop,p_pi,p_g1,p_g2,Rsqrt);

%% Run LCV

% Last 3 arguments not used if crosstrait_intercept=ldsc_intercept=1
 [ zsc_asym,gcp_est,gcp_err,rho,rho_err] = ...
    run_LCV( ell,Z1,Z2,crosstrait_intercept,ldsc_intercept,weights,sig_threshold,...
    no_blocks,cross_int,n1,n2);

LCV_pval=tcdf(-abs(zsc_asym),no_blocks-2)*2;% Two-tailed

fprintf('P-value for partial causality: %.2g\n',LCV_pval)
fprintf('gcp estimate: %.2f(%.2f)\n',gcp_est,gcp_err)