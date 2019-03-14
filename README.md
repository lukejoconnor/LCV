# LCV (Latent Causal Variable model)
LCV is a method for inferring genetically causal relationships using GWAS data.

LCV is implemented in Matlab and R. In order to run LCV, you will need LD scores (non-stratified, with ancestry matching your GWAS data), which can be downloaded [here] (https://data.broadinstitute.org/alkesgroup/LDSCORE/). You can also compute your own LD scores using the [LDSC software](https://github.com/bulik/ldsc). You will also need signed summary statistics: either effect size estimates (in units of per-normalized-genotype effect size) or Z scores.

Usage of each function is described within the source code. There are example simulation scripts in Matlab and R, and an example script to run on real data in R.

## Contents:

Matlab/

example_script.m: Example script to generate data under the LCV model and run LCV on the  simulated data

simulate_LCV.m: Simulates causal effect sizes and summary statistics under the LCV model

run_LCV.m: Runs LCV on summary statistics for two traits

estimate_k4.m, weighted_mean.m, weighted_regression.m: subroutines of run_LCV.m

run_LCV_parallel.m: Runs LCV on summary statistics for two traits with parallelization across jackknife blocks
estimate_k4.m, weighted_mean.m, weighted_regression.m: Functions to compute sample moments used by LCV

R/

RunLCV.R: Runs LCV on summary statistics for two traits. Calls functions defined within MomentFunctions.R as subroutines.

MomentFunctions.R: Functions to compute sample moments used by LCV

ExampleRealdataScript.R: Example script to run LCV on real data

ExampleSimulationScript.R: Example simulations script

SimulateLCV.R: Generates simulated summary statistics following the LCV model

## Details and potential issues

1. The summary statistics and LD scores must be sorted by genomic position, as LCV uses a block-jackknife procedure to compute standard errors; if consecutive SNPs are not approximately contiguous, standard errors will be underestimated.
2. Your datasets should have approximately the same ancestry as each other, and with the LD scores. For example, it would be fine to use one UK Biobank dataset and one dataset which is a European meta-analysis, but don't try to use a European dataset with an East Asian one.
3. We recommend using SNPs with allele frequency greater than 0.05; adding additional SNPs will probably cause decreased power unless you assign them lower regression weights.
4. We recommend removing the MHC region in all analyses.

## Changes in most recent update
1. A bug in the R implementation (specifically, the WeightedRegression function) was fixed. Previously this implementation would give incorrect estimates whenever the weights vector is not uniform.
2. The sign of the Z scores no longer depend on the sign of the genetic correlation. If you run RunLCV(LDscores,Z.1,Z.2), you will now get the same Z score (but opposite genetic correlation) as if you run RunLCV(LDscores,Z.1,-Z.2).
3. Error handling and warnings are now the same for both implementations.
4. The Matlab implementation now outputs a single data structure rather than a long list of output arguments. The R data structure output was also modified.
5. There is a new R example script - thank you to Katie Siewert for supplying it.
6. Removed run_LCV_parallel.m because it only runs a few times faster.

Reference:

O'Connor, L.J. and A.L. Price. "Distinguishing genetic correlation from causation across 52 diseases and complex traits." Nature genetics (2018).

Non-paywalled link: https://rdcu.be/bajzC
