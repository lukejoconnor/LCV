# LCV
Preliminary software implementing the Latent Causal Variable Model in Matlab and R. Usage of each function is described within the source code. If you have any questions, feel free to contact Luke O'Connor at loconnor@g.harvard.edu.

Contents:

=======
Matlab/
example_script.m: Example script to generate data under the LCV model and run LCV on the  simulated data

simulate_LCV.m: Simulates causal effect sizes and summary statistics under the LCV model

run_LCV.m: Runs LCV on summary statistics for two traits

estimate_k4.m, weighted_mean.m, weighted_regression.m: subroutines of run_LCV.m

run_LCV_parallel.m: Runs LCV on summary statistics for two traits with parallelization across jackknife blocks
estimate_k4.m, weighted_mean.m, weighted_regression.m: Functions to compute sample moments used by LCV

R/
RunLCV.R: Runs LCV on summary statistics for two traits
MomentFunctions.R: Functions to compute sample moments used by LCV

Reference:
O'Connor, L.J. and A.L. Price. "Distinguishing genetic correlation from causation across 52 diseases and complex traits." Nature genetics (2018).

Non-paywalled link: https://rdcu.be/bajzC
