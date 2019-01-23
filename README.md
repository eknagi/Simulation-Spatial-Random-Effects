# Simulating-Spatial-Dependency-Random-Effects

The motivation for this topic arose due to my interest in scenarios where multilevel models could not account
for the independence of residuals: ie, when the ICC is not same for two individual members of a group. This
can happen due to different distances between members in an area associated with a group. Additionally,
correlation does not necessarily stop at boundaries, but can spill over. Therefore, modelling such data as
though they are independent will lead to biased standard errors.
My goal is to simulate datasets with multilevel random effects and spatial dependencies incorporated; and
assess and compare the fit and parameter estimation of different (pure) multilevel models, (pure) spatial
models, and models that are a hybrid between the two.


Simulation Analysis: Distance-based (strength of correlation between two locations decreases as the distance between them
increases)

The simulation will consist of 2-levels: one at the individual level and the other at a block level.
1) “N” points from an exponential spatial correlation structure (Gaussian random field) are simulated on
a square grid.
2) A 2-level multilevel model is simulated with random effects for the slope and intercept: One individuallevel (N observations) covariate “x1” and a block-level covariate “x2” (for i blocks) are generated from
normal distributions respectively. Random intercept are incorporated for block-level effects and random
slope for the block-level covariate was also generated from a normal distribution according to the
specifications below.
3) For each observation, the fixed effect covariates, random effects, and spatial effects are added to arrive
at the outcome variable as specified below:
Yij = b0 + b1jx1ij + b2x2j + ζj + ηij + ij
where b0j = b0 + ζj is the random intercept, b1j is the random slope, ηij is the spatial correlation term,
x1ij ∼ N(0, 1), x2j ∼ N(0, 1), ζj ∼ N(0, 1), (block-level random intercept effect), ηij ∼ N(0, σ2η) (spatial
correlation)
εijk is the residual term (independent of the random effect term and spatial correlation term above)
x1 is individual-level covariate x2 is group covariate (same value for every person in a group) b0 is true
intercept, b1 is true slope of individual covariate, b2 is the true slope for the group covariate 


Run different hybrids of multilevel and spatial models to assess the level of parameter estimation.
