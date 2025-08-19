library(brms)
library(ape)

d_pca = read_csv("Data/d_pca.csv") #pca scores for native and non-native urban plant species

# scale variables
d_pca$ave_tmean = scale(d_pca$ave_tmean)[,1]
d_pca$ave_precip = scale(d_pca$ave_precip)[,1]

phy = read.tree("Data/phy.tre")#phylogeny for species in the PCA data set


# Step 1: Create the phylogenetic correlation matrix
phylo_cov_matrix <- ape::vcv(phy, corr = TRUE)


priors <- c(
  # Regression coefficients for both mixture components
  prior(normal(0, 2), class = "b", dpar = "mu1"),
  prior(normal(0, 2), class = "b", dpar = "mu2"),
  
  # Intercepts for both components
  prior(normal(0, 5), class = "Intercept", dpar = "mu1"),
  prior(normal(0, 5), class = "Intercept", dpar = "mu2"),
  
  # Skewness parameters
  prior(normal(0, 2), class = "alpha", dpar = "alpha1"),
  prior(normal(0, 2), class = "alpha", dpar = "alpha2"),
  
  # Standard deviations
  prior(exponential(1), class = "sigma", dpar = "sigma1"),
  prior(exponential(1), class = "sigma", dpar = "sigma2"),
  
  # Mixing proportions
  prior(dirichlet(2, 2), class = "theta"),
  
  # Group-level standard deviations (phylogenetic random effect)
  prior(exponential(1), class = "sd", group = "tip", dpar = "mu1"),
  prior(exponential(1), class = "sd", group = "tip", dpar = "mu2")
)



# Step 2: Fit the Bayesian model
model_bayes <- brm(
  PC1 ~ provenance + ave_tmean + ave_precip +
    ave_tmean:provenance + ave_precip:provenance +
    I(ave_tmean^2) + I(ave_precip^2) +
    (1 | gr(tip, cov = phylo_cov_matrix)),
  data = d_pca,
  data2 = list(phylo_cov_matrix = phylo_cov_matrix),
  family = mixture(skew_normal(), skew_normal()),
  chains = 4, cores = 4, iter = 500,
  control = list(adapt_delta = 0.99, max_treedepth = 15)
)

summary(model_bayes)


brm(
  formula = ...,
  data = ...,
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  iter = 4000,
  chains = 4
)

prior_summary(model_bayes)

