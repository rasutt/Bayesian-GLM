# Code to illustrate Bayesian GLM

# Load libraries
library(MASS)
library(plot3D)

# Load, summarize, and plot data
precip_coal_14Ma <- read.csv("precip_coal_14Ma.csv")
summary(precip_coal_14Ma)
plot(Coal ~ Precip, data = precip_coal_14Ma)

# Fit, summarize, and plot fitted values of binomial GLM
coal_glm <- glm(
  formula = Coal ~ Precip, family = binomial(link="logit"),
  data = precip_coal_14Ma, method = "glm.fit"
)
summary(coal_glm)
points(precip_coal_14Ma$Precip, fitted(coal_glm), col = "red")

# Find mean and covariance of MLE distribution
beta_hat <- coef(coal_glm)
beta_cov <- vcov(coal_glm)

# Set number of samples to draw
n_samps <- 1e4

# Simulate proposals for regression coefficients using MLEs
beta_props <- mvrnorm(n = n_samps, beta_hat, beta_cov)

# Set break points for histograms
breaks <- t(beta_hat + outer(sqrt(diag(beta_cov)), -10:10 * 3 / 10))

# Plot simulated values of regression coefficients
z_props <- table(cut(beta_props[, 1], breaks[, 1]), 
                 cut(beta_props[, 2], breaks[, 2]))
hist3D(
  z = z_props, 
  border = 'black', col = 'light blue', theta = 135,
  xlab = 'beta_0', ylab = 'beta_1', zlab = 'frequency',
  main = 'Regression coefficients simulated from MLEs'
)

# Function to run MCMC with a given prior variance on the weights and plot the
# results
MCMC_plot <- function(c) {
  # Function proportional to log posterior of weights given N(0, cI) prior and
  # Bernoulli likelihood with logit link
  y = precip_coal_14Ma[["Coal"]]
  X = cbind(1, precip_coal_14Ma[["Precip"]])
  log_post_shape <- function(beta) {
    sum(X %*% beta * y - log(1 + exp(X %*% beta))) - t(beta) %*% beta / (2 * c)
  }
  
  # Function proportional to log proposal probability of weights
  log_prop_shape <- function(beta) {
    t(beta - beta_hat) %*% solve(beta_cov) %*% (beta - beta_hat)
  }
  
  # Draw uniform RVs for acceptance of MCMC proposals
  U <- runif(n = n_samps)
  
  # Create matrix for samples from posterior including burn-in period of chain
  betas <- matrix(nrow = n_samps + 2e3, ncol = 2)
  
  # Run MCMC to draw from posterior
  for (i in 2:n_samps) {
    # Find acceptance probability
    alpha <- min(
      1, 
      exp(log_post_shape(beta_props[i, ]) + log_prop_shape(beta_props[i - 1, ]) 
          - log_post_shape(beta_props[i - 1, ]) - 
            log_prop_shape(beta_props[i, ]))
    )
    
    # If draw from uniform RV smaller than acceptance probability take step,
    # otherwise repeat sample
    if (U[i] < alpha) betas[i, ] <- beta_props[i, ]
    else betas[i, ] <- betas[i - 1, ]
  }
  
  # Discard samples from burn-in period of chain
  betas <- betas[-(1:2e3), ]
  
  # Plot distributions of weights
  z_post <- table(cut(betas[, 1], breaks[, 1]), 
                  cut(betas[, 2], breaks[, 2]))
  hist3D(
    z = z_props, 
    border = 'black', col = 'light blue', theta = 135,
    xlab = 'beta_0', ylab = 'beta_1', zlab = 'frequency',
    main = 'Simulated regression coefficients', sub = paste('c =', c),
    zlim = c(0, max(z_props, z_post))
  )
  hist3D(z = z_post, border = 'black', col = 'light pink', add = T)
  legend("topleft", legend = c('MLE', 'Posterior'), pch = rep(20, 2),
         col = c('light blue', 'light pink'))
}

# Run MCMC and plot results for priors of different strengths determined by
# their variances

# Weak prior gives posterior similar to MLE
MCMC_plot(c = 100)

# Strong prior means posterior shifted from MLE
MCMC_plot(c = 0.001)
