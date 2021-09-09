### This code was written by: Brian Halstead ###
### If you have any questions, please email: bhalstead@usgs.gov ###



### JAGS model code for full multi-scale occupancy model with variable selection and time-to-detection ###

model{
  ### Priors and constraints
  # site occupancy
  mu_psi ~ dbeta(1, 1) # mean site occupancy
  delta0 <- logit(mu_psi) # convert to logit-scale intercept
  sigma_delta ~ dt(0, 1, 1)T(0,) # SD for hierarchical shrinkage prior on site occupancy coefficients
  tau_delta <- pow(sigma_delta, -2) # convert SD to precision
  for(k in 1:n.delta){ # loop over k predictors for site occupancy
    delta[k] ~ dt(0, tau_delta, 1) # coefficients
    v[k] ~ dbern(0.5) # indicator variables
  }
  for(i in 1:n.site){ # loop over i sites
    logit(psi[i]) <- delta0 + sum(v[1:4] * delta[1:4] * psi_pred[i,1:4]) +
      v[4] * v[5] * delta[5] * psi_pred[i,5] # linear deterministic component of model for site occupancy (coefficient 5 is quadratic effect of slope constrained to include linear effect in model)
  }
  # pool occupancy
  mu_theta ~ dbeta(1, 1) # mean pool occupancy (conditional on site occupancy)
  gamma0 <- logit(mu_theta) # convert to logit-scale intercept
  sigma_gamma ~ dt(0, 1, 1)T(0,) # SD for hierarchical shrinkage prior on pool occupancy coefficients
  tau_gamma <- pow(sigma_gamma, -2) # convert SD to precision
  for(k in 1:n.gamma){ # loop over k predictors for pool occupancy
    gamma[k] ~ dt(0, tau_gamma, 1) # coefficients
    u[k] ~ dbern(0.5) # indicator variables
  }
  for(i in 1:n.site){ # loop over i sites
    for(j in 1:n.pool[i]){ # loop over j pools within site i
      # priors for missing predictor variables
      theta_pred[i,j,1] ~ dnorm(0, 1) 
      theta_pred[i,j,2] ~ dnorm(0, 1)
      theta_pred[i,j,3] ~ dnorm(0, 1)
      theta_pred[i,j,4] ~ dnorm(0, 1)
      theta_pred[i,j,5] ~ dnorm(0, 1)
      theta_pred[i,j,6] ~ dnorm(0, 1)
      theta_pred[i,j,7] ~ dnorm(0, 1)
      # linear deterministic component of model for pool occupancy
      logit(theta[i,j]) <- gamma0 + sum(u[] * gamma[] * theta_pred[i,j,])
    }
  }
  ## time-to-detection
  mu_mu ~ dunif(1, msd) # mean time to detection
  beta0 <- log(mu_mu) # convert to intercept on log scale
  sigma_beta ~ dt(0, 1, 1)T(0,) # SD for hierarchical shrinkage prior on
  # time-to-detection coefficients
  tau_beta <- pow(sigma_beta, -2) # convert SD to precision
  for(k in 1:n.beta){ # loop over k predictors for time-to-detection
    beta[k] ~ dt(0, tau_beta, 1) # coefficients
    o[k] ~ dbern(0.5) # indicator variables
  }
  for(i in 1:n.site){ # loop over i sites
    for(j in 1:n.pool[i]){ # loop over j pools in site i
      # priors for missing predictor variables
      mu_pred[i,j,1] ~ dnorm(0, 1)
      mu_pred[i,j,2] ~ dnorm(0, 1)
      mu_pred[i,j,3] ~ dnorm(0, 1)
      mu_pred[i,j,4] ~ dnorm(0, 1)
      # linear deterministic part of time-to-detection model (coefficient 6 is quadratic effect of date constrained to include linear effect in model)
      log(mu[i,j]) <- beta0 + sum(o[1:5] * beta[1:5] * mu_pred[i,j,1:5]) +
        o[5] * o[6] * beta[6] * mu_pred[i,j,6]
      lambda[i,j] <- 1 / mu[i,j] # convert time-to-detection to detection rate
    }
  }
  ### Multi-scale occupancy model
  ## Occupancy of region
  for(i in 1:n.site){ # loop over i sites
    z[i] ~ dbern(psi[i]) # site i occupied with probability psi[i]
    ## Use of pool for breeding
    for(j in 1:n.pool[i]){ # loop over j pools in site i
      theta_eff[i, j] <- theta[i,j] * z[i] # effective probability of pool occupancy (i.e., conditional on occurrence at site i)
      w[i,j] ~ dbern(theta_eff[i,j]) # pool i,j occupied with probability theta_eff[i,j]
      ## Time to detection
      # exponentially distributed detection times
      t.det[i,j] ~ dexp(lambda[i,j]) # time to detection at pool j in site i
      # model for censoring of data (i.e., detection did not occur during survey) 
      omega[i,j] <- w[i,j] * step(t.det[i,j] - t.survey[i,j]) + (1 - w[i,j])
      d[i,j] ~ dbern(omega[i,j]) # time-to-detection data censored at pool i,j with probability omega[i,j]
      ## GOF ##
      w.new[i,j] ~ dbern(theta_eff[i,j]) # simulated pool occupancy
    }
    eval[i] <- sum(theta_eff[i,1:n.pool[i]]) # expected number of occupied pools at site i
    n_occ_obs[i] <- sum(w[i,1:n.pool[i]]) # number of observed pools estimated to be occupied at site i
    E[i] <- pow((n_occ_obs[i] - eval[i]), 2) / (eval[i] + 0.5) # chi-square observed
    n_occ_sim[i] <- sum(w.new[i,1:n.pool[i]]) # number of simulated pools estimated to be occupied at site i
    E.new[i] <- pow((n_occ_sim[i] - eval[i]), 2) / (eval[i] + 0.5) # chi-# square simulated
    ## derived parameter: finite sample inference for number of sampled pools occupied at each site
    n_occ_pool[i] <- sum(w[i,1:n.pool[i]]) # estimated number of occupied pools at site i
  }
  ## GOF (continued)
  fit <- sum(E[]) # sum of chi-square observed
  fit.new <- sum(E.new[]) # sum of chi-square simulated
  test <- step(fit.new - fit) # simulated data deviate from expectations more than observed data?
  ## derived parameter: finite sample inference for number of sampled sites occupied
  n_occ_site <- sum(z[]) 
} # model stop
