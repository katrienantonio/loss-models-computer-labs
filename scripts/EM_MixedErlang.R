######################################################################
## EM algorithm mixtures of Erlangs for censored and truncated data ##
######################################################################

# Date: 16/12/2014
# Author: Roel Verbelen

## Packages

packages <- c("MASS", "survival", "doParallel")
packages <- lapply(packages, FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
        install.packages(x)
        library(x, character.only = TRUE)
    }
})

## Initial values

# supply number of shapes M and spread factor s

ME_initial <- function(lower, upper, trunclower=0, truncupper=Inf, M=10, s=1){
  # data for initial step: treat left and right censored data as observed and take mean for interval censored data
  uc <- lower[lower==upper & !is.na(lower==upper)]
  lc <- upper[is.na(lower)]
  rc <- lower[is.na(upper)]
  ic <- (lower[lower!=upper & !is.na(lower!=upper)] + upper[lower!=upper & !is.na(lower!=upper)]) / 2
  initial_data = c(uc, lc, rc, ic)
  theta <- max(initial_data)/M
  shape <- seq(1, M)
  alpha <- rep(0,M)
  for (i in 1:M){
    alpha[i] <- sum(initial_data <= i*theta & initial_data > (i-1)*theta)
  }    
  shape <- shape[alpha>0]
  alpha <- alpha[alpha>0]/sum(alpha)
  # spread out shapes and adjust theta
  shape <- s*shape
  theta <- theta/s
  # alpha to beta
  t_probabilities <- pgamma(truncupper, shape, scale=theta) - pgamma(trunclower, shape, scale=theta)  
  beta = alpha * t_probabilities / sum(alpha*t_probabilities)
  list(theta=theta, shape=shape, alpha=alpha, beta=beta) 
}

## Log likelihood

ME_loglikelihood <- function(x_densities, c_probabilities, beta, t_probabilities, no_censoring, censoring){
  likelihood_contribution = numeric(0)
  if(no_censoring){  
    # matrix containing alpha*density (uncensored)
    x_components <-  t(t(x_densities)*beta/t_probabilities)
    # likelihood contribution (uncensored)
    likelihood_contribution <- rowSums(x_components) 
  }   
  if(censoring){  
    # matrix containing alpha*probabilities (censored)
    c_components <- t(t(c_probabilities)*beta/t_probabilities)
    # likelihood contribution (censored)
    likelihood_contribution <- c(likelihood_contribution, rowSums(c_components)) 
  }   
  loglikelihood_contribution <- ifelse(likelihood_contribution>0, log(likelihood_contribution), -1000)
  # loglikelihood
  sum(loglikelihood_contribution)
}

## ^{u}z_{ij}^{(k)}: posterior probabilities (uncensored)

ME_u_z <-function(x_densities, beta, t_probabilities, M){      
  x_components <- t(t(x_densities)*beta/t_probabilities)
  u_z <- x_components / rowSums(x_components)
  # in case all ^{u}z_{ij}^{k} for j=1,...,M are numerically 0
  u_z[is.nan(u_z)] = 1/M
  u_z
}

## ^{c}z_{ij}^{(k)}: posterior probabilities (censored)

ME_c_z <-function(c_probabilities, beta, t_probabilities, M){      
  c_components <- t(t(c_probabilities)*beta/t_probabilities)
  c_z <- c_components / rowSums(c_components)
  # in case all ^{c}z_{ij}^{k} for j=1,...,M are numerically 0
  c_z[is.nan(c_z)] = 1/M
  c_z
}
 
## Expected value of censored observations

ME_expected_c <-function(lower, upper, shape, theta, c_z){  
  c_exp <- theta * (outer(upper, shape+1, pgamma, scale=theta) - outer(lower, shape+1, pgamma, scale=theta))/(outer(upper, shape, pgamma, scale=theta) - outer(lower, shape, pgamma, scale=theta))
  c_exp <- t(t(c_exp)*shape)
  # replace numerical 0/0 (NaN) or Inf by correct expected value
  c_exp <- ifelse(is.nan(c_exp) | c_exp==Inf, ifelse(outer(lower, shape*theta, ">"), lower, upper), c_exp)
  c_exp <- rowSums(c_z * c_exp)  
}

## T^{(k)}: correction term due to truncation

ME_T <- function(trunclower, truncupper, shape, theta, beta){
  # avoid NaN
  if(truncupper==Inf){ 
    # take log first for numerical stability (avoid Inf / Inf)
    deriv_trunc_log_1 <- shape*log(trunclower)-trunclower/theta - (shape-1)*log(theta) - lgamma(shape) - log(1 - pgamma(trunclower, shape, scale=theta))
    deriv_trunc <- exp(deriv_trunc_log_1)    
  } else{    
    deriv_trunc_log_1 <- shape*log(trunclower)-trunclower/theta - (shape-1)*log(theta) - lgamma(shape) - log(pgamma(truncupper, shape, scale=theta) - pgamma(trunclower, shape, scale=theta))
    deriv_trunc_log_2 <- shape*log(truncupper)-truncupper/theta - (shape-1)*log(theta) - lgamma(shape) - log(pgamma(truncupper, shape, scale=theta) - pgamma(trunclower, shape, scale=theta))
    deriv_trunc <- exp(deriv_trunc_log_1)-exp(deriv_trunc_log_2)
  }
  sum(beta*deriv_trunc)
}

## Auxiliary functions used to maximize theta in the M-step

theta_nlm_u_c <- function(theta, x, c_exp, n, beta, shape, trunclower, truncupper){
  T <- ME_T(trunclower, truncupper, shape, theta, beta)
  (theta - ((sum(x)+sum(c_exp))/n-T)/sum(beta*shape))^2
}

theta_nlm_u <- function(theta, x, n, beta, shape, trunclower, truncupper){
  T <- ME_T(trunclower, truncupper, shape, theta, beta)
  (theta - (sum(x)/n-T)/sum(beta*shape))^2
}

theta_nlm_c <- function(theta, c_exp, n, beta, shape, trunclower, truncupper){
  T <- ME_T(trunclower, truncupper, shape, theta, beta)
  (theta - (sum(c_exp)/n-T)/sum(beta*shape))^2
}

## EM algorithm (censored and truncated data)

ME_em <- function(lower, upper, trunclower=0, truncupper=Inf, theta, shape, beta, eps=1e-03, print=TRUE){
  n <- length(lower)
  M <- length(shape)
  # separate uncensored and censored observations
  uncensored <- (lower==upper & !is.na(lower==upper))
  # Uncensored observations
  x <- lower[uncensored]
  # Censored observations
  lower[is.na(lower)] <- trunclower
  upper[is.na(upper)] <- truncupper
  lower <- lower[!uncensored]
  upper <- upper[!uncensored]  
  # Boolean for having uncensored and censored observations
  no_censoring <- (length(x) != 0)  
  censoring <- (length(lower) != 0)  
  iteration <- 1
  if(no_censoring){  
    # matrix containing densities (uncensored)
    x_densities <- outer(x,shape,dgamma, scale=theta)
  }   
  if(censoring){  
    # matrix containing censoring probabilities (censored)
    c_probabilities <- outer(upper,shape,pgamma, scale=theta)-outer(lower,shape,pgamma, scale=theta)
  } 
  # truncation probabilities
  t_probabilities <- pgamma(truncupper, shape, scale=theta) - pgamma(trunclower, shape, scale=theta)
  loglikelihood <- ME_loglikelihood(x_densities, c_probabilities, beta, t_probabilities, no_censoring, censoring)
  old_loglikelihood <- -Inf
  history_loglikelihood <- loglikelihood
  history_theta <- theta
  while(loglikelihood - old_loglikelihood > eps){
    old_loglikelihood <- loglikelihood
    # E step
    if(no_censoring & censoring){
      u_z <- ME_u_z(x_densities, beta, t_probabilities, M)
      c_z <- ME_c_z(c_probabilities, beta, t_probabilities, M)
      c_exp <- ME_expected_c(lower, upper, shape, theta, c_z)      
    } else if(no_censoring){
      u_z <- ME_u_z(x_densities, beta, t_probabilities, M)
    } else{
      c_z <- ME_c_z(c_probabilities, beta, t_probabilities, M)
      c_exp <- ME_expected_c(lower, upper, shape, theta, c_z)      
    }
    # M step
    if(no_censoring & censoring){
      beta <- (colSums(u_z)+colSums(c_z))/n 
      theta <- nlm(theta_nlm_u_c, theta, x, c_exp, n, beta, shape, trunclower, truncupper)$estimate
    } else if(no_censoring){
      beta <- colSums(u_z)/n   
      theta <- nlm(theta_nlm_u, theta, x, n, beta, shape, trunclower, truncupper)$estimate
    } else{
      beta <- colSums(c_z)/n     
      theta <- nlm(theta_nlm_c, theta, c_exp, n, beta, shape, trunclower, truncupper)$estimate
    }
    iteration <- iteration + 1
    if(no_censoring){  
      # matrix containing densities (uncensored)
      x_densities <- outer(x,shape,dgamma, scale=theta)
    }   
    if(censoring){  
      # matrix containing censoring probabilities (censored)
      c_probabilities <- outer(upper,shape,pgamma, scale=theta)-outer(lower,shape,pgamma, scale=theta)
    } 
     # truncation probabilities
    t_probabilities <- pgamma(truncupper, shape, scale=theta) - pgamma(trunclower, shape, scale=theta)
    loglikelihood <- ME_loglikelihood(x_densities, c_probabilities, beta, t_probabilities, no_censoring, censoring)
    if(print) print(loglikelihood)
    history_loglikelihood <- c(history_loglikelihood, loglikelihood)
    history_theta <- c(history_theta, theta)
  }
  # beta to alpha
  alpha_tilde <- beta / t_probabilities
  alpha <- alpha_tilde / sum(alpha_tilde)
 list(alpha = alpha, beta = beta, shape = shape, theta = theta, loglikelihood = loglikelihood, history_loglikelihood = history_loglikelihood, iteration = iteration, AIC=-2*loglikelihood+2*(2*length(alpha)),BIC=-2*loglikelihood+(2*length(alpha))*log(n), history_theta=history_theta) 
}

## Shape adjustments

ME_shape_adj <- function(lower, upper, trunclower=0, truncupper=Inf, theta, shape, beta, eps=1e-03, print=TRUE){
  shape <- shape[beta>0]
  beta <- beta[beta>0]/sum(beta)
  M <- length(shape)
  fit <- ME_em(lower, upper, trunclower, truncupper, theta, shape, beta, eps, print=FALSE)
  loglikelihood <- fit$loglikelihood
  theta <- fit$theta
  beta <- fit$beta
  alpha <- fit$alpha
  # before and after are the loglikelihoods used in the outer while loop
  before_loglikelihood <- -Inf
  after_loglikelihood <- loglikelihood    
  iteration <- 1
  while(after_loglikelihood > before_loglikelihood + eps){    
    if(print) cat("iteration = ", iteration, "\n")
    before_loglikelihood <- after_loglikelihood
    # Try increasing the shapes
    for(i in M:1){
      improve <- TRUE
      while( (improve==TRUE) && (i == M || (shape[i] < shape[i+1]-1))) {
        new_shape <- shape
        new_shape[i] <- new_shape[i]+1        
        fit <- ME_em(lower, upper, trunclower, truncupper, theta, new_shape, beta, eps, print=FALSE)
        new_loglikelihood <- fit$loglikelihood
        if(new_loglikelihood > loglikelihood + eps){
          loglikelihood <- new_loglikelihood
          shape <- new_shape
          theta <- fit$theta
          beta <- fit$beta 
          alpha <- fit$alpha
          if(print) cat("loglikelihood = ", loglikelihood, ", shape = ", shape, "\n", "theta = ", theta, ", alpha = ", alpha, "\n")      
        } else{improve <- FALSE}
      }
    }
    # Try decreasing the shapes
    for(i in 1:M){
      improve <- TRUE
      while( (improve==TRUE) && ( (i == 1) || shape[i] > shape[i-1]+1 ) && shape[i]>1){
        new_shape <- shape
        new_shape[i] <- new_shape[i]-1
        fit <- ME_em(lower, upper, trunclower, truncupper, theta, new_shape, beta, eps, print=FALSE)
        new_loglikelihood <- fit$loglikelihood
        if(new_loglikelihood > loglikelihood + eps){
          loglikelihood <- new_loglikelihood
          shape <- new_shape
          theta <- fit$theta
          beta <- fit$beta 
          alpha <- fit$alpha
          if(print) cat("loglikelihood = ", loglikelihood, ", shape = ", shape, "\n", "theta = ", theta, ", alpha = ", alpha, "\n")      
        } else{improve <- FALSE}
      }          
    }
    after_loglikelihood <- loglikelihood  
    iteration <- iteration + 1
  }
  list(alpha = alpha, beta = beta, shape = shape, theta = theta, loglikelihood = loglikelihood, AIC=-2*loglikelihood+2*(2*length(alpha)),BIC=-2*loglikelihood+(2*length(alpha))*log(length(data))) 
}

## Reduction of M based on an information criterium: AIC and BIC implemented

ME_shape_red <- function(lower, upper, trunclower=0, truncupper=Inf, theta, shape, beta, criterium="AIC", eps=1e-03, print=TRUE){
  n <- length(lower)
  fit <- ME_shape_adj(lower, upper, trunclower, truncupper, theta, shape, beta, eps, print=FALSE)
  loglikelihood <- fit$loglikelihood
  IC <- fit[[criterium]]
  shape <- fit$shape
  theta <- fit$theta
  beta <- fit$beta   
  alpha <- fit$alpha
  M <- length(shape)
  if(print) cat("M = ", M, ", ", criterium, " = ", IC, ", shape = ", shape, "\n", "theta = ", theta, ", alpha = ", alpha, "\n")
  improve <- TRUE
  while((improve==TRUE) && length(shape) > 1){    
    new_shape <- shape[beta != min(beta)]
    new_beta <- beta[beta != min(beta)]
    new_beta <- new_beta/sum(new_beta)
    fit <- ME_shape_adj(lower, upper, trunclower, truncupper, theta, new_shape, new_beta, eps, print=FALSE)
    new_IC <- fit[[criterium]]
    if(new_IC < IC){ 
      IC <- new_IC
      loglikelihood <- fit$loglikelihood  
      shape <- fit$shape
      theta <- fit$theta
      beta <- fit$beta  
      alpha <- fit$alpha
      M <- length(shape)
      if(print) cat("M = ", M, ", ", criterium, " = ", IC, ", shape = ", shape, "\n", "theta = ", theta, ", alpha = ", alpha, "\n")
    } else{improve <- FALSE}        
  }
  list(M = M, alpha = alpha, beta = beta, shape = shape, theta = theta, loglikelihood = loglikelihood, AIC=-2*loglikelihood+2*(2*length(alpha)),BIC=-2*loglikelihood+(2*length(alpha))*log(n)) 
}

## Calibration procedure for mixtures of Erlangs by repeatedly using the EM algorithm while adjusting and reducing the shape parameters based on an information criterium (AIC and BIC implemented)

# Specify lower and upper censoring points (lower, upper), lower and upper truncation points (trunclower, truncupper).
# The censoring status is determined as follows:
# Uncensored: lower and upper are both present and equal.
# Left Censored: lower is missing (NA), but upper is present.
# Right Censored: lower is present, but upper is missing (NA).
# Interval Censored: lower and upper are present and different.
# e.g.: lower=c(1,NA,3,4); upper=c(1,2,NA,5); specifies an observed event at 1, left censoring at 2, right censoring at 3, and interval censoring at [4,5],
# By default no truncation: trunclower=0, truncupper=Inf
# alpha = beta in case of no truncation

ME_fit <- function(lower, upper = lower, trunclower = 0, truncupper = Inf, M = 10, s = 1, criterium="AIC", eps=1e-03, print=TRUE){
  initial <- ME_initial(lower, upper, trunclower, truncupper, M, s)
  fit <- ME_shape_red(lower, upper, trunclower, truncupper, initial$theta, initial$shape, initial$beta, criterium, eps, print)
  list(alpha = fit$alpha, beta = fit$beta, shape = fit$shape, theta = fit$theta, loglikelihood = fit$loglikelihood, AIC=fit$AIC, BIC=fit$BIC, M = fit$M, M_initial = M, s = s) 
}

## Tune the initialising parameters M and s using a grid search over the supplied parameter ranges

ME_tune <- function(lower, upper = lower, trunclower = 0, truncupper = Inf, M = 10, s = 1, nCores = detectCores(), criterium = "AIC", eps = 1e-03, print=TRUE, file="log.txt"){
  tuning_parameters = expand.grid(M, s)
  cl <- makePSOCKcluster(nCores)
  registerDoParallel(cl)  
  if(print) writeLines(c(""), file)
  all_model <- foreach(i = 1:nrow(tuning_parameters), .export=c("ME_initial", "ME_loglikelihood", "ME_u_z", "ME_c_z", "ME_expected_c", "ME_T", "theta_nlm_u_c", "theta_nlm_u", "theta_nlm_c", "ME_em", "ME_shape_adj", "ME_shape_red", "ME_fit"), .errorhandling = 'remove') %dopar% {
    if(print) cat(paste("M = ", tuning_parameters[i, 1], ", s = ", tuning_parameters[i, 2], "\n"), file = file, append = TRUE)
    ME_fit(lower, upper, trunclower, truncupper, M = tuning_parameters[i, 1], s = tuning_parameters[i, 2], criterium, eps, FALSE)
  }
  stopCluster(cl)
  performances <- data.frame(tuning_parameters[,1], tuning_parameters[,2], sapply(all_model, function(x) with(x, get(criterium))), sapply(all_model, with, M))
  colnames(performances) = c('M_initial', 's', criterium, 'M')
  best_index <- which.min(performances[, criterium])
  best_model <- all_model[[best_index]]  
  list(best_model = best_model, performances = performances, all_model = all_model)
}

## Density function

ME_density <- function(x, theta, shape, alpha, trunclower = 0, truncupper = Inf, log = FALSE){
  f <- outer(x, shape, dgamma, scale = theta)
  d <- rowSums(t(t(f)*alpha))
  if(!(trunclower==0 & truncupper==Inf)){
    d <- d / (ME_cdf(truncupper, theta, shape, alpha) - ME_cdf(trunclower, theta, shape, alpha)) * ((trunclower <= x) & (x <= truncupper))
  }
  if(log){
    d <- log(d)
  }
  d
}

## Cumulative distribution function

ME_cdf <- function(x, theta, shape, alpha, trunclower = 0, truncupper = Inf, lower.tail = TRUE, log.p = FALSE){     
  cdf <- outer(x, shape, pgamma, scale=theta)
  p <- rowSums(t(t(cdf)*alpha))
  if(!(trunclower==0 & truncupper==Inf)){
    l <- ME_cdf(trunclower, theta, shape, alpha)
    u <- ME_cdf(truncupper, theta, shape, alpha)
    p <- ((p - l) / (u - l)) ^ {(x <= truncupper)} * (trunclower <= x)
  }
  if(!lower.tail){
    p <- 1 - p
  }
  if(log.p){
    p <- log(p)
  }
  p
} 

## Noncentral moments of order k (k can be vector)

ME_moment <- function(k, theta, shape, alpha){
  alpha_factorials <- t(exp(t(lgamma(outer(k,shape,'+'))) - lgamma(shape)) * alpha)
  rowSums(alpha_factorials)*theta^k
}

## Value-at-Risk (VaR) or quantile function

ME_VaR <- function(p, theta, shape, alpha, trunclower = 0, truncupper = Inf, interval = if(trunclower == 0 & truncupper == Inf){c(qgamma(p, shape = min(shape), scale = theta), qgamma(p, shape = max(shape), scale = theta))}else{c(trunclower, min(truncupper, trunclower + qgamma(p, shape = max(shape), scale = theta)))}, start = qgamma(p, shape = shape[which.max(alpha)], scale = theta)){
  if(p==1){
   return(Inf) 
  }    
  if(length(shape) == 1 & trunclower == 0 & truncupper == Inf){
    VaR <- qgamma(p, shape = shape, scale = theta)
  }
  else{
    objective <- function(x){return(10000000*(ME_cdf(x, theta, shape, alpha, trunclower, truncupper)-p)^2)}    
    VaR_nlm <- nlm(f = objective, p = start)
    VaR_optimize <- optimize(f = objective, interval = interval)
    VaR <- ifelse(VaR_nlm$minimum < VaR_optimize$objective, VaR_nlm$estimate, VaR_optimize$minimum)    
  }
  if(objective(VaR)>1e-06){ # in case optimization fails, retry with more different starting values
    alpha <- alpha[order(shape)]
    shape <- shape[order(shape)]
    VaR_nlm <-  vector("list", length(shape))
    VaR_optimize <-  vector("list", length(shape))
    interval <- c(0, qgamma(p, shape, scale = theta))
    for(i in 1:length(shape)){
      VaR_nlm[[i]] <- nlm(f = objective, p = qgamma(p, shape = shape[i], scale = theta))    
      VaR_optimize[[i]] <- optimize(f = objective, interval = interval[c(i, i+1)])
    }
    VaR_nlm <- sapply(VaR_nlm, with, estimate)[which.min(sapply(VaR_nlm, with, minimum))]
    VaR_optimize <- sapply(VaR_optimize, with, minimum)[which.min(sapply(VaR_optimize, with, objective))]
    VaR <- ifelse(objective(VaR_nlm) < objective(VaR_optimize), VaR_nlm, VaR_optimize)  
  }  
  VaR  
}

ME_VaR <- Vectorize(ME_VaR, vectorize.args = c("p", "start"))

## Tail-Value-at-Risk (TVaR)

# using a Newton-type algorithm

ME_TVaR <- function(p, theta, shape, alpha, interval = if(trunclower == 0 & truncupper == Inf){c(qgamma(p, shape = min(shape), scale = theta), qgamma(p, shape = max(shape), scale = theta))}else{c(trunclower, min(truncupper, trunclower + qgamma(p, shape = max(shape), scale = theta)))}, start = qgamma(p, shape = shape[which.max(alpha)], scale = theta)){     
  VaR <- ME_VaR(p, theta, shape, alpha, interval = interval, start = start)
  M <- max(shape)
  shapes <- seq(1, M)
  alphas <- rep(0, M)
  alphas[shape] <- alpha
  A <- rev(cumsum(rev(alphas)))
  # note: A_n = A[n+1]
  AA <- rev(cumsum(rev(A)))
  # note: AA_n = AA[n+1]
  VaR + theta^2/(1-p)*sum(AA*dgamma(VaR, shapes, scale=theta))  
}

ME_TVaR <- Vectorize(ME_TVaR, vectorize.args = c("p", "start"))

## Expected Shortfall (ESF)

ME_ESF <- function(R, theta, shape, alpha){ 
  M <- max(shape)
  shapes <- seq(1, M)
  alphas <- rep(0, M)
  alphas[shape] <- alpha
  A <- rev(cumsum(rev(alphas)))
  # note: A_n = A[n+1]
  AA <- rev(cumsum(rev(A)))
  # note: AA_n = AA[n+1]
  theta^2*sum(AA*dgamma(R, shapes, scale=theta))  
}

## Stop loss moments of order delta (delta doesn't have to be an integer), given lower and upper truncation bounds
## Equal to ESF for delta = 1, trunclower = 0, truncupper = Inf

ME_SLM <- function(R, delta = 1, theta, shape, alpha, trunclower = 0, truncupper = Inf){ 
  M <- max(shape)
  shapes <- seq(1, M)
  alphas <- rep(0, M)
  alphas[shape] <- alpha
  coeff <- rep(0, M)
  for(n in 1:M){
    coeff[n] <- sum( alphas[0:(M-n)+n] * exp(lgamma(delta+0:(M-n)+1) - lgamma(0:(M-n)+1)) * (pgamma(truncupper-R, shape = delta+0:(M-n)+1, scale = theta) - pgamma(max(trunclower, R)-R, shape = delta+0:(M-n)+1, scale = theta)) )
  }  
  theta^(delta+1) / (ME_cdf(truncupper, theta, shape, alpha) - ME_cdf(trunclower, theta, shape, alpha)) * sum( coeff * dgamma(R, shapes, scale=theta) )  
}

## Excess-of-loss reinsurance premium: L xs R

## Excess-of-loss reinsurance premium: C xs R (Retention R, Cover C, Limit L = R+C)

ME_XL <- function(R, C, theta, shape, alpha){ 
  M <- max(shape)
  shapes <- seq(1, M)
  alphas <- rep(0, M)
  alphas[shape] <- alpha
  coeff <- rep(0, M)
  for(n in 1:M){
    coeff[n] <- sum( alphas[0:(M-n)+n] * (0:(M-n)+1) * pgamma(C, shape = 0:(M-n)+2, scale = theta) )
  }
  if(C == Inf){
    XL <- theta^2 * sum( coeff * dgamma(R, shapes, scale=theta) )
  }else{
    XL <- theta^2 * sum( coeff * dgamma(R, shapes, scale=theta) ) +  C * mix.erlang.cdf(R+C, theta, shape, alpha, lower.tail = FALSE)    
  }
  XL  
}

## Excess loss or residual lifetime distribution

ME_excess_loss <- function (R, theta, shape, alpha) {
  M <- max(shape)
  alphas <- rep(0, M)
  alphas[shape] <- alpha
  # note: A_n = A[n+1]
  A <- rev(cumsum(rev(alphas)))
  shape_el <- 1:M
  alpha_el <- rep(0, M)
  f <- dgamma(R, shape_el, scale = theta)
  denom <- sum(A*f)  
  for(i in 1:M){
    alpha_el[i] <- sum(alphas[i:M]*f[1:(M-i+1)]) / denom
  }
  list(theta = theta, shape = shape_el, alpha = alpha_el)
}



## Random generation

ME_random <- function(n, theta, shape, alpha){   
  rgamma(n, shape=sample(shape, size=n, replace=TRUE, prob=alpha), scale=theta)  
} 

## Plot density of mixture of erlangs

ME_plot_density <- function(theta, shape, alpha, trunclower = 0, truncupper = Inf, xlim = c(0, 10), ...) {
  x <-  seq(xlim[1], xlim[2], by=(xlim[2]-xlim[1])/10000)
  y <- ME_density(x, theta, shape, alpha, trunclower, truncupper)  
  plot(x, y, type="l", xlim = xlim, ...)
}

## Plot cdf of mixture of erlangs

ME_plot_cdf <- function(theta, shape, alpha, trunclower = 0, truncupper = Inf, lower.tail = TRUE, xlim = c(0, 10), ...) {
  x <-  seq(xlim[1], xlim[2], by=(xlim[2]-xlim[1])/10000)
  y <- ME_cdf(x, theta, shape, alpha, trunclower, truncupper, lower.tail)  
  plot(x, y, type="l", xlim = xlim, ...)
}

## Plot density of mixture of erlangs and histogram of data

ME_plot_data <- function(data, shape, alpha, theta, trunclower = 0, truncupper = Inf, nbins = 50, xlim = c(max(c(min(data)-0.1*(max(data)-min(data)),0)), max(data)+0.1*(max(data)-min(data))), xlab = "", legend = TRUE, lwd = 2, ...) {
  x <-  seq(xlim[1], xlim[2], by=(xlim[2]-xlim[1])/10000)
  y <- ME_density(x, theta, shape, alpha, trunclower, truncupper)  
  truehist(data, h = (xlim[2]-xlim[1])/nbins, x0 = trunclower, xlim = xlim, xlab = xlab, ... )
  lines(x, y, col = "red", lwd = lwd)  
  if(legend){
    legend('topright', legend = if(trunclower==0 & truncupper==Inf) c("Fitted Density Function", "Observed Relative Frequency") else c("Fitted Truncated Density Function", "Observed Relative Frequency"), col = c("red","cyan"), pch=c(NA,15), pt.cex=2, lty = c(19,NA), lwd=c(lwd,NA))
  }
}  

## Plot density of mixture of erlangs, histogram of simulated data and true density

ME_plot_sim_data <- function(data, dens, dens_param, shape, alpha, theta, trunclower = 0, truncupper = Inf, nbins = 50, xlim = c(max(c(min(data)-0.1*(max(data)-min(data)),0)), max(data)+0.1*(max(data)-min(data))), xlab = "", legend = TRUE, lwd = 2, ...) {
  x <-  seq(xlim[1], xlim[2], by=(xlim[2]-xlim[1])/10000)
  y <- ME_density(x, theta, shape, alpha, trunclower, truncupper)  
  yy <- do.call(dens, c(list(x = x), dens_param))
  truehist(data, h = (xlim[2]-xlim[1])/nbins, x0 = trunclower, xlim = xlim, xlab = xlab, ... )
  lines(x, y, col = "red", lwd = lwd)
  lines(x, yy, col = "blue", lwd = lwd, lty=2)
  if(legend){
    legend('topright', legend = if(trunclower==0 & truncupper==Inf) c("Fitted Density Function", "True Density Function", "Observed Relative Frequency") else c("Fitted Truncated Density Function", "True Density Function", "Observed Relative Frequency"), col = c("red", "blue", "cyan"), pch=c(NA, NA,15), pt.cex=2, lty = c(1,2,NA), lwd=c(lwd, lwd,NA))
  }
}