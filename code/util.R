### adapt the ebnm function for prior (1-pi)*delta_0 + pi*N+(mu, sigma^2), where N+ represents truncated normal on [0, \infty)
ebnm_npmle_truncnorm <- function(x, s, g_init, output=ebnm::output_default()){
  s2 <- s^2
  pi <- g_init$pi[2]
  mu <- g_init$mean[2]
  sigma2 <- g_init$sd[2]^2
  
  ############################## calculate the posterior which is a mixture of delta_0 and truncated normal ##############################
  # calculate mean and variance of posterior corresponding to the truncated normal component in the mixture prior
  wgts <- s2/(s2 + sigma2)
  post_mean <- wgts*mu + (1-wgts)*x
  post_var <- 1/(1/sigma2 + 1/s2)
  post_sd <- pmax(sqrt(post_var), 1e-15)
  
  # calculate the marginal log likelihood for each mixture component
  llik_mat <- matrix(NA, nrow=length(x), ncol=2)
  llik_mat[,1] <- -0.5*(log(s2) + x^2/s2)
  llik_mat[,2] <- -0.5*(log(sigma2+s2) + (x-mu)^2/(sigma2+s2)) + 
    pnorm(-post_mean/post_sd, lower.tail = FALSE, log.p = TRUE) - pnorm(-mu/sqrt(sigma2), lower.tail = FALSE, log.p = TRUE)
  llik_norms <- apply(llik_mat, 1, max)
  L_mat <- exp(llik_mat - llik_norms)
  
  # calculate the posterior weight for truncated normal in the posterior
  zeta_mat <- t(t(L_mat)*g_init$pi)
  zeta_mat <- zeta_mat*(1/rowSums(zeta_mat)) 
  zeta <- zeta_mat[,2]
  
  # create the list for return
  ebnm_res <- list()
  
  # return the prior mixture
  if("fitted_g" %in% output)
    ebnm_res$fitted_g <- g_init
  
  # return the posterior probabilities zeta
  if("posterior_prob" %in% output)
    ebnm_res$zeta <- zeta
  
  # calculate the log marginal likelihood 
  if ("log_likelihood" %in% output) 
    ebnm_res$log_likelihood <- sum(log(L_mat %*% g_init$pi)) + sum(llik_norms) - 0.5*length(x)*log(2*3.141593)
  
  # calculate posterior of theta 
  if ("posterior_mean" %in% output) {
    tmp1 <- truncnorm::etruncnorm(a=0, mean=post_mean, sd=post_sd)
    tmp1[is.nan(tmp1)] <- 0
    tmp1[tmp1 < 0] <- 0
    
    tmp2 <- truncnorm::vtruncnorm(a=0, mean=post_mean, sd=post_sd)
    tmp2[is.nan(tmp2)] <- 0
    tmp2[tmp2 < 0] <- 0
    
    ebnm_res$posterior$mean <- zeta*tmp1
    ebnm_res$posterior$second_moment <- zeta*(tmp1^2 + tmp2)
    ebnm_res$posterior$sd <- sqrt(pmax(0, ebnm_res$posterior$second_moment - ebnm_res$posterior$mean^2))
    ebnm_res$posterior$lfsr <- 1-zeta
  }
  
  return(ebnm_res)
}


### define the generalized binary prior (1-pi)*delta_0 + pi*N+(mu, sigma^2), where N+ represents truncated normal on [0, \infty), mu/sigma = 10
### split the optimization problem into several subproblems with different range of pi to optimize over
ebnm_binary_general <- function(x, s, g_init, fix_g, output) {
  if (!fix_g) {
    s2 <- s^2
    
    ### separately optimize over parameters in different intervals for pi and take the best solution
    wlist <- c(1e-5, 1e-2, 1e-1, 2.5e-1, 9e-1)
    opt_list <- list(NULL)
    val_list <- rep(NA, length(wlist)-1)
    for(k in 1:(length(wlist)-1)){
      ### initialize g
      w <- wlist[k+1]/2
      mu <- 1
      
      ### update g
      for(iter in 1:50){
        g_init <- ashr::normalmix(c(1-w, w), c(0, mu), c(0, mu/10))
        zeta <- ebnm_npmle_truncnorm(x=x, s=s, g_init=g_init, output="posterior_prob")$zeta
        
        # update g given the expectation of latent variables z
        w_new <- pmin(pmax(mean(zeta), wlist[k]), wlist[k+1])
        
        # define the objective function of mu to maximize
        opt_fn <- function(par){
          par_sigma2 <- (par/10)^2
          # calculate mean and variance of posterior corresponding to the truncated normal component in the mixture prior
          wgts <- s2/(s2 + par_sigma2)
          post_mean <- wgts*par + (1-wgts)*x
          post_var <- 1/(1/par_sigma2 + 1/s2)
          post_sd <- pmax(sqrt(post_var), 1e-15)
          llik <- -0.5*(log(par_sigma2+s2) + (x-par)^2/(par_sigma2+s2)) + 
            pnorm(-post_mean/post_sd, lower.tail = FALSE, log.p = TRUE) - pnorm(-par/sqrt(par_sigma2), lower.tail = FALSE, log.p = TRUE)
          # calculate objective function
          return(sum(zeta*llik))
        }
        
        # update mu
        mu_new <- optim(par=1, fn=opt_fn, lower=1e-3, upper=32, method="L-BFGS-B", control=list(fnscale=-1))$par
        
        # check for stopping criterion
        if(abs(w_new-w) < 1e-3 & abs(mu_new - mu) < 1e-3*mu)
          break
        mu <- mu_new
        w <- w_new
      }
      
      ### return the estimated g and likelihood  
      g_init <- ashr::normalmix(c(1-w, w), c(0, mu), c(0, mu/10))
      opt_list[[k]] <- g_init
      val_list[k] <- ebnm_npmle_truncnorm(x, s, g_init = g_init, output = "log_likelihood")$log_likelihood
    }
    
    ### take the optimal g among all intervals
    g_init <- opt_list[[which.max(val_list)]]
  }
  
  return(ebnm_npmle_truncnorm(x, s, g_init = g_init, output = output))
}


### initialize the nonnegative fit to covariance matrix XX' s.t. E[XX'] = LL'+ D based on unconstrained estimate of L
init.snn.cov <- function(fl, kset=1:ncol(fl$flash.fit$EF[[1]])) {
  LL <- fl$flash.fit$EF[[1]][, kset, drop = FALSE]
  FF <- fl$flash.fit$EF[[2]][, kset, drop = FALSE]
  
  LL <- cbind(pmax(LL, 0), pmax(-LL, 0))
  LL <- cbind(fl$flash.fit$EF[[1]][, -kset, drop = FALSE], LL)
  FF <- cbind(pmax(FF, 0), pmax(-FF, 0))
  FF <- cbind(fl$flash.fit$EF[[2]][, -kset, drop = FALSE], FF)
  
  to.keep <- (colSums(LL) > .Machine$double.eps) & (colSums(FF) > .Machine$double.eps)
  LL <- LL[, to.keep, drop = FALSE]
  FF <- FF[, to.keep, drop = FALSE]
  
  return(list(LL, FF))
}


### apply flash to covariance matrix XX' s.t. E[XX'] = LL'+ D, where D = sigma2*I with pre-determined initialization
fit.ebcovmf <- function(dat, fl, prior, extrapolate=TRUE, maxiter=500, epsilon=0, verbose=1){
  s2 <- max(0, mean(diag(dat) - diag(fitted(fl))))
  s2_diff <- Inf
  
  # Alternate between estimating s2 and backfitting until convergence.
  while(s2 > 0 && abs(s2_diff - 1) > 1e-3) {
    dat_minuss2 <- dat - diag(rep(s2, ncol(dat)))
    kset <- which(fl$pve > epsilon)
    kmax <- which.max(fl$pve)
    kset <- kset[-c(kmax)]
    fl <- flash.init(dat_minuss2) %>% flash.init.factors(init = lapply(fl$flash.fit$EF, function(x) x[, kmax, drop = FALSE]),
                     ebnm.fn = ebnm::ebnm_unimodal_nonnegative) %>% flash.init.factors(
                     init = lapply(fl$flash.fit$EF, function(x) x[, kset, drop = FALSE]), ebnm.fn = prior) %>% flash.backfit(
                     extrapolate = extrapolate, maxiter = maxiter, verbose = verbose)
    old_s2 <- s2
    s2 <- max(0, mean(diag(dat) - diag(fitted(fl))))
    s2_diff <- s2 / old_s2
  }
  
  return(list(dat=dat, fl=fl, s2=s2))
}


### run flash to the covariance matrix by adding all factors greedily, then backfitting all of them (initialization using point laplace)
### the default is to use the generalized binary prior, but the function can also be called to apply other priors defined in ebnm package
fit.gbcd <- function(dat, Kmax, prior=ebnm_binary_general, extrapolate=TRUE, maxiter=500, verbose=1){
  ### fit unconstrained flash with point laplace prior to XXt without considering the diagonal component for now
  fit.cov <- flash.init(dat, var.type = 0) %>% flash.add.greedy(Kmax = 1, ebnm.fn = ebnm::ebnm_unimodal_nonnegative, 
                        init.fn = function(f) init.fn.default(f, dim.signs = c(1, 1)), verbose = 1) %>% flash.add.greedy(
                        Kmax = Kmax-1, ebnm.fn = ebnm::ebnm_point_laplace, verbose = 1) %>% flash.backfit(maxiter = maxiter, verbose = 1)
  
  ### fit unconstrained flash again with the diagonal component
  fit.cov <- fit.ebcovmf(dat = dat, fl = fit.cov, prior = ebnm::ebnm_point_laplace, maxiter = maxiter, verbose = 1)$fl
  
  ### initialize the nonnegative fit based on the unconstrained flash fit 
  snn.cov.init <- init.snn.cov(fit.cov)
  kmax <- which.max(colSums(snn.cov.init[[1]]))
  fit.cov <- flash.init(dat, var.type = 0) %>% flash.init.factors(init = lapply(snn.cov.init, function(x) x[, kmax, drop = FALSE]), 
                        ebnm.fn = ebnm::ebnm_unimodal_nonnegative) %>% flash.init.factors(
                        init = lapply(snn.cov.init, function(x) x[, -c(kmax), drop = FALSE]), ebnm.fn = prior
                        ) %>% flash.backfit(extrapolate = extrapolate, maxiter = 25, verbose = verbose) 
  
  ### fit flash using generalized binary prior with the diagonal component
  fit.cov <- fit.ebcovmf(dat = dat, fl = fit.cov, prior = prior, extrapolate = extrapolate, maxiter = 25, verbose = verbose)$fl
  
  ### select a subset of factors and refit
  kset <- (length(fit.cov$pve) - rank(fit.cov$pve) < Kmax) & (fit.cov$pve > 0)
  kall <- 1:fit.cov$n.factors
  if(!all(kset))
    fit.cov <- flash.remove.factors(fit.cov, kset=kall[!kset])
  fit.cov <- fit.ebcovmf(dat = dat, fl = fit.cov, prior = prior, extrapolate = extrapolate, maxiter = maxiter, verbose = verbose)$fl
  
  return(fit.cov)
}


### define the function to plot estimated L from GBCD fit
plot.gbcd <- function(fit.cov, rho=0.8, title=""){
  ### obtain the consistent factors from covariance decomposition
  k.order <- order(fit.cov$pve, decreasing = TRUE)
  fit.L <- fit.cov$L.pm[, k.order]
  fit.F <- fit.cov$F.pm[, k.order]
  corr <- diag(cor(fit.L, fit.F))
  fit.L <- t(t(fit.L)/apply(fit.L, 2, max))
  fit.L <- fit.L[, corr > rho]
  colnames(fit.L) <- paste0("k", 1:ncol(fit.L))
  
  ### define the color map
  cols <- colorRampPalette(c("gray96", "red"))(99)
  brks <- seq(0, 1, length=100)
  
  ### plot the heatmap for L
  plt <- pheatmap(fit.L, show_rownames = FALSE, cluster_rows = FALSE, cluster_cols = FALSE, color = cols, breaks = brks, main = title)
  return(plt)
}