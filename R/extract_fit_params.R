
#' Returns the fitted parameters for a model 1 or model 2 object.
#'
#' Note that in all cases these return the median value from the post-warmup draws.
#'
#' @export
#' @param fit fitted rstan model from model1 or model2 fit
#' @param ndim optional parameter to specify number of categories (e.g. male, female), defaults to 2
#' @return named list of fit parameters
extract_fit_params <- function(fit, ndim=2){
  .my_assert("error. model is not of class 1 or 2", fit@model_name %in% c("model1", "model2"))

  p <- get_proportions(fit)

  if (fit@model_name == "model1"){
    Sigma <- get_var_covar_matrix(fit, ndim)
    rg <- get_gen_cor(fit, ndim)
    return(list("prop"=p, "var_cov"=Sigma, "gen_cor"=rg))
  }
  if (fit@model_name == "model2"){
    sigmasq <- get_vars(fit)
    return(list("prop"=p, "vars"=sigmasq))
  }
}

#' Get the proportions associated with a fit
#'
#' @export
#' @param fit fitted rstan model from model1 or model2 fit
#' @return vector containing the proportions for the fit
get_proportions <- function(fit){
  fit_summ_pi <- rstan::summary(fit, pars=c("pi"), probs=c(0.05, 0.50, 0.95))
  p <-fit_summ_pi$summary[,c("50%")]
  return(p)
}

#' Get the variance-covariance matrix associated with a model1 fit
#'
#' @export
#' @param fit1 fitted rstan model from model1 fit
#' @param ndim optional parameter to specify number of categories (e.g. male, female), defaults to 2
#' @return variance covariance matrix of size ndim x ndim
get_var_covar_matrix <- function(fit1, ndim=2){
  .my_assert("error. please provide a model1 fit", fit1@model_name=="model1")
  fit_summ_S <- rstan::summary(fit1, pars=c("Sigma"), probs=c(0.05, 0.50, 0.95))
  Sigma <- matrix(fit_summ_S$summary[,c("50%")], ndim, ndim)
  return(Sigma)
}

#' Get the genetic correlation associated with a model1 fit
#'
#' @export
#' @param fit1 fitted rstan model from model1 fit
#' @param ndim optional parameter to specify number of categories (e.g. male, female), defaults to 2
#' @return list of pairwise genetic correlations ndim-1 long (for ndim==2 it is a single value)
get_gen_cor <- function(fit1, ndim=2){
  .my_assert("error. please provide a model1 fit", fit1@model_name=="model1")

  fit_summ_R <- rstan::summary(fit1, pars=c("Omegacor"), probs=c(0.05, 0.50, 0.95))
  rg <- matrix(fit_summ_R$summary[,c("50%")], ndim, ndim)
  return(rg[upper.tri(rg, diag=FALSE)])
}

#' Get a 95% confidence interval for the genetic correlation
#'
#' @export
#' @param fit1 fitted rstan model from model1 fit
#' @param ndim optional parameter to specify number of categories (e.g. male, female), defaults to 2
#' @return lists of lower (l) and upper (u) CIs for the each of the pairwise genetic correlations
#'        each sublist will be ndim-1 long (for ndim==2 it is a single value), ordered by the beta pairs
get_gen_cor_ci <- function(fit1, ndim=2){
  .my_assert("error. please provide a model1 fit", fit1@model_name=="model1")

  fit_summ_R <- rstan::summary(fit1, pars=c("Omegacor"), probs=c(0.025, 0.975))
  pairs <- combn(ndim, 2)
  list.labels <- apply(pairs, 2, function(x) sprintf("Omegacor[%s,%s]", x[1], x[2]))
  conf.l <- fit_summ_R$summary[list.labels,"2.5%"]
  conf.u <- fit_summ_R$summary[list.labels,"97.5%"]
  return(list("l"=conf.l, "u"=conf.u))
}



#' Get the variances from the m2 fit
#'
#' @param fit2 the output of an m2
#' @return a list of the variances estimated from the m2 fit
get_vars <- function(fit2){
  .my_assert("error. please provide a model2 fit", fit2@model_name=="model2")

  fitS <- rstan::summary(fit2, pars=c("sigmasq"), probs=c(0.05, 0.50, 0.95))
  sigmasq <- as.vector(fitS$summary[,c("50%")])
  return(sigmasq)
}




