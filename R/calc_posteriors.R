#' Calculate posteriors for the variants.
#'
#' @export
#' @param fit the output of model1 or model2
#' @param dat a list with B and SE elements, variant IDs; from prep_gwas_input()
#' @return a dataframe with the posteriors for each of the variants:
#'         for m1 this contains the variant name and only the non-null posterior,
#'         for m2 this contains the variant name and all four posteriors
calc_posteriors <- function(fit, dat){
  .my_assert("error. model is not of class 1 or 2", fit@model_name %in% c("model1", "model2"))
  .check_in_dat_format(dat)

  # grab the data and fit output
  B_dat <- dat$B
  SE_dat <- dat$SE
  N <- nrow(B_dat)
  p <- get_proportions(fit)

  if (fit@model_name=="model1"){
    Sigma <- get_var_covar_matrix(fit)
    posteriors <- sapply(1:N, function(i)
      .calc_posterior_variant_m1(B_dat[i,], SE_dat[i,], p, Sigma))
    posterior.df <- data.frame(posteriors)
    colnames(posterior.df) <- c("p2")
  }

  if (fit@model_name=="model2"){
    sigmasq <- get_vars(fit)
    posteriors <- lapply(1:N, function(i)
      .calc_posterior_variant_m2(B.dat[i,], SE.dat[i,], p, sigmasq, Sigma))
    posterior.df <- data.frame(do.call(rbind, posteriors))
    colnames(posterior.df) <- c("p1", "p2", "p3", "p4")
  }
  posterior.df$ID <- dat$id
  return(posterior.df)
}


#' Calculate the posterior probability for an m1 variant
#'
#' @param B betas from the dat object
#' @param SE ses from the dat object
#' @param p proportion vector from the m1 fit
#' @param Sigma var-covar matrix from the m1 fit
#' @return the posterior probability of belonging to the non-null component
.calc_posterior_variant_m1 <- function(B, SE, p, Sigma){
  zeros <- rep(0, length(SE))
  SE_mat <- diag(SE)
  p_1 = p[1]*mnormt::dmnorm(B, zeros, SE_mat)
  p_2 = p[2]*mnormt::dmnorm(B, zeros, SE_mat + Sigma)
  prob_1 = log(p_1) - log(p_1 + p_2)
  prob_2 = log(p_2) - log(p_1 + p_2)
  return(exp(prob_2))
}

#' Calculate the posterior probability for an m2 variant
#'
#' @param B betas from the dat object
#' @param SE ses from the dat object
#' @param p proportion vector from the m2 fit
#' @param sigmasq variances from the m2 fit
#' @return a vector with the posterior probability of belonging to each component
.calc_posterior_variant_m2 <- function(B, SE, p, sigmasq){
  zeros <- c(0,0)
  SE_mat <- matrix(c(SE[1], 0, 0, SE[2]), 2, 2)
  Sigma <- matrix(c(sigmasq[3], 0, 0, sigmasq[4]),2,2)
  p_1 = p[1]*mnormt::dmnorm(B, zeros, SE_mat)
  p_2 = p[2]*mnormt::dmnorm(B, zeros, SE_mat + matrix(c(sigmasq[1], 0, 0, 0),2, 2))
  p_3 = p[3]*mnormt::dmnorm(B, zeros, SE_mat + matrix(c(0, 0, 0, sigmasq[2]),2,2))
  p_4 = p[4]*mnormt::dmnorm(B, zeros, SE_mat + Sigma)
  p_tot = p_1 + p_2+ p_3 + p_4
  prob_1 = exp(log(p_1) - log(p_tot))
  prob_2 = exp(log(p_2) - log(p_tot))
  prob_3 = exp(log(p_3) - log(p_tot))
  prob_4 = exp(log(p_4) - log(p_tot))
  return(list(prob_1, prob_2, prob_3, prob_4))
}







