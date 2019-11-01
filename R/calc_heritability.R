#' Calculate the heritability for a trait.
#' Note: current code does not return the confidence interval;
#' this will be added in the next version.
#'
#' @export
#' @param fit1 a model1 fit
#' @param dat a list with B and SE elements, variant IDs; from prep_gwas_input()
#' @param posterior.df optional: provide the posterior data.frame if you already
#'        have it to avoid recalculating
#' @param cutoff the cutoff to assign to the non-null component, defaults to 0.8
#' @return a list of the heritabilities
calc_heritability <- function(fit1, dat, posterior.df=NULL, cutoff=0.8){
  .my_assert("error. model is not of class 1", fit1@model_name=="model1")
  .check_in_dat_format(dat)

  # calculate posterior probability and assign to categories
  if (is.null(posterior.df)){
    posterior.df <- calc_posteriors(fit1, dat)
  }
  var.categories <- assign_to_components(posterior.df, cutoff)$component
  se.p2 <- dat$SE[var.categories==2,]

  # grab parameters
  Sigma <- get_var_covar_matrix(fit1)
  p <- get_proportions(fit1)

  # calculate heritability
  list.h <- sapply(1:nrow(Sigma), function(idx)
    .calc_herit_one(se.p2, Sigma, p, idx))
  names(list.h) <- paste(rep("h", length(dat$cols)), dat$cols, sep=".")
  return(list.h)
}

#' Helper function to calculate the heritability for one group (e.g. males or females)
#'
#' @param se.p2 vector of standard errors from the non-null component
#' @param Sigma the variance covariance matrix from an m1 fit
#' @param p the proportion vector from an m1 fit
#' @param idx the index in the list of groups
#' @return the heritability of the trait in that group
.calc_herit_one <- function(se.p2, Sigma, p, idx){
  n <- nrow(se.p2)
  num_i <- n*(p[2])*Sigma[idx,idx]
  h_i <- num_i/(num_i + sum(se.p2[,idx], na.rm=TRUE))
  return(h_i)
}
