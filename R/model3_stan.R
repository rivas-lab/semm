#' Identify sex-specific variants using SEMM.
#'
#' This function takes as input paired betas and standard errors from GWAS summary
#' statistics and outputs a fitted four component model that can be used to assign
#' variants to components.
#'
#' @export
#' @param dat input data object, produced from prep_gwas_input()
#' @param B optional, if dat is not provided,
#'          a data frame containing beta estimates from GWAS summary statistics.
#'          each row corresponds to a variant and each column corresponds to
#'          a group (e.g. males and females)
#' @param SE optional, if dat is not provided,
#'          a data frame containing SE estimates from GWAS summary statistics,
#'           each row/column refers to the same variant and group as the beta values
#' @return An object of class `stanfit` returned by `rstan::sampling` fitted to model2
#'
model3_stan <- function(dat=NULL, B=NULL, SE=NULL, ...){

  # check the input
  .my_assert("please provide either a dataset -or- betas and SEs",
             (is.null(dat) & !is.null(B) & !is.null(SE)) |
               (!is.null(dat) & is.null(B) & is.null(SE)))
  if (is.null(dat)){
    dat=list("B"=B, "SE"=SE)
  }
  .check_in_dat_format(dat)
  .my_assert("model3 only works for data with two categories", ncol(dat$B)==2)

  # fit the model
  standata <- list(B=dat$B, SE=dat$SE, N=nrow(dat$B), M = ncol(dat$B), K = 4)
  out <- rstan::sampling(stanmodels$model2, data = standata, ...)
  return(out)
}
