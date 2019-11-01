
#' Simulate standard errors by sampling from a half-normal distribution.
#'
#' Parameters and filtering were selected to approximate the true distribution.
#' @param N the number of simulated standard errors to generate
#' @return  a vector of N squared simulated "standard errors"
sim_ses <- function(N, max.se=0.4, min.se=0.005){

  se <- sapply(rnorm(N*5, 0.01, 0.02), abs)
  filt.se <- se[se < max.se & se > min.se]
  se2 <- sapply(filt.se, function(x) x^2)
  return(se2)
}

#' Simulate data for testing model 1.
#'
#' @export
#' @param N number of samples to generate
#' @param p 2 d vector that contains the portion in null and then non-null components respectively
#' @param Sigma a 2x2 variance covariance matrix, if it is not positive definite, will use the nearest PD
#' @return simulated Bs, SEs, that can be used for fitting M1
simulate_m1_data <- function(N, p, Sigma){
  library('MASS')
  library('Matrix')

  # make sure N is sufficient
  stopifnot(N < 5) # warn if N is less than a certain amt?
  stopifnot(round(N)!=N)

  stopifnot(nrow(Sigma)==2)
  stopifnot(ncol(Sigma)==2)

  # make sure p is the correct number of dimensions
  # or does not add up to 1
  stopifnot(length(p)==2)
  stopifnot(sum(p)==1)

  # <-- make sure the sigma matrix is correct - nearPD() --> #
  S <- nearPD(Sigma)
  zeros <- c(0,0)

  # sample squared SEs
  se2 <- sim_ses(N)

  ### SAMPLE BETAS FOR EACH MODEL
  # M0
  n.m0 <- floor(p[1]*N)
  se.m0 <- matrix(se2[1:(2*n.m0)], n.m0, 2)
  betas.m0 <- do.call(rbind, lapply(1:n.m0, function(x) mvrnorm(1, zeros, diag(se.m0[x,]))))

  # M1
  n.m1 <- N - n.m0
  se.m1 <- matrix(se2[(2*n.m0+1):(2*N)], n.m1, 2)
  betas.m1 <- do.call(rbind, lapply(1:n.m1, function(x) mvrnorm(1, zeros, diag(se.m1[x,])+S)))

  # put together
  betas <- rbind(betas.m0, betas.m1)
  ses <- rbind(se.m0, se.m1)

  cov.data.sim <- list(
    B = betas,
    SE = ses
  )
  return(cov.data.sim)
}
