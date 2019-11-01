#' Custom assertion function, prints formatted error messages when they fail.
#'
#' Prints an error message if the assertion fails, otherwise does nothing.
#' Error message is of the format (and does not include the function call):
#'    "Error: <error_msg>"
#'
#' @param error_msg error message to print on failure
#' @param test the value to test
.my_assert <- function(error_msg, test){
  if (!test){
    stop(error_msg, call.=FALSE)
  }
}


#' Custom warning function, prints formatted warning messages when they fail.
#'
#' Prints an warning message if the assertion fails, otherwise does nothing.
#' Warning message is of the format (and does not include the function call):
#'    "Warning: <warn_msg>"
#'
#' @param warn_msg warning message to print on failure
#' @param test the value to test
.my_warn <- function(warn_msg, test){
  if (!test){
    warning(warn_msg, call.=FALSE)
  }
}


# ---------- DATA PREP CHECKS ------------ #


.check_gwas_f_format <- function(dat, ID, BETA, SE, CHR, P){
  .my_assert("GWAS file has insufficient columns and/or rows",
             nrow(dat) > 1 & ncol(dat) >= 5)
  .my_assert("the required columns are not present, please rename or specify",
             length(setdiff(c(ID, BETA, SE, CHR, P), colnames(dat)))==0)
}

.check_read_in_params <- function(gwas.files, f.labels, se.cut){

  # ---- check required params ---- #
  # make sure there are a reasonable number of files (2+ and warn if a lot)
  ndim <- length(gwas.files)
  .my_assert("please provide multiple input files", ndim > 1)
  .my_warn("this has not been tested for more than 4 categories", ndim > 4)

  # make sure the files exist
  f_checks <- lapply(gwas.files, function(gwas.file){
    .my_assert("gwas file not found", file_exists(gwas.file))
  })

  # ---- check optional params ---- #
  # if category names are specified, there must be the same number as files
  if (!is.null(f.labels)){
    .my_assert("please provide the same number of labels as files", length(f.labels) == ndim)
  }

  # make sure the standard error specified is a number
  .my_assert("the standard error cutoff must be number between 0 and 1",
             is.numeric(se.cut) & se.cut > 0 & se.cut < 1)

}

.check_var_to_keep <- function(var.keep, data.ids){
  length_overlap <- length(intersect(var.keep, data.ids))
  .my_assert("the variants to keep must overlap the variant ids in the GWAS file",
             length_overlap == 0)
}

.check_overlap_rows <- function(rows.to.keep){
  .my_assert("there are no overlapping variants",
             rows.to.keep == 0)
  .my_warn("you have less than 100 variants after filtering",
           length(rows.to.keep) < 100)
}

# ----------------------------------------- #

# ---------- STAN INPUT CHECKS ------------ #

.all_numeric <- function(df){
  all(apply(df, c(1,2), is.numeric))
}

.check_in_dat_format <- function(dat){

  # make sure the required columns are present
  req_cols <- c("B", "SE")
  .my_assert("Data does not have betas or standard errors. Use prep_gwas_input() to set up input correctly." ,
    all(req_cols %in% names(dat)))
  .my_warn("data is missing id, pval, or chr columns" ,
             "id" %in% names(dat))
  .my_warn("data is missing pval, or chr columns" ,
           all(c("p", "chr") %in% names(dat)))

  # make sure the betas and standard errors are the same size
  .my_assert("Betas and se matrices are not the same size. Use prep_gwas_input() to set up input correctly.",
             nrow(dat$B)==nrow(dat$SE) & ncol(dat$B)==ncol(dat$SE))

  .my_assert("betas and se matrices must only contain numberic values", .all_numeric(dat$B) & .all_numeric(dat$SE))
}

# ----------------------------------------- #

.check_posterior_df <- function(posterior.df){
  .my_assert("Posterior df formatted incorrectly. Please use calc_posteriors() to set up posterior.df correctly.",
    ncol(posterior.df)>=2 & "ID" %in% colnames(posterior.df) & nrow(posterior.df) > 1 )
}

