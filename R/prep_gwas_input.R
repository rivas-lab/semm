#' Reads in multiple GWAS files and reformats to a data object that can be used for model fitting.
#'
#' @export
#' @param gwas.files a list of files containing GWAS summary statistics,
#'       this program uses the "ID", "CHR", BETA","SE", and "P" columns for the analysis.
#'       Please rename your columns to match this or specify this as parameters.
#'
#' ---- optional params ----
#' @param cat.labels optionally provide a vector of labels for the files (e.g. c("female", "male")),
#'       defaults to cat1, ..., catn
#' @param se.cut optionally specify the max standard error allowed, defaults to 0.2
#' @param var.keep optional list of variants to keep, often used for LD or other filtering,
#'       defaults to keeping all
#' @param B optional column name with betas, default "B"
#' @param SE optional column name with ses, default "SE"
#' @param ID optional column name with variant ids, default "ID"
#' @param CHR optional column name with chromosome number, default "CHR"
#' @param P optional column name with p-vals, default "P"
#' ------------------------
#' @return a data object containing betas and standard errors formatted
#'        for model fitting, also contains variant ids
prep_gwas_input <- function(gwas.files, cat.labels=NULL, se.cut=0.2, var.keep=NULL,
                         B="BETA", ID="ID", CHR="CHR", SE="SE", P="P"){

  ndim <- length(gwas.files)
  .check_read_in_params(gwas.files, cat.labels, se.cut)
  # create category labels if they are not provided
  if(is.null(cat.labels)) {cat.labels <-paste(rep("cat", ndim), c(1:ndim), sep="")}

  list.ds <- lapply(gwas.files, function(gwas.file) {

    # read in an input file and rename the columns
    dat.0 <- readr::read_tsv(gwas.file)
    .check_gwas_f_format(dat.0, ID, B, SE, CHR, P)

    # add in p-val and chr columns as NAs if not present
    if (!CHR %in% colnames(dat.0)){ dat.0$CHR <- NA}
    if (!P %in% colnames(dat.0)){ dat.0$P <- NA }
    dat.1 <- dplyr::rename(dat.0, "ID"=ID, "B"=B, "SE"=SE, "CHR"=CHR, "P"=P)
    rownames(dat.1) <- dat.1$ID

    # remove NAs
    dat.2 <- dat.1[!is.na(dat.1$SE),]

    # filter based on the standard error cutoff
    dat.3 <- dat.2[dat.2$SE < se.cut,]

    # filter to a subset of variants
    if (!is.null(var.keep)){
      .check_var_to_keep(var.keep,dat.3$ID)
      dat.3 <- dat.3[dat.3$ID %in% var.keep,]
    }

    return(dat.3)
  })

  # extract overlapping rows
  list.ds2 <- .extract_overlapping_rows(list.ds)

  # reformat for stan
  stan.dat <- .prep_dat_stan(list.ds2)
  print(cat.labels)
  stan.dat$cols <- cat.labels # the columns for the data
  return(stan.dat)
}

#' Helper function to extract overlapping rows from a list of datasets.
#' This finds the rows present in all datasets and reorders all datasets
#' to have the same order.
#'
#' @param list.ds a list of datasets with betas and ses
#' @return the reordered and filtered set of datasets
.extract_overlapping_rows <- function(list.ds){
  rows.to.keep <- rownames(list.ds[[1]])
  for (i in 2:length(list.ds)){
    dat <- list.ds[[i]]
    rows.to.keep <- intersect(rows.to.keep, rownames(dat))
  }
  .check_overlap_rows(rows.to.keep)
  print(sprintf("After filtering, we have %s variants remaining", length(rows.to.keep)))

  list.ds2 <- lapply(list.ds, function(x) x[rows.to.keep,])

  return(list.ds2)
}

#' Helper function to prep the dataset for stan
#' @param list.ds a list of datasets, already filtered so the rows are the same
#' @return the data put together such that there is one column per dataset
.prep_dat_stan <- function(list.ds){

  # put together betas and ses, sq se for SE matrix
  betas <- do.call(cbind, lapply(list.ds, function(x) x$B))
  ses <- do.call(cbind, lapply(list.ds, function(x) x$SE))
  pvals <- do.call(cbind, lapply(list.ds, function(x) x$P))
  se2 <- apply(ses, c(1,2), function(x) x^2)

  # grab other information to keep
  ids <- sapply(list.ds[[1]]$ID, as.character)
  chr <- sapply(list.ds[[1]]$CHR, as.character)
  dat <- list("B"=betas, "SE"=se2, "id"=ids, "chr"=chr, "p"=pvals)

  return(dat)
}

