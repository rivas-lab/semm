#' Reads in multiple GWAS files and reformats to a data object that can be used for model fitting.
#'
#' @export
#' @param gwas.files a list of locations of files containing GWAS summary statistics,
#'       this program uses the "ID", "CHR", BETA","SE", and "P" columns for the analysis.
#'       Please rename your columns to match this or specify this as parameters.
#'
#' ---- optional params ----
#' @param f.labels optionally provide a vector of labels for the files (e.g. c("female", "male")),
#'       defaults to cat1, ..., catn
#' @param se.cut optionally specify the max standard error allowed, defaults to 0.2
#' @param var.keep optional list of variants to keep, often used for LD or other filtering,
#'       defaults to keeping all
#' @param BETA optional column name with betas, default "BETA"
#' @param SE optional column name with ses, default "SE"
#' @param ID optional column name with variant ids, default "ID"
#' @param CHR optional column name with chromosome number, default "CHR"
#' @param P optional column name with p-vals, default "P"
#' ------------------------
#' @return a data object containing betas and standard errors formatted
#'        for model fitting, also contains variant ids
prep_gwas_input <- function(gwas.files, f.labels="", se.cut=0.2, var.keep="",
                         BETA="BETA", ID="ID", CHR="CHR", SE="SE", P="P"){

  .check_read_in_params(gwas.files, f.labels, se.cut)

  # create category labels if they are not provided
  f.labels <- ifelse(f.labels=="",paste(rep("cat", ndim), c(1:ndim), sep=""), f.labels)

  list.ds <- lapply(gwas.files, function(gwas.file) {

    # read in an input file and rename the columns
    dat.0 <- read.table(gwas.file, header=TRUE, sep=" ")
    .check_gwas_f_format(dat.0, ID, BETA, SE, CHR, P)

    dat.1 <- dat.0 %>% rename("ID"=ID, "BETA"=BETA, "SE"=SE, "CHR"=CHR, "P"=P)
    rownames(dat.1) <- dat.1$ID

    # remove NAs
    dat.2 <- dat.1[!is.na(dat.1$SE),]

    # filter based on the standard error cutoff
    dat.3 <- dat.2[dat.2$SE < se.cut,]

    # filter to a subset of variants
    if (var.keep != ""){
      .check_var_to_keep(var.keep, dat3$ID)
      dat.3 <- dat.3[dat.3$ID %in% var.keep,]
    }

    return(dat.3)
  })

  # extract overlapping rows
  list.ds2 <- .extract_overlapping_rows(list.ds)

  # reformat for stan
  stan.dat <- .prep_dat_stan(list.ds2)
  stan.dat$cols <- f.labels # the columns for the data
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

  list.ds2 <- lapply(list.ds, function(x) x[rows.to.keep,])

  return(list.ds2)
}

#' Helper function to prep the dataset for stan
#' @param list.ds a list of datasets, already filtered so the rows are the same
#' @return the data put together such that there is one column per dataset
.prep_dat_stan <- function(list.ds){

  # put together betas and ses, sq se for SE matrix
  betas <- do.call(cbind, lapply(list.ds, function(x) x$BETA))
  ses <- do.call(cbind, lapply(list.ds, function(x) x$SE))
  pvals <- do.call(cbind, lapply(list.ds, function(x) x$P))
  se2 <- apply(ses, c(1,2), function(x) x^2)

  # grab other information to keep
  ids <- sapply(list.ds[[1]]$ID, as.character)
  chr <- sapply(list.ds[[1]]$CHR, as.character)
  dat <- list("B"=betas, "SE"=se2, "id"=ids, "chr"=chr, "p"=pvals)

  return(dat)
}

