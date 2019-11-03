#' Make a variant table to easily view the output of m2
#'
#' @export
#' @param dat a list with B and SE elements, variant IDs; from prep_gwas_input()
#' @param posterior.df posterior.df from calc_posteriors()
#' @param cutoff the minimum cutoff to assign to a non-null component, default 0.8
#' @param include.null whether to return the table for the null variants, defaults to false
#' @return a table specifying the variants, their betas, pvals, and posteriors
make_var_table <- function(dat, posterior.df, cutoff=0.8, include.null=FALSE){
  component.df <- assign_to_components(posterior.df, cutoff)

  # make the categories readable names
  component.df$component[component.df$component==0] <- "NULL"
  component.df$component[component.df$component==1] <- dat$cols[[1]]
  component.df$component[component.df$component==2] <- dat$cols[[2]]
  component.df$component[component.df$component==3] <- "shared"
  if (!include.null){
    component.df <- dplyr::filter(component.df, component!="NULL")
  }
  component.df$component <- factor(component.df$component)

  df <- data.frame(do.call(cbind, list(dat$id, dat$chr, dat$B, dat$SE, dat$p)))
  df[,2:ncol(df)] <- apply(df[,2:ncol(df)], c(1,2), as.numeric)
  colnames(df) <- c("ID", "CHR",paste(c("B", "B", "SE", "SE","P", "P"), dat$cols, sep="."))
  if (!include.null) {
    component.df <- dplyr::filter(component.df, component != 0)
  }

  var.table <- dplyr::right_join(df, component.df, by="ID")
  return(var.table)
}
