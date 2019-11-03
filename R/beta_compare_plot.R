
#' Make a plot to compare the betas for variants assigned to different
#' components. Each point is one variant, each axis is the category. Error
#' bars are standard errors along that axis.
#'
#' @export
#' @param var.df variant data frame from `make_var_table()`
#' @param dat a list with B and SE elements, variant IDs; from `prep_gwas_input()`
#' @param include.null whether to return the table for the null variants, defaults to false
beta_compare_plot <- function(var.df, dat, include.null=FALSE){
  library(ggplot2)
  # rename for plotting
  cat1 <-dat$cols[[1]]
  cat2 <-dat$cols[[2]]
  var.df <- dplyr::rename(var.df, "B1"=sprintf("B.%s", cat1),
                      "B2"=sprintf("B.%s", cat2),
                      "SE1"=sprintf("SE.%s", cat1),
                      "SE2"=sprintf("SE.%s", cat2))

  # plot
  ggplot(var.df, aes(x=B1, y=B2, color=component))+
  geom_hline(yintercept=0)+
  geom_vline(xintercept=0)+
  geom_point(aes(color=component),alpha=0.7)+
  geom_errorbar(aes(ymin=B2-SE2, ymax=B2+SE2, color=component))+
  geom_errorbarh(aes(xmin=B1-SE1, xmax=B1+SE1, color=component))+
  theme(panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_blank())+
  xlab(sprintf("Estimated effect size in %s", cat1))+
  ylab(sprintf("Estimated effect size in %s", cat2))
}
