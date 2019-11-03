
#' Calculate the false discovery rate for non-null components given the data and at
#' a particular posterior cutoff. Returns NA for a component if there are fewer than
#' five variants in it.
#'
#' @export
#' @param posterior.df posterior.df from calc_posteriors()
#' @param cutoff the minimum cutoff to assign to a non-null component, default 0.8
#' @return list of false discovery rates at that cutoff for each non-null component
calc_fdr <- function(posterior.df, cutoff=0.8){
  require('dplyr')
  # average null posterior probabilities for variants assigned to a component

  comp.df <- assign_to_components(posterior.df, cutoff)

  # add a p0 column  for model1
  if (ncol(comp.df==2)){
    comp.df <- comp.df %>% dplyr::mutate(p0=1-p1)
  }

  # calculate the FDRs
  fdrs <- comp.df %>%
    filter(component != 0) %>%
    dplyr::group_by(component) %>%
    dplyr::summarize("fdr"=mean(p0), n=n()) %>%
    dplyr::mutate(fdr=ifelse(n < 5, NA, fdr)) %>%
    dplyr::ungroup(component) %>%
    dplyr::select(component, fdr)
  fdr <- fdrs$fdr
  names(fdr) <- paste("p", fdrs$component, sep="")
  return(fdr)
}
