
#' Assign variants to components
#'
#' @param posterior.df posterior.df from calc_posteriors()
#' @param cutoff the minimum cutoff to assign to a non-null component
#' @return the posterior df with an additional 'component' column
assign_to_components <- function(posterior.df, cutoff=0.8){

  if (ncol(posterior.df)==2){
    components <- sapply(posterior.df$p2, function(post)
      ifelse(post > cutoff, 2, 1))
  } else{
    components <- sapply(posterior.df %>% dplyr::select(-ID), function(post){
      max_group <- which.max(post)
      if (post[max_group] <= 0.8){
        max_group <- 1
      }
    })
  }
  return(cbind(posterior.df, "component"=components))
}


