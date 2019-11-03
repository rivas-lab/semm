
#' Assign variants to components
#'
#' @param posterior.df posterior.df from calc_posteriors()
#' @param cutoff the minimum cutoff to assign to a non-null component, default 0.8
#' @return the posterior df with an additional 'component' column
assign_to_components <- function(posterior.df, cutoff=0.8){

  if (ncol(posterior.df)==2){
    components <- sapply(posterior.df$p1, function(post)
      ifelse(post > cutoff, 1, 0))
  } else{
    components <- apply(posterior.df, 1, function(post){
      max_group <- c(0:3)[which.max(post[1:4])]
      if (post[max_group+1] <= 0.8){
        max_group <- 0
      }
      return(max_group)
    })
  }
  return(cbind(posterior.df, "component"=unlist(components)))
}


