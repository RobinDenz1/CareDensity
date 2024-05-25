
## function that projects a bipartite graph into one of its possible
## projections. faster and more RAM efficient than the one in igraph
#' @export
project_to_one_mode <- function(g, mode, sparse=TRUE) {
  
  requireNamespace("igraph")
  requireNamespace("MatrixExtra")
  
  # input checks
  if (!inherits(g, "igraph")) {
    stop("'g' must be an igraph object.")
  } else if (!(length(mode)==1 && is.character(mode) &&
               mode %in% c("cols", "rows"))) {
    stop("'mode' must be either 'cols' or 'rows'.")
  }
  
  # projecting
  bip_mat <- igraph::as_incidence_matrix(g, sparse=sparse)
  if (mode=="cols") {
    out <- MatrixExtra::t(MatrixExtra::crossprod(bip_mat))
  } else if (mode=="rows") {
    out <- MatrixExtra::tcrossprod(bip_mat)
  }
  diag(out) <- 0
  
  return(out)
}
