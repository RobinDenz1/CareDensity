
## given the original input data + an edge list with weights, create
## a data.table with all weights
#' @importFrom data.table fifelse
#' @importFrom data.table as.data.table
#' @importFrom data.table rbindlist
#' @importFrom data.table merge.data.table
#' @importFrom data.table :=
#' @importFrom utils combn
get_d_combs <- function(data, dat_g) {
  
  # silence devtools::check()
  . <- PatID <- ArztID <- combs <- weight <- weight_1 <- NULL
  
  # get all combinations
  d_combs <- data[, .(combs = list(as.data.table(t(combn(ArztID, 2))))),
                  by=PatID]
  d_combs <- d_combs[, rbindlist(combs), by=PatID]
  
  # merge with weight data one way
  colnames(d_combs) <- c("PatID", "from", "to")
  d_combs <- merge.data.table(d_combs, dat_g, by=c("from", "to"), all.x=TRUE)
  
  # merge with weight data the other way
  colnames(d_combs) <- c("to", "from", "PatID", "weight_1")
  d_combs <- merge.data.table(d_combs, dat_g, by=c("from", "to"), all.x=TRUE)
  
  # put together in one weight variable
  d_combs[, weight := fifelse(is.na(weight_1), weight, weight_1)]
  d_combs[, weight_1 := NULL]
  
  return(d_combs)
}

## calculate the classic care density as defined by Pollack
#' @importFrom data.table setDT
#' @importFrom data.table setcolorder
#' @importFrom data.table setkey
#' @importFrom data.table .N
#' @importFrom data.table %chin%
#' @importFrom data.table :=
#' @importFrom igraph graph.data.frame
#' @importFrom igraph bipartite.mapping
#' @importFrom igraph V
#' @importFrom igraph V<-
#' @importFrom igraph graph.adjacency
#' @importFrom igraph get.data.frame
#' @export
care_density <- function(data, pat_col=1, as_data_frame=TRUE) {
  
  requireNamespace("data.table")
  requireNamespace("igraph")
  
  # silence devtools::check()
  . <- PatID <- ArztID <- n <- weight <- sum_weights <- NULL
  
  check_inputs_care_density(data=data, pat_col=pat_col,
                            as_data_frame=as_data_frame)
  
  # as data.table to make it faster
  setDT(data)
  
  # put in right format
  if (pat_col == 1) {
    colnames(data) <- c("PatID", "ArztID")
  } else if (pat_col == 2) {
    colnames(data) <- c("ArztID", "PatID")
    setcolorder(data, c("PatID", "ArztID"))
  }
  
  # only characters work well
  data[, PatID := as.character(PatID)]
  data[, ArztID := as.character(ArztID)]
  data <- unique(data)
  
  # test if unique (would otherwise not create a bipartite graph)
  if (any(data$PatID %in% data$ArztID)) {
    stop("There are patient IDs that are equal to provided IDs, which is not",
         " permissible. Please set unique IDs for patients and providers",
         " and rerun this function.")
  } else if (any(data$AztrID %in% data$PatID)) {
    stop("There are provider IDs that are equal to patient IDs, which is not",
         " permissible. Please set unique IDs for patients and providers",
         " and rerun this function.")
  }
  
  # create bipartite graph from data
  g <- graph.data.frame(data, directed=FALSE)
  V(g)$type <- bipartite.mapping(g)$type
  
  # project to provider mode
  g_prov <- graph.adjacency(adjmatrix=project_to_one_mode(g, mode="cols"),
                            mode="upper", weighted=TRUE, diag=FALSE)
  
  # remove patients with just one connection
  d_counts <- data[, .(n = .N), by=PatID]
  data0 <- data[PatID %chin% d_counts[n == 1]$PatID]
  data <- data[PatID %chin% d_counts[n > 1]$PatID]
  
  # extract weighted edge list
  dat_g <- get.data.frame(g_prov)
  
  # weights of all occured connections
  d_combs <- get_d_combs(data, dat_g)
  
  # put together
  out <- d_combs[, .(sum_weights = sum(weight)), by=PatID]
  out <- merge(out, d_counts, by="PatID", all.x=TRUE)
  
  # apply formula
  out[, care_density := sum_weights / (n * (n - 1) / 2)]
  
  # add patients with just one contact back
  data0[, ArztID := NULL]
  data0[, sum_weights := NA]
  data0[, n := 1]
  data0[, care_density := NA]
  out <- rbind(out, data0)
  setkey(out, "PatID")
  
  if (as_data_frame) {
    out <- as.data.frame(out)
  }
  
  return(out)
}
