
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

## given a data.frame of weigts per connections, expand it to include
## all equal connections (direction does not matter)
#' @importFrom data.table melt.data.table
expand_connection_weights <- function(weights) {
  
  con <- from <- to <- variable <- con_rev <- NULL
  
  setDT(weights)
  
  # get connection types
  weights[, con := paste0(from, " - ", to)]
  weights[, con_rev := paste0(to, " - ", from)]
  weights[, from := NULL]
  weights[, to := NULL]
  
  # put in long format and remove duplicates
  weights <- melt.data.table(weights, id.vars="weight")
  weights[, variable := NULL]
  weights <- unique(weights)
  colnames(weights) <- c("con_weight", "connection")
  
  return(weights)
}

## general function to calculate either the care density or the fragmented
## care density
#' @importFrom data.table setDT
#' @importFrom data.table setcolorder
#' @importFrom data.table setkey
#' @importFrom data.table .N
#' @importFrom data.table %chin%
#' @importFrom data.table :=
#' @importFrom igraph graph_from_data_frame
#' @importFrom igraph bipartite_mapping
#' @importFrom igraph V
#' @importFrom igraph V<-
#' @importFrom igraph graph_from_adjacency_matrix
#' @importFrom igraph as_data_frame
care_density.fit <- function(data, pat_col=1, weights, type,
                             by_connection=FALSE, data_frame=TRUE,
                             estimator) {
  
  requireNamespace("data.table")
  requireNamespace("igraph")
  
  # silence devtools::check()
  . <- PatID <- ArztID <- n <- weight <- sum_weights <- NULL
  connection <- type_from <- type_to <- con_weight <- NULL
  
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
  if (any(data$PatID %chin% data$ArztID)) {
    stop("There are patient IDs that are equal to provider IDs, which is not",
         " permissible. Please set unique IDs for patients and providers",
         " and rerun this function.")
  }
  
  # create bipartite graph from data
  g <- graph_from_data_frame(data, directed=FALSE)
  V(g)$type <- bipartite_mapping(g)$type
  
  # project to provider mode
  mat <- project_to_one_mode(g, mode="cols")
  g_prov <- graph_from_adjacency_matrix(adjmatrix=mat, mode="upper",
                                        weighted=TRUE, diag=FALSE)
  
  # remove patients with just one connection
  d_counts <- data[, .(n = .N), by=PatID]
  data0 <- data[PatID %chin% d_counts[n == 1]$PatID]
  data <- data[PatID %chin% d_counts[n > 1]$PatID]
  
  # extract weighted edge list
  dat_g <- as_data_frame(g_prov)
  
  # weights of all occurred connections
  d_combs <- get_d_combs(data, dat_g)
  
  if (estimator=="simple") {
    # put together
    out <- d_combs[, .(sum_weights = sum(weight)), by=PatID]
    out <- merge(out, d_counts, by="PatID", all.x=TRUE)
    
    # apply formula
    out[, care_density := sum_weights / (n * (n - 1) / 2)]
    
  } else if (estimator=="fragmented") {
    
    # merge types to it
    colnames(type) <- c("from", "type_from")
    d_combs <- merge(d_combs, type, by="from", all.x=TRUE)
    
    colnames(type) <- c("to", "type_to")
    d_combs <- merge(d_combs, type, by="to", all.x=TRUE)
    
    # create connection type ID
    d_combs[, connection := paste0(type_from, " - ", type_to)]
    d_combs[, type_from := NULL]
    d_combs[, type_to := NULL]
    
    # calculate sums per connection type
    out <- d_combs[, .(sum_weights = sum(weight)), by=list(PatID, connection)]
    
    # merge with n_p data
    out <- merge(out, d_counts, by="PatID", all.x=TRUE)
    
    # calculate care density per connection type
    out[, care_density := sum_weights / (n * (n - 1) / 2)]
    
    # calculate final fragmented care density
    if (!by_connection) {
      
      # expand to get redundant connections as well
      d_weights <- expand_connection_weights(weights)
      
      # check if some are missing
      missing_con <- out$connection[!out$connection %chin% d_weights$connection]
      missing_con <- unique(missing_con)
      
      if (length(missing_con) > 0) {
        stop("The following connection types are not defined in the 'weights'",
             " argument: ", paste0(missing_con, collapse=", "), ".\n Please",
             " add those to the 'weight' data.frame and rerun this function.")
      }
      
      # merge with connection weights
      out <- merge(out, d_weights, by="connection", all.x=TRUE)
      
      # calculate FCD
      out <- out[, .(fragmented_care_density = sum(care_density*con_weight)),
                 by=PatID]
    }
  }
  
  # add patients with just one contact back
  data0[, ArztID := NULL]
  
  if (estimator=="simple") {
    data0[, sum_weights := NA]
    data0[, n := 1]
    data0[, care_density := NA]
  } else if (by_connection) {
    data0[, connection := NA_character_]
    data0[, sum_weights := NA]
    data0[, n := 1]
    data0[, care_density := NA]
  } else {
    data0[, fragmented_care_density := NA]
  }
  
  out <- rbind(out, data0)
  setkey(out, "PatID")
  
  if (data_frame) {
    out <- as.data.frame(out)
  }
  
  return(out)
}

## calculate the classic care density as defined by Pollack et al. (2013)
#' @export
care_density <- function(data, pat_col=1, data_frame=TRUE) {
  
  check_inputs_care_density(data=data, pat_col=pat_col,
                            data_frame=data_frame)
  
  out <- care_density.fit(data=data, pat_col=pat_col,
                          data_frame=data_frame,
                          by_connection=FALSE, weights=NULL, type=NULL,
                          estimator="simple")
  return(out)
}

## calculate the fragmented care density as defined by Engels et al. (2024)
#' @export
fragmented_care_density <- function(data, pat_col=1, weights, type,
                                    by_connection=FALSE, data_frame=TRUE) {
  
  check_inputs_fcd(data=data, pat_col=pat_col, weights=weights, type=type,
                   by_connection=by_connection, data_frame=data_frame)
  
  out <- care_density.fit(data=data, pat_col=pat_col,
                          data_frame=data_frame,
                          by_connection=by_connection, weights=weights,
                          type=type, estimator="fragmented")
  return(out)
}
