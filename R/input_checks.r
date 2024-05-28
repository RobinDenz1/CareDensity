
## input checks for care_density() function
check_inputs_care_density <- function(data, pat_col, data_frame) {
  
  if (!inherits(data, "data.frame")) {
    stop("'data' must be a data.frame like object.")
  } else if (ncol(data) != 2) {
    stop("'data' must contain exactly two columns.")
  } else if (nrow(data) == 0) {
    stop("'data' may not be empty.")
  } else if (!(length(pat_col)==1 && is.numeric(pat_col) &&
               pat_col %in% c(1, 2))) {
    stop("'pat_col' must be either 1 or 2.")
  } else if (!(length(data_frame)==1 && is.logical(data_frame))) {
    stop("'data_frame' must be either TRUE or FALSE.")
  }
}

## input checks for fragmented_care_density()
#' @importFrom data.table %chin%
#' @importFrom data.table fifelse
check_inputs_fcd <- function(data, pat_col, data_frame, type,
                             weights, by_connection) {
  
  check_inputs_care_density(data=data, pat_col=pat_col,
                            data_frame=data_frame)
  
  prov_col <- fifelse(pat_col==1, 2, 1)
  data <- as.data.frame(data)
  
  if (!(length(by_connection)==1 && is.logical(by_connection))) {
    stop("'by_connection' must be either TRUE or FALSE.")
  } else if (!(is.data.frame(type) && ncol(type)==2 &&
               all(colnames(type)==c("ID", "Type")) &&
               is.character(type$ID) && is.character(type$Type))) {
    stop("'type' must be a data.frame with only two columns called",
         " 'ID' and 'Type', both including only characters.")
  } else if (!all(data[,prov_col][[1]] %chin% type$ID)) {
    missing_prov <- data[,prov_col][[1]][!data[,prov_col][[1]] %chin% type$ID]
    stop("The 'type' data.frame must contain information on all available",
         " providers. Information not included for: ", missing_prov)
  } else if (!by_connection && is.null(weights)) {
    stop("'weights' must be a data.frame containing connection-specific",
         " weights if by_connection=FALSE.")
  } else if (!(by_connection || (is.data.frame(weights) && ncol(weights)==3 &&
               all(colnames(weights)==c("from", "to", "weight")) &&
               is.character(weights$from) && is.character(weights$to) &&
               is.numeric(weights$weight)))) {
    stop("'weights' must be a data.frame with exactly three columns called",
         " 'from' (character), 'to' (character) and 'weight' (numeric).")
  }
}
