
## input checks for care_density() function
check_inputs_care_density <- function(data, pat_col, as_data_frame) {
  
  if (!inherits(data, "data.frame")) {
    stop("'data' must be a data.frame like object.")
  } else if (ncol(data) != 2) {
    stop("'data' must contain exactly two columns.")
  } else if (nrow(data) == 0) {
    stop("'data' may not be empty.")
  } else if (!(length(pat_col)==1 && is.numeric(pat_col) &&
               pat_col %in% c(1, 2))) {
    stop("'pat_col' must be either 1 or 2.")
  } else if (!(length(as_data_frame)==1 && is.logical(as_data_frame))) {
    stop("'as_data_frame' must be either TRUE or FALSE.")
  }
}
