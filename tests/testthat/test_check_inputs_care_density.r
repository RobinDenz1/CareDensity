
data <- data.frame(PatID=c("1", "1", "1", "2", "2", "3", "3", "4", "5"),
                   ArztID=c("A", "C", "D", "A", "D", "A", "D", "D", "C"))

test_that("wrong data", {
  expect_error(care_density(data="nah"),
               "'data' must be a data.frame like object.")
  expect_error(care_density(data=data.frame(ID=data$PatID)),
               "'data' must contain exactly two columns.")
  
  empty_dat <- data.frame(matrix(ncol=2, nrow=0))
  colnames(empty_dat) <- c("PatID", "ArztID")
  expect_error(care_density(data=empty_dat), "'data' may not be empty.")
})

test_that("wrong pat_col", {
  expect_error(care_density(data, pat_col="PatID"),
               "'pat_col' must be either 1 or 2.")
})

test_that("wrong data_frame", {
  expect_error(care_density(data, data_frame="table"),
               "'data_frame' must be either TRUE or FALSE.")
})

test_that("PatID in ArztID", {
  data <- data.frame(PatID=c("1", "1", "1", "2", "2", "3", "3", "4", "5"),
                     ArztID=c("A", "C", "D", "A", "D", "A", "D", "D", "1"))
  expect_error(care_density(data),
               paste0("There are patient IDs that are equal to provider ",
                      "IDs, which is not permissible. Please set unique ",
                      "IDs for patients and providers and rerun ",
                      "this function."))
})
