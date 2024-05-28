
expected <- data.frame(PatID=c("1", "2", "3", "4", "5"),
                       sum_weights=c(5, 3, 3, NA, NA),
                       n=c(3, 2, 2, 1, 1),
                       care_density=c(1.66666667, 3, 3, NA, NA))

test_that("general test case", {
  # some arbitrary patient-provider contact data
  data <- data.frame(PatID=c("1", "1", "1", "2", "2", "3", "3", "4", "5"),
                     ArztID=c("A", "C", "D", "A", "D", "A", "D", "D", "C"))
  
  # calculate the care densities
  out <- care_density(data)
  expect_equal(out, expected)
})

test_that("with pat_col", {
  data <- data.frame(ArztID=c("A", "C", "D", "A", "D", "A", "D", "D", "C"),
                     PatID=c("1", "1", "1", "2", "2", "3", "3", "4", "5"))
  
  out <- care_density(data, pat_col=2)
  expect_equal(out, expected)
})

test_that("returning a data.table", {
  data <- data.frame(PatID=c("1", "1", "1", "2", "2", "3", "3", "4", "5"),
                     ArztID=c("A", "C", "D", "A", "D", "A", "D", "D", "C"))
  
  out <- care_density(data, data_frame=FALSE)
  expect_equal(out, as.data.table(expected, key="PatID"))
})
