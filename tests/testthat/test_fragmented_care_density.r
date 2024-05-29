
expected <- data.frame(PatID=c("1", "2", "3", "4", "5"),
                       sum_weights=c(5, 3, 3, NA, NA),
                       n=c(3, 2, 2, 1, 1),
                       care_density=c(1.66666667, 3, 3, NA, NA))

data <- data.frame(PatID=c("1", "1", "1", "2", "2", "3", "3", "4", "5"),
                   ArztID=c("A", "C", "D", "A", "D", "A", "D", "D", "C"))

d_type <- data.frame(ID=c("A", "C", "D"),
                     Type=c("GP", "GP", "Psychiatrist"))

test_that("with equal weights, equal to care_density", {
  d_weights <- data.frame(from=c("GP", "GP", "Psychiatrist"),
                          to=c("GP", "Psychiatrist", "Psychiatrist"),
                          weight=c(1, 1, 1))
  
  out <- fragmented_care_density(data, type=d_type, weights=d_weights)
  expect_equal(out$fragmented_care_density, care_density(data)$care_density)
})

test_that("with unequal weights", {
  d_weights <- data.frame(from=c("GP", "GP", "Psychiatrist"),
                          to=c("GP", "Psychiatrist", "Psychiatrist"),
                          weight=c(0.3, 1.6, 1.1))
  out <- fragmented_care_density(data, type=d_type, weights=d_weights)
  expected <- data.frame(PatID=c("1", "2", "3", "4", "5"),
                         fragmented_care_density=c(2.23333333, 4.8, 4.8, NA,
                                                   NA))
  expect_equal(out, expected)
})

test_that("using byconnection=TRUE", {
  d_weights <- data.frame(from=c("GP", "GP", "Psychiatrist"),
                          to=c("GP", "Psychiatrist", "Psychiatrist"),
                          weight=c(0.3, 1.6, 1.1))
  out <- fragmented_care_density(data, type=d_type, weights=NULL,
                                 by_connection=TRUE)
  expected <- data.frame(PatID=c("1", "1", "2", "3", "4", "5"),
                         connection=c("GP - GP", "Psychiatrist - GP",
                                      "Psychiatrist - GP", "Psychiatrist - GP",
                                      NA_character_, NA_character_),
                         sum_weights=c(1, 4, 3, 3, NA, NA),
                         n=c(3, 3, 2, 2, 1, 1),
                         care_density=c(0.33333333, 1.333333333, 3, 3, NA, NA))
  expect_equal(out, expected)
})
