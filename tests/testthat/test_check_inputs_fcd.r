
# some arbitrary patient-provider contact data
data <- data.frame(PatID=c("1", "1", "1", "2", "2", "3", "3", "4", "5"),
                   ArztID=c("A", "C", "D", "A", "D", "A", "D", "D", "C"))

# defining the provider types
d_type <- data.frame(ID=c("A", "C", "D"),
                     Type=c("GP", "GP", "Psychiatrist"))

# defining the connection-specific weights
d_weights <- data.frame(from=c("GP", "GP", "Psychiatrist"),
                        to=c("GP", "Psychiatrist", "Psychiatrist"),
                        weight=c(1.1, 0.8, 1.3))

test_that("wrong by_connection", {
  expect_error(fragmented_care_density(data, type=d_type, weights=d_weights,
                                       by_connection="GP - GP"),
               "'by_connection' must be either TRUE or FALSE.")
})

test_that("wrong type structure", {
  d_type2 <- d_type
  colnames(d_type2) <- c("a", "b")
  expect_error(fragmented_care_density(data, type=d_type2, weights=d_weights,
                                       by_connection=FALSE),
               paste0("'type' must be a data.frame with only two columns ",
                      "called 'ID' and 'Type', both including ",
                      "only characters."))
})

test_that("missing type for some providers", {
  d_type2 <- subset(d_type, ID != "A")
  expect_error(fragmented_care_density(data, type=d_type2, weights=d_weights,
                                       by_connection=FALSE),
               paste0("The 'type' data.frame must contain information on ",
                      "all available providers. Information not ",
                      "included for: A"))
})

test_that("no weights with by_connection=FALSE", {
  expect_error(fragmented_care_density(data, type=d_type, weights=NULL,
                                       by_connection=FALSE),
               paste0("'weights' must be a data.frame containing ",
                      "connection-specific weights if by_connection=FALSE."))
})

test_that("wrong weights structure", {
  d_weights2 <- d_weights
  colnames(d_weights2) <- c("a", "b", "c")
  expect_error(fragmented_care_density(data, type=d_type, weights=d_weights2,
                                       by_connection=FALSE),
               paste0("'weights' must be a data.frame with exactly three ",
                      "columns called 'from' (character), 'to' ",
                      "(character) and 'weight' (numeric)."), fixed=TRUE)
})

test_that("undefined connection weights", {
  d_weights2 <- subset(d_weights, !(from=="GP" & to=="GP"))
  expect_error(fragmented_care_density(data, type=d_type, weights=d_weights2,
                                       by_connection=FALSE),
               paste0("The following connection types are not defined in ",
               "the 'weights' argument: GP - GP.
 Please add those to the 'weight' data.frame and rerun this function."))
})
