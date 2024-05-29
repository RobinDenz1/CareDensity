
# some arbitrary patient-provider contact data
data <- data.frame(PatID=c("1", "1", "1", "2", "2", "3", "3", "4", "5"),
                   ArztID=c("A", "C", "D", "A", "D", "A", "D", "D", "C"))

# create graph
g <- graph_from_data_frame(data, directed=FALSE)

# add type
V(g)$type <- bipartite_mapping(g)$type

test_that("row mode", {
  out <- project_to_one_mode(g, mode="rows")
  out2 <- project_to_one_mode(g, mode="rows", sparse=FALSE)
  expect_true(ncol(out)==5)
  expect_true(ncol(out2)==5)
  expect_true(nrow(out)==5)
  expect_true(nrow(out2)==5)
  expect_equal(as.matrix(out), out2)
})

test_that("column mode", {
  out <- project_to_one_mode(g, mode="cols")
  out2 <- project_to_one_mode(g, mode="cols", sparse=FALSE)
  expect_true(ncol(out)==3)
  expect_true(ncol(out2)==3)
  expect_true(nrow(out)==3)
  expect_true(nrow(out2)==3)
  expect_equal(as.matrix(out), out2)
})

test_that("wrong g", {
  expect_error(project_to_one_mode("nope", mode="rows"))
})

test_that("wrong mode", {
  expect_error(project_to_one_mode(g, mode="row"))
})
