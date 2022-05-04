
test_that("average_fre function works", {

  x <- c(1, 1, 1, 2)
  a <- average_fre(x)

  #Square everything and then divide by total number
  b <- (3^2 + 1^1) / 4
  expect_equal(a, b)
})


test_that("average_fre2 function works", {

  x <- c(3, 1)
  a <- average_fre2(x, 4)

  #Square everything and then divide by total number
  b <- (3^2 + 1^1) / 4
  expect_equal(a, b)
})



test_that("g2 function works with 0 probabilities", {

  actual <- c(0.4, 0.2, 0.2, 0)
  expected <- c(0.5, 0.25, 0.25, 0)
  a <- g2_stat(50, actual, expected)

  expect_false(is.na(a))
})


test_that("chisq function works with 0 probabilities", {

  actual <- c(0.4, 0.2, 0.2, 0)
  expected <- c(0.5, 0.25, 0.25, 0)
  a <- chisq_stat(50, actual, expected)

  expect_false(is.na(a))
})

test_that("ft function works with 0 probabilities", {

  actual <- c(0.4, 0.2, 0.2, 0)
  expected <- c(0.5, 0.25, 0.25, 0)
  a <- ft_stat(50, actual, expected)

  expect_false(is.na(a))
})


