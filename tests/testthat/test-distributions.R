

test_that("data_in produces same results as R", {

  set.seed(4)

  a <- data_in(distribution = 1,
            set_n = 50,
            set_mean = 100,
            set_sd = 1,
            duplicates = 0)

  set.seed(4)

  b <- rnorm(n = 50,
             mean = 100,
             sd = 1)

  expect_equal(a, b)

})


test_that("data_in with duplicates produces same results as R", {

  set.seed(4)

  a <- data_in(distribution = 1,
               set_n = 50,
               set_mean = 100,
               set_sd = 1,
               duplicates = 10)

  set.seed(4)

  b <- rnorm(n = 50,
             mean = 100,
             sd = 1)

  c <- c(b, b[1:10])

  expect_equal(a, c)

})
