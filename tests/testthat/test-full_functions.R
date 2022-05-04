# Golden tests


test_that("perm basic produces same results", {

  set.seed(40)

  expect_snapshot_output(

    perm_basic(distribution = 1,
               duplicates = 10,
               set_n = 90,
               set_mean = 100,
               set_sd = 2,
               decimals = 2,
               reps = 500,
               times = 100,
               tolerance = 0)
  )
})



