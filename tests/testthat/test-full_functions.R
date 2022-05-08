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


  test_that("td_indepedence produces same result as a test
            based on r2dtables", {

              set.seed(490)

              x <- rnorm(300, mean = 10, sd = 14)

              #For one decimal

              sam <- as.integer(x * 10)
              int <- as.integer(x)
              dec <- sam - int * 10


              tab <- table(sam)
              tab_int <- table(int)
              tab_dec <- table(dec)

              u_int <- sort(unique(int))
              u_dec <- sort(unique(dec))

              #This creates a vector of counts for all possible cells

              observed <- data.frame(v = 1:(length(u_int) * length(u_dec)))

              count <- 1

              for(i in u_int)  {

                for(j in u_dec)  {

                  num <- as.character((i * 10) + j)

                  observed$num[count] <- num
                  observed$value[count] <- as.integer(tab[num])

                  if(is.na(observed$value[count])) {

                    observed$value[count] <- 0

                  }

                  count <- count + 1

                }}

              mat <- matrix(observed$value, byrow = TRUE, nrow = length(u_int))

              #Running tests:

              set.seed(202)

              a <- td_independence(x = sam * 0.1,
                                   decimals = 1,
                                   reps = 10000)

              b <- chisq.test(mat, simulate.p.value = TRUE, B = 10000)

              expect_equal(a, b, tolerance = 0.01)
            })





