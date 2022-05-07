# Function to use across several tests.

out_vector_r <- function(c_sum)  {

  dat <- NULL

  for(i in 1:length(c_sum)) {

    x <- (rep(i, c_sum[i]))

    dat <- append(dat, x)

  }

  return(dat)
}


test_that("full_vec = example_1", {

  x <- c(20, 40, 20, 20, 40, 90)
  boot <- c(0:9)

  a <- full_vec(x, boot)

  b <- c(20:29, 40:49, 90:99)

  expect_equal(a, b)
})


test_that("full_vec = correct length", {

  x <- c(1:8)
  boot <- c(1, 2, 5, 7, 8)

  a <- length(full_vec(x, boot))

  b <- 8 * 5

  expect_equal(a, b)
})

test_that("vector_out = example_1", {

  x <- c(1:7)

  #First we produce the function in R. Each number we generate
  #here represents a column sum. So if the first column has a sum
  #of 8, then there will be 8 1's, if the second column has a sum
  #of 4, then there will be 4 2's, etc.

  a <- out_vector_r(x)

  b <- out_vector_cpp(x) + 1 #b/c this indexes starting at 0 instead of 1

  expect_equal(a, b)

  })


test_that("expected_cells = 1", {

  m <- matrix(1:12, nrow = 3, ncol = 4)

  r_sum <- rowSums(m) / sum(rowSums(m))
  c_sum <- colSums(m) / sum(colSums(m))

  expect_equal(1, sum(expected_cells(r_sum, c_sum)))

})

test_that("expected_cells = example 1", {

  x <- c(1.0, 1.1, 1.2)

  int <- as.integer(x)
  dec <- (x - int) * 10

  r_frac <- table(int)/sum(table(int))
  c_frac <- table(dec)/sum(table(dec))

  a <- c(1/3, 1/3, 1/3)

  b <- expected_cells(r_frac, c_frac)

  expect_equal(a, b)

})


test_that("expected_cells = example 2", {

  x <- c(2.0, 1.0, 1.1, 1.2)

  int <- as.integer(x)
  dec <- (x - int) * 10

  r_frac <- table(int)/sum(table(int))
  c_frac <- table(dec)/sum(table(dec))

  a <- c(3/4 * c(1/2, 1/4, 1/4),
         1/4 * c(1/2, 1/4, 1/4))

  b <- expected_cells(r_frac, c_frac)

  expect_equal(a, b)

})


test_that("expected_cells = example 3", {

  #Only works if 10 decimals are present.

  #x <- c(1.0, 1.2, 1.3, 2.1, 2.8, 2.9, 2.9, 2.4, 2.5, 2.6, 2.7)
  x <- c(20:30, 21, 31, 22) / 10

  int <- as.integer(x)
  dec <- (x - int) * 10

  r_frac <- table(int)/sum(table(int))
  c_frac <- table(dec)/sum(table(dec))

  #a <- bunch_prob(dec, int)

  b <- expected_cells(r_frac, c_frac)

  d <- as.numeric(c(r_frac[1] * c_frac, r_frac[2] * c_frac))

  expect_equal(d, b)

})



test_that("perm_vector samples approximate r2dtable",  {

  chisq_s <- function(x)  {

    n <- sum(x)
    nr <- as.integer(nrow(x))
    nc <- as.integer(ncol(x))
    sr <- rowSums(x)
    sc <- colSums(x)
    E <- outer(sr, sc, "*") / n

    STATISTIC <- sum(abs((x - E))^2 / E)

    return(STATISTIC)

  }


  #r_sum <- c(10, 20, 30, 40)
  #c_sum <- c(40, 30, 20, 10)

  r_sum <- c(1:4)
  c_sum <- c(4:1)

  v <- out_vector_cpp(c_sum)

  my_chi <- 1:10000
  r_chi <- 1:10000

  set.seed(404)

  for(i in 1:10000)  {

    a <- perm_vector(v + 1, r_sum, c_sum)

    b <- matrix(a, byrow = TRUE, ncol = 4)

    my_chi[i] <- chisq_s(b)

    d <- r2dtable(1, r_sum, c_sum)

    e <- matrix(unlist(d), ncol = 4)

    r_chi[i] <- chisq_s(e)

  }

  p <- chisq.test(my_chi, r_chi)$p.value

  #Looking at whole distribution

  expect_gt(p, 0.05)

})



test_that("perm_vector equal rowSums", {

  mat <- matrix(data = 1:16, ncol = 4)

  r_sum <- rowSums(mat)
  c_sum <- rowSums(mat)

  v <- out_vector_cpp(c_sum)
  a <- perm_vector(v + 1, r_sum, c_sum)
  b <- matrix(a, ncol = 4)
  c <- rowSums(b)

  expect_equal(r_sum, c)

})

test_that("int_dec works: one decimal", {

  x <- c(1.12345, 2.12345, 3.82345)

  a <- int_dec(x, decimals = 1)

  #Sample
  b <- c(11, 21, 38)

  # Integer
  c <- c(10, 20, 30)

  # Decimals
  d <- c(1, 1, 8)

  expect_equal(a$s_f, b)
  expect_equal(a$int_f, c)
  expect_equal(a$dec_f, d)

})

test_that("int_dec works: two decimals", {

  x <- c(1.12345, 2.12345, 3.85345)

  a <- int_dec(x, decimals = 2)

  #Sample
  b <- c(112, 212, 385)

  # Integer
  c <- c(110, 210, 380)

  # Decimals
  d <- c(2, 2, 5)

  expect_equal(a$s_f, b)
  expect_equal(a$int_f, c)
  expect_equal(a$dec_f, d)

})

test_that("int_dec works: three decimals", {

  x <- c(1.12345, 2.12745, 3.85345)

  a <- int_dec(x, decimals = 3)

  #Sample
  b <- c(1123, 2127, 3853)

  # Integer
  c <- c(1120, 2120, 3850)

  # Decimals
  d <- c(3, 7, 3)

  expect_equal(a$s_f, b)
  expect_equal(a$int_f, c)
  expect_equal(a$dec_f, d)

})


test_that("actual_frac works", {

  p <- 1:10

  q <- c(1, 1, 3, 4, 5, 6, 6, 6)

  a <- actual_frac(p, q, new_n = 20)

  b <- c(0.1, 0, 0.05, 0.05, 0.05, 0.15, 0, 0, 0, 0)

  expect_equal(a, b)

})


test_that("observed_vec gives expect result", {

  u_int <- c(-2, 3, 4, 13)
  u_dec <- c(1, 2, 3, 4)
  u_sam <- c(-21, -24, 33, 42, 131, 132, 134)
  tab_sam <- c(1, 1, 1, 2, 1, 1, 1)

  #For vec: -21, -22, -23, -24,
  #         31, 32, 33, 34,
  #         41, 42, 43, 44,
  #         81, 82, 83, 84

  a <- c(1, 0, 0, 1, 0, 0, 1, 0, 0, 2, 0, 0, 1, 1, 0, 1)

  b <- observed_vec(u_int, u_dec, u_sam, tab_sam)

  expect_equal(a, b)


})


test_that("observed_vec gives same result as in R", {

              set.seed(490)

              x <- rnorm(300, mean = 54, sd = 14)

              #For one decimal

              sam <- as.integer(x * 10)
              int <- as.integer(x)
              dec <- sam - int * 10

              tab_sam <- table(sam)
              tab_int <- table(int)
              tab_dec <- table(dec)

              u_int <- sort(unique(int))
              u_dec <- sort(unique(dec))
              u_sam <- sort(unique(sam))

              #This creates a vector of counts for all possible cells

              observed <- data.frame(v = 1:(length(u_int) * length(u_dec)))

              count <- 1

              for(i in u_int)  {

                for(j in u_dec)  {

                  num <- as.character((i * 10) + j)

                  observed$num[count] <- num
                  observed$value[count] <- as.integer(tab_sam[num])

                  if(is.na(observed$value[count])) {

                    observed$value[count] <- 0

                  }

                  count <- count + 1

                }}

              a <- observed$value

              b <- observed_vec(u_int, u_dec, u_sam, tab_sam)

              expect_equal(a, b)

            })

