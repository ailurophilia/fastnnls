
test_that("Entirely positive example works", {
  set.seed(0)
  d <- matrix(c(1,2,3),ncol = 1)
  Z <- matrix(rnorm(300),ncol = 3)
  x <- rnorm(100) + Z%*%d
  expect_equal(as.numeric(fast_nnls(ZTx = t(Z)%*%x,
                         ZTZ = t(Z)%*%Z)),
               optim(c(0,0,0),
                     function(d) sum((x - Z%*%matrix(d,ncol = 1))^2),
                     method = "L-BFGS-B",
                     lower = 0,
                     upper = 1000)$par)

})

test_that("Mixed example works", {
  set.seed(0)
  d <- matrix(c(-1,2,3),ncol = 1)
  Z <- matrix(rnorm(300),ncol = 3)
  x <- rnorm(100) + Z%*%d
  expect_equal(as.numeric(fast_nnls(ZTx = t(Z)%*%x,
            ZTZ = t(Z)%*%Z)),
  optim(c(0,0,0),
        function(d) sum((x - Z%*%matrix(d,ncol = 1))^2),
        method = "L-BFGS-B",
        lower = 0,
        upper = 1000)$par)

})


test_that("Another mixed example works", {
  set.seed(0)
  d <- matrix(c(-1,-2,3),ncol = 1)
  Z <- matrix(rnorm(300),ncol = 3)
  x <- rnorm(100) + Z%*%d
  expect_equal(as.numeric(fast_nnls(ZTx = t(Z)%*%x,
                         ZTZ = t(Z)%*%Z)),
               optim(c(0,0,0),
                     function(d) sum((x - Z%*%matrix(d,ncol = 1))^2),
                     method = "L-BFGS-B",
                     lower = 0,
                     upper = 1000)$par)

})

test_that("Mixed example with small positive solutions works", {
  set.seed(0)
  d <- matrix(c(0,0,4*1e-7,5*1e-9),ncol = 1)
  Z <- matrix(rnorm(400),ncol = 4)
  x <-  Z%*%d
  expect_equal(as.numeric(fast_nnls(ZTx = t(Z)%*%x,
                                    ZTZ = t(Z)%*%Z)),
              as.numeric(d))

})

test_that("Mixed example with many small positive solutions works", {
  set.seed(0)
  d <- matrix(c(rep(0,50),(1e-9)*sample(1:10000,50)),ncol = 1)
  Z <- matrix(rnorm(100^2),ncol = 100)
  x <-  Z%*%d
  expect_equal(as.numeric(fast_nnls(ZTx = t(Z)%*%x,
                                    ZTZ = t(Z)%*%Z)),
               as.numeric(d))

})


test_that("Entirely negative example works", {
  set.seed(0)
  d <- matrix(c(-1,-2,-3),ncol = 1)
  Z <- matrix(rnorm(300),ncol = 3)
  x <- rnorm(100) + Z%*%d
  expect_equal(as.numeric(fast_nnls(ZTx = t(Z)%*%x,
            ZTZ = t(Z)%*%Z)),

  optim(c(0,0,0),
        function(d) sum((x - Z%*%matrix(d,ncol = 1))^2),
        method = "L-BFGS-B",
        lower = 0,
        upper = 1000)$par)

})


test_that("One-dimensional example works", {
  set.seed(0)
  d <- matrix(c(-1),ncol = 1)
  Z <- matrix(rnorm(100),ncol = 1)
  x <- rnorm(100) + Z%*%d

  expect_equal(as.numeric(fast_nnls(ZTx = t(Z)%*%x,
            ZTZ = t(Z)%*%Z)),

  optim(c(0),
        function(d) sum((x - Z%*%matrix(d,ncol = 1))^2),
        method = "L-BFGS-B",
        lower = 0,
        upper = 1000)$par)

})

test_that("Higher-dimensional example works", {
  set.seed(0)
  d <- matrix(rnorm(50),ncol = 1)
  Z <- matrix(rnorm(5000),ncol = 50)
  x <- rnorm(100) + Z%*%d

  expect_equal(as.numeric(fast_nnls(ZTx = t(Z)%*%x,
                                    ZTZ = t(Z)%*%Z)),

               optim(rep(0,50),
                     function(d) sum((x - Z%*%matrix(d,ncol = 1))^2),
                     method = "L-BFGS-B",
                     lower = 0,
                     upper = 1000)$par,
               tolerance = 1e-4)

})

test_that("Giving starting active set doesn't break nnls", {
  set.seed(0)
  d <- matrix(rnorm(50),ncol = 1)
  Z <- matrix(rnorm(5000),ncol = 50)
  x <- rnorm(100) + Z%*%d
  active_start <- sample(1:50,10)

  expect_equal(as.numeric(fast_nnls(ZTx = t(Z)%*%x,
                                    ZTZ = t(Z)%*%Z,
                                    active = active_start)),

               optim(rep(0,50),
                     function(d) sum((x - Z%*%matrix(d,ncol = 1))^2),
                     method = "L-BFGS-B",
                     lower = 0,
                     upper = 1000)$par,
               tolerance = 1e-4)

})

test_that("Giving starting active set doesn't break nnls", {
  set.seed(0)
  d <- matrix(rnorm(50),ncol = 1)
  Z <- matrix(rnorm(5000),ncol = 50)
  x <- rnorm(100) + Z%*%d
  active_start <- sample(1:50,10)

  expect_equal(as.numeric(fast_nnls(ZTx = t(Z)%*%x,
                                    ZTZ = t(Z)%*%Z,
                                    active = active_start)),

               optim(rep(0,50),
                     function(d) sum((x - Z%*%matrix(d,ncol = 1))^2),
                     method = "L-BFGS-B",
                     lower = 0,
                     upper = 1000)$par,
               tolerance = 1e-4)

})

test_that("Underdetermined system doesn't stall", {
  set.seed(0)
  d <- 1:5
  Z <- cbind(matrix(rnorm(400),ncol = 4),0)
  x <- rnorm(100) + Z%*%d


  expect_equal(as.numeric(fast_nnls(ZTx = t(Z)%*%x,
                                    ZTZ = t(Z)%*%Z)),

               optim(rep(0,5),
                     function(d) sum((x - Z%*%matrix(d,ncol = 1))^2),
                     method = "L-BFGS-B",
                     lower = 0,
                     upper = 1000)$par,
               tolerance = 1e-4)

})


test_that("Can set some variables as unconstrained", {
  set.seed(0)
  d <- matrix(rnorm(50),ncol = 1)
  Z <- matrix(rnorm(5000),ncol = 50)
  x <- rnorm(100) + Z%*%d
  active_start <- sample(6:50,10)

  expect_equal(as.numeric(fast_nnls(ZTx = t(Z)%*%x,
                                    ZTZ = t(Z)%*%Z,
                                    active = active_start,
                                    unconstrained = 1:5)),

               optim(rep(0,50),
                     function(d) sum((x - Z%*%matrix(d,ncol = 1))^2),
                     method = "L-BFGS-B",
                     lower = c(rep(-1000,5),
                               rep(0,45)),
                     upper = 1000)$par,
               tolerance = 1e-4)

})

test_that("Large example works", {
  set.seed(0)
  d <- matrix(rnorm(500),ncol = 1)
  Z <- matrix(rnorm(500000),ncol = 500)
  x <- rnorm(1000) + Z%*%d

  expect_equal(as.numeric(fast_nnls(ZTx = t(Z)%*%x,
                                    ZTZ = t(Z)%*%Z)),

               optim(rep(0,500),
                     function(d) sum((x - Z%*%matrix(d,ncol = 1))^2),
                     method = "L-BFGS-B",
                     lower = 0,
                     upper = 1000)$par,
               tolerance = 1e-4)

})
