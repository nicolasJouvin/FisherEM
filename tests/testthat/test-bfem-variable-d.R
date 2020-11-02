set.seed(12345)
n = 100
simu = simu_bfem(n, which = "Chang1983")
Y = simu$Y
p = ncol(Y)
K = 2
d = 4

test_that("bfem method svd", {
  res.bfem = bfem(Y, K, d, model="DB", init = 'kmeans', method = 'svd', nstart = 2, mc.cores = 2)
  expect_equal(res.bfem$K, K)
  expect_equal(res.bfem$d, d)
  expect_equal(dim(res.bfem$P), c(n, K))
  expect_equal(dim(res.bfem$var_param$Varmeank), c(d, K))
  expect_equal(dim(res.bfem$var_param$Varcovk), c(d, d, K))
  expect_equal(dim(res.bfem$U), c(p, d))
  expect_equal(dim(res.bfem$proj), c(n, d))
  expect_equal(length(res.bfem$elbos), res.bfem$n_ite+1)
})

test_that("bfem method gs", {
  res.bfem = bfem(Y, K, d, model="DB", init = 'kmeans', method = 'gs', nstart = 2, mc.cores = 2)
  expect_equal(res.bfem$K, K)
  expect_equal(res.bfem$d, d)
  expect_equal(dim(res.bfem$P), c(n, K))
  expect_equal(dim(res.bfem$var_param$Varmeank), c(d, K))
  expect_equal(dim(res.bfem$var_param$Varcovk), c(d, d, K))
  expect_equal(dim(res.bfem$U), c(p, d))
  expect_equal(dim(res.bfem$proj), c(n, d))
  expect_equal(length(res.bfem$elbos), res.bfem$n_ite+1)
})

test_that("bfem method reg", {
  res.bfem = bfem(Y, K, d, model="DB", init = 'kmeans', method = 'reg', nstart = 2, mc.cores = 2)
  expect_equal(res.bfem$K, K)
  expect_equal(res.bfem$d, d)
  expect_equal(dim(res.bfem$P), c(n, K))
  expect_equal(dim(res.bfem$var_param$Varmeank), c(d, K))
  expect_equal(dim(res.bfem$var_param$Varcovk), c(d, d, K))
  expect_equal(dim(res.bfem$U), c(p, d))
  expect_equal(dim(res.bfem$proj), c(n, d))
  expect_equal(length(res.bfem$elbos), res.bfem$n_ite+1)
})

test_that("test all models", {
  res.bfem = bfem(Y, K, d, model="all", init = 'kmeans', nstart = 1, mc.cores = 2)
  expect_equal(res.bfem$K, K)
  expect_equal(res.bfem$d, d)
  expect_equal(dim(res.bfem$P), c(n, K))
  expect_equal(dim(res.bfem$var_param$Varmeank), c(d, K))
  expect_equal(dim(res.bfem$var_param$Varcovk), c(d, d, K))
  expect_equal(dim(res.bfem$U), c(p, d))
  expect_equal(dim(res.bfem$proj), c(n, d))
  expect_equal(length(res.bfem$elbos), res.bfem$n_ite+1)
})

test_that("test 3 models, K grid", {
  res.bfem = bfem(Y, K = 2:6, d=d, model=c('DkBk', 'AkjBk', 'AB'), init = 'kmeans', nstart = 1, 
                  maxit.em = 10, eps.em = 1e-3, maxit.ve = 3, mc.cores = 2)
  expect_equal(dim(res.bfem$var_param$Varmeank), c(res.bfem$d, res.bfem$K))
  expect_equal(dim(res.bfem$var_param$Varcovk), c(res.bfem$d, res.bfem$d, res.bfem$K))
  expect_equal(dim(res.bfem$U), c(p, res.bfem$d))
  expect_equal(dim(res.bfem$proj), c(n, res.bfem$d))
  expect_equal(length(res.bfem$elbos), res.bfem$n_ite+1)
})


Y = iris[,-5]
n = nrow(Y)
p = ncol(Y)
K = 3
d = 3

test_that("bfem method gs", {
  res.bfem = bfem(Y, K, d, model="DB", init = 'kmeans', method = 'gs', nstart = 2, mc.cores = 2)
  expect_equal(res.bfem$K, K)
  expect_equal(res.bfem$d, d)
  expect_equal(dim(res.bfem$P), c(n, K))
  expect_equal(dim(res.bfem$var_param$Varmeank), c(d, K))
  expect_equal(dim(res.bfem$var_param$Varcovk), c(d, d, K))
  expect_equal(dim(res.bfem$U), c(p, d))
  expect_equal(dim(res.bfem$proj), c(n, d))
  expect_equal(length(res.bfem$elbos), res.bfem$n_ite+1)
})

test_that("test 3 models, K grid", {
  res.bfem = bfem(Y, K = 2:6, d=d, model=c('DkBk', 'AkjBk', 'AB'), init = 'kmeans', nstart = 1, 
                  maxit.em = 10, eps.em = 1e-3, maxit.ve = 3, mc.cores = 2)
  expect_equal(dim(res.bfem$var_param$Varmeank), c(res.bfem$d, res.bfem$K))
  expect_equal(dim(res.bfem$var_param$Varcovk), c(res.bfem$d, res.bfem$d, res.bfem$K))
  expect_equal(dim(res.bfem$U), c(p, res.bfem$d))
  expect_equal(dim(res.bfem$proj), c(n, res.bfem$d))
  expect_equal(length(res.bfem$elbos), res.bfem$n_ite+1)
})