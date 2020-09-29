set.seed(12345)

n = 100
simu = simu_bfem(n, which = "section4.3", snr = 10)
Y = simu$Y
p = ncol(Y)
K = 3
d = K - 1

test_that("bfem with n>>p and method = 'gs' returns an error", {
  expect_error(res.bfem <- bfem(Y, K, model="DB", init = 'kmeans', nstart = 2, mc.cores = 2,
                  kernel = 'linear'))
})

# ------------ Linear kernel
test_that("linear kernel + svd", {
  res.bfem <- bfem(Y, K, model="DB", init = 'kmeans', nstart = 2, mc.cores = 2,
                                method = 'svd', kernel = 'linear')
  expect_equal(res.bfem$K, K)
  expect_equal(res.bfem$d, d)
  expect_equal(dim(res.bfem$P), c(n, K))
  expect_equal(dim(res.bfem$var_param$Varmeank), c(d, K))
  expect_equal(dim(res.bfem$var_param$Varcovk), c(d, d, K))
  expect_equal(dim(res.bfem$U), c(p, d))
  expect_equal(dim(res.bfem$proj), c(n, d))
  expect_equal(length(res.bfem$elbos), res.bfem$n_ite+1)
})

test_that("linear kernel + reg", {
  res.bfem <- bfem(Y, K, model="DB", init = 'kmeans', nstart = 2, mc.cores = 2,
                   method = 'reg', kernel = 'linear')
  expect_equal(res.bfem$K, K)
  expect_equal(res.bfem$d, d)
  expect_equal(dim(res.bfem$P), c(n, K))
  expect_equal(dim(res.bfem$var_param$Varmeank), c(d, K))
  expect_equal(dim(res.bfem$var_param$Varcovk), c(d, d, K))
  expect_equal(dim(res.bfem$U), c(p, d))
  expect_equal(dim(res.bfem$proj), c(n, d))
  expect_equal(length(res.bfem$elbos), res.bfem$n_ite+1)
})

# ------------ Sigmoid kernel
test_that("sigmoid kernel + svd", {
  res.bfem <- bfem(Y, K, model="DB", init = 'kmeans', nstart = 2, mc.cores = 2,
                   method = 'svd', kernel = 'sigmoid')
  expect_equal(res.bfem$K, K)
  expect_equal(res.bfem$d, d)
  expect_equal(dim(res.bfem$P), c(n, K))
  expect_equal(dim(res.bfem$var_param$Varmeank), c(d, K))
  expect_equal(dim(res.bfem$var_param$Varcovk), c(d, d, K))
  expect_equal(dim(res.bfem$U), c(p, d))
  expect_equal(dim(res.bfem$proj), c(n, d))
  expect_equal(length(res.bfem$elbos), res.bfem$n_ite+1)
})

test_that("sigmoid kernel + reg", {
  res.bfem <- bfem(Y, K, model="DB", init = 'kmeans', nstart = 2, mc.cores = 2,
                   method = 'reg', kernel = 'sigmoid')
  expect_equal(res.bfem$K, K)
  expect_equal(res.bfem$d, d)
  expect_equal(dim(res.bfem$P), c(n, K))
  expect_equal(dim(res.bfem$var_param$Varmeank), c(d, K))
  expect_equal(dim(res.bfem$var_param$Varcovk), c(d, d, K))
  expect_equal(dim(res.bfem$U), c(p, d))
  expect_equal(dim(res.bfem$proj), c(n, d))
  expect_equal(length(res.bfem$elbos), res.bfem$n_ite+1)
})

# ------------ RBF kernel
test_that("RBF kernel + svd", {
  res.bfem <- bfem(Y, K, model="DB", init = 'kmeans', nstart = 2, mc.cores = 2,
                   method = 'svd', kernel = 'rbf')
  expect_equal(res.bfem$K, K)
  expect_equal(res.bfem$d, d)
  expect_equal(dim(res.bfem$P), c(n, K))
  expect_equal(dim(res.bfem$var_param$Varmeank), c(d, K))
  expect_equal(dim(res.bfem$var_param$Varcovk), c(d, d, K))
  expect_equal(dim(res.bfem$U), c(p, d))
  expect_equal(dim(res.bfem$proj), c(n, d))
  expect_equal(length(res.bfem$elbos), res.bfem$n_ite+1)
})

test_that("RBF kernel + reg", {
  res.bfem <- bfem(Y, K, model="DB", init = 'kmeans', nstart = 2, mc.cores = 2,
                   method = 'reg', kernel = 'rbf')
  expect_equal(res.bfem$K, K)
  expect_equal(res.bfem$d, d)
  expect_equal(dim(res.bfem$P), c(n, K))
  expect_equal(dim(res.bfem$var_param$Varmeank), c(d, K))
  expect_equal(dim(res.bfem$var_param$Varcovk), c(d, d, K))
  expect_equal(dim(res.bfem$U), c(p, d))
  expect_equal(dim(res.bfem$proj), c(n, d))
  expect_equal(length(res.bfem$elbos), res.bfem$n_ite+1)
})

# fem.ari(res.bfem, simu$cls)
# plot(res.bfem)
