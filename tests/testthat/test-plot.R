library(ggplot2)
set.seed(12345)
n= 200
simu = simu_bfem(n, which = "Chang1983")
Y = simu$Y
p = ncol(Y)
K = 2
d= K - 1
res.bfem = bfem(Y, K, model="DB", init = 'kmeans', nstart = 2)

test_that("Plot bound", {
  expect_true(is.ggplot(plot.bfem(res.bfem, type = "elbo")))
})

test_that("Plot crit", {
  expect_true(is.ggplot(plot.bfem(res.bfem, type = "crit", crit = 'bic')))
  expect_true(is.ggplot(plot.bfem(res.bfem, type = "crit", crit = 'icl')))
  expect_true(is.ggplot(plot.bfem(res.bfem, type = "crit", crit = 'aic')))
})