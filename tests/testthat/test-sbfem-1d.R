set.seed(12345)
n = 300
simu = simu_bfem(n, which = "Chang1983")
Y = simu$Y
p = ncol(Y)
K = 2
d = K - 1

test_that("bfem kmeans init", {
  res.sbfem = sbfem(Y, K, model="DB", init = 'kmeans', nstart = 10)
  expect_equal(res.sbfem$K, K)
  expect_equal(res.sbfem$d, d)
  expect_equal(dim(res.sbfem$P), c(n, K))
  expect_equal(dim(res.sbfem$var_param$Varmeank), c(d, K))
  expect_equal(dim(res.sbfem$var_param$Varcovk), c(d, d, K))
  expect_equal(dim(res.sbfem$U), c(p, d))
  expect_equal(dim(res.sbfem$proj), c(n, d))
  expect_equal(length(res.sbfem$elbos), res.sbfem$n_ite+1)
})
# 
# test_that("bfem random init", {
#   res.sbfem = bfem(Y, K, model="DB", init = 'random', nstart = 2, mc.cores = 2)
#   expect_equal(res.sbfem$K, K)
#   expect_equal(res.sbfem$d, d)
#   expect_equal(dim(res.sbfem$P), c(n, K))
#   expect_equal(dim(res.sbfem$var_param$Varmeank), c(d, K))
#   expect_equal(dim(res.sbfem$var_param$Varcovk), c(d, d, K))
#   expect_equal(dim(res.sbfem$U), c(p, d))
#   expect_equal(dim(res.sbfem$proj), c(n, d))
#   expect_equal(length(res.sbfem$elbos), res.sbfem$n_ite+1)
# })
# 
# test_that("bfem user init", {
#   Tinit = t(rmultinom(n, 1, prob = rep(1/2,2)))
#   res.sbfem = bfem(Y, K, model="DB", init = 'user', Tinit = Tinit, mc.cores = 2)
#   expect_equal(res.sbfem$K, K)
#   expect_equal(res.sbfem$d, d)
#   expect_equal(dim(res.sbfem$P), c(n, K))
#   expect_equal(dim(res.sbfem$var_param$Varmeank), c(d, K))
#   expect_equal(dim(res.sbfem$var_param$Varcovk), c(d, d, K))
#   expect_equal(dim(res.sbfem$U), c(p, d))
#   expect_equal(dim(res.sbfem$proj), c(n, d))
#   expect_equal(length(res.sbfem$elbos), res.sbfem$n_ite+1)
# })
