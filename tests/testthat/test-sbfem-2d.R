set.seed(12345)
n = 900
simu = simu_bfem(n, which = "section4.2", p = 50, noise = 1e-2)
Y = simu$Y %*% simu$W
p = ncol(Y)
K = 3
d = K - 1

# method = 'reg'
# res.fem = fem(Y, K, model='DB', init = 'kmeans', nstart = 10, method = method)
# fem.ari(res.fem, simu$cls)
# 
# res.sfem = sfem(Y, K,obj = res.fem)
# plot(res.sfem, 3)
# fem.ari(res.sfem, simu$cls)
# 
# res.bfem = bfem(Y, K, model = 'DB',  init = 'kmeans', nstart = 10, method = "gs")
# fem.ari(res.bfem, simu$cls)
# 
# res.sbfem = sbfem(Y, K, obj = res.bfem, maxit.sparse = 1)
# fem.ari(res.sbfem, simu$cls)
# plot(res.sbfem)

test_that("bfem kmeans init", {
  res.sbfem = sbfem(Y, K, model="DB", init = 'kmeans', nstart = 10, maxit.sparse = 1)
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
