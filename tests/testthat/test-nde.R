context("rNDE bias corrections")

test_that("approximately no correction when bias is truly 0", {
  set.seed(123)
  u_ei <- 0
  am_intx <- 0
  dl <- simulate_data(n = 50000, 
                      u_ei = u_ei, am_intx = am_intx, 
                      yu_strength = 0, mu_strength = 0)
  small_dl <- simulate_data(n = 50000, 
                            u_ei = u_ei, am_intx = am_intx, 
                            yu_strength = 0, mu_strength = 0)
  
  uc <- mean(run_uncorrected(dl$df, am_intx, mean = FALSE))
  dg <- run_delta_gamma(dl$df, small_df = small_dl$df)
  ix <- run_intx_corr(dl$df, small_df = small_dl$df, am_intx = am_intx)
  
  # Less than 1 percent deviation from truth
  expect_true(abs(uc - dg)/uc < 0.01, label = "delta-gamma applies no correction")
  expect_true(abs(uc - ix)/uc < 0.01, label = "intx applies no correction")

})



context("true rNDE calculations")

test_that("rNDE is truly zero when no direct effect", {
  params <- return_dgp_parameters(u_ei = 0, am_intx = 0, 
                                  yu_strength = 0, mu_strength = 0)
  params$alpha["a"] <- 0
  nder <- calculate_nder(params = NULL, 
                         alpha = params$alpha,
                         beta = params$beta,
                         gamma = params$gamma,
                         u_ei = 0, am_intx = 0,
                         z1 = 0:1, z2 = 0:1,
                         yu_strength = 0, mu_strength = 0)
  expect_equal(nder, rep(0, 4))
})

test_that("rNDE non-zero even when only A effect on Y is interaction", {
  params <- return_dgp_parameters(u_ei = 0, am_intx = 1, 
                                  yu_strength = 0, mu_strength = 0)
  params$alpha["a"] <- 0
  nder <- calculate_nder(params = NULL, 
                         alpha = params$alpha,
                         beta = params$beta,
                         gamma = params$gamma,
                         u_ei = 0, am_intx = 1,
                         z1 = 0:1, z2 = 0:1,
                         yu_strength = 0, mu_strength = 0)
  expect_gt(min(abs(nder)), 0)
})

