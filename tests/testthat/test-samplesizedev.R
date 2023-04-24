test_that("Error message for missing prevalence", {
  expect_error(samplesizedev(S=0.9, c = 0.85, n.predictors = 10,  nsim = 100, parallel=FALSE))
})

test_that("Error message for missing C-statistic", {
  expect_error(samplesizedev(S=0.9, p = 0.2, n.predictors = 10,  nsim = 100, parallel=FALSE))
})

test_that("Error message for missing numbers of predictors", {
  expect_error(samplesizedev(S=0.9, c = 0.85, p = 0.2,  nsim = 100, parallel=FALSE))
})


# test_that("Sample size calculations are correct", {
#   a <- sampsizeval(p=0.057, c=0.77, se_c=0.025, se_cs =0.15, se_cl = 0.15)
#   expect_equal(a$size_c_statistic,1599)
#   expect_equal(a$size_calibration_slope,934)
#   expect_equal(a$size_calibration_large,897)
#   expect_equal(a$size_recommended,1599)
# })









