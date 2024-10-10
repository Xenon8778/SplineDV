test_that("HVG works with matrix",{
  X <- abs(matrix(round(rnorm(2000*500,0,10),0), nrow=2000, ncol=500))
  rownames(X) <- paste0('g', as.character(1:2000))
  res <- splineHVG(X)
  expect_s3_class(res, "data.frame")
  expect_identical(ncol(res), 14L)
})
