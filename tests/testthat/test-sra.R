test_that("Cimple computation works", {

  res1 <-sra(matrix(cbind(1:8,c(1,2,3,5,6,7,4,8),c(1,5,3,4,2,8,7,6)),ncol=3))

  expect_equal(length(res1), 8)
  expect_equal(as.numeric(res1[1]), 0)
  expect_equal(as.numeric(res1[2]), 4/3)
  expect_equal(as.numeric(res1[5]), 4/3)	
})
