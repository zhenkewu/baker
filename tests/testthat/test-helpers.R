context("simple helpers")

test_that("utilities for checking lists",{
  l = list(a = 5, b = 1:2)
  expect_false(is_length_all_one(l))
  l = list(a = 5, b = 1)
  expect_true(is_length_all_one(l))
})
