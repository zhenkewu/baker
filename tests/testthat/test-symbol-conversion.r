context("symbol conversions")

test_that("nplcm utilities for converting symbols
          to 0/1 coding",{
           expect_equal(symb2I("A",c("A","B","C")),matrix(c(1,0,0),nrow=1))
           expect_equal(symb2I("A+B",c("A","B","C")),matrix(c(1,1,0),nrow=1))
           expect_equal(symb2I("NoA",c("A","B","C")),matrix(c(0,0,0),nrow=1))
           expect_equal(symb2I(c("A","B+C"),c("A","B","C")),matrix(c(1,0,0,0,1,1),
                                                                   nrow=2,
                                                                   byrow=TRUE))
          })