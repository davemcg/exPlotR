library(geyser)

context("Check that shiny app is generated")

load(system.file('extdata/tiny_rse.Rdata', package = 'geyser'))


test_that("Shiny app is generated", {
  expect_is(geyser(tiny_rse,'Test RSE'), "shiny.appobj")
})
