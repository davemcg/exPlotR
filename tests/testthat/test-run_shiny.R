
test_that("Shiny app is generated", {
  load(system.file('extdata/tiny_rse.Rdata', package = 'geyser'))
  expect_is(geyser(tiny_rse,'Test RSE'), "shiny.appobj")
})
