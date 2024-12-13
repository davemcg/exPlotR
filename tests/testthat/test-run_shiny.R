
test_that("Shiny app is generated", {
  load(system.file('extdata/tiny_rse.Rdata', package = 'geyser'))
  expect_type(geyser(tiny_rse,'Test RSE'), "list")
})
