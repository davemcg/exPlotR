
test_that("Shiny app is generated", {
  expect_is(geyser(tiny_rse,'Test RSE'), "shiny.appobj")
})
