
test_that("Exp Plot is generated without error", {
  input <- list()
  input$genes <- c("TYRP1 (ENSG00000107165.12)","OPN1LW (ENSG00000102076.9)")
  input$groupings <- c('disease')
  input$slot <- 'counts'
  input$expression_scale <- TRUE
  expect_error(geyser:::.exp_plot(input, 'tiny_rse', 'counts')$plot, NA)
})

test_that("Heatmap is generated without error", {
  input <- list()
  input$genes <- c("TYRP1 (ENSG00000107165.12)","OPN1LW (ENSG00000102076.9)")
  input$groupings <- c('disease')
  input$slot <- 'counts'
  input$expression_scale <- TRUE
  input$row_clust <- TRUE
  input$col_clust <- TRUE
  expect_error(geyser:::.hm_plot(input, 'tiny_rse', 'counts')$plot, NA)
})