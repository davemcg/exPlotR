

test_that("Exp Plot is generated without error", {
  load(system.file('extdata/tiny_rse.Rdata', package = 'geyser'))
  
  input <- list()
  input$genes <- c("TYRP1 (ENSG00000107165.12)","OPN1LW (ENSG00000102076.9)")
  input$groupings <- c('disease')
  input$slot <- 'counts'
  input$expression_scale <- TRUE
  input$color_by <- 'tissue'
  plot <- geyser:::.exp_plot(input, tiny_rse, 'counts')$plot
  expect_s3_class(plot, 'gg')
})

test_that("Heatmap is generated without error", {
  load(system.file('extdata/tiny_rse.Rdata', package = 'geyser'))
  
  input <- list()
  input$genes <- c("TYRP1 (ENSG00000107165.12)","OPN1LW (ENSG00000102076.9)")
  input$groupings <- c('disease')
  input$slot <- 'counts'
  input$expression_scale <- TRUE
  input$row_clust <- TRUE
  input$col_clust <- TRUE
  plot <- geyser:::.hm_plot(input, tiny_rse, 'counts')$plot
  expect_s4_class(plot, "Heatmap")
})
