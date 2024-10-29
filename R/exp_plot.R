#' @title exp_plot
#'
#' @description draws the expression box plot
#'
#' @keywords internal
#'
#' @import SummarizedExperiment
#' @import ggplot2
#' @importFrom ggbeeswarm geom_beeswarm
#' @importFrom dplyr filter left_join mutate pull row_number 
#' @importFrom magrittr "%>%"
#' @importFrom tidyr all_of pivot_longer pivot_wider unite
#' @importFrom tibble rownames_to_column
#' 
#' @param input From ui.R
#' @param rse The rse object
#' @param slot which slot to pull the count data from the rse assay
#'
#' @details
#'
#' Makes the box plot for the geyser Shiny app
#'
#' @author David McGaughey
#'
#' @returns 
#' 
#' Returns a list with the $plot slot holding ggplot object and $grouping_length contains
#' the number of features to scale the plot
#' 
#' @examples
#'
#' load(system.file('extdata/tiny_rse.Rdata', package = 'geyser'))
#' input <- list()
#' input$genes <- c("TYRP1 (ENSG00000107165.12)","OPN1LW (ENSG00000102076.9)")
#' input$groupings <- c('disease')
#' input$slot <- 'counts'
#' input$expression_scale <- TRUE
#' input$color_by <- 'tissue'
#' geyser:::.exp_plot(input, tiny_rse, 'counts')$plot

.exp_plot <- function(input, rse, slot){
  Gene <- rowid <- group <- counts <- geyser_color_by <- NULL
  
  genes <- input$genes
  groupings <- input$groupings
  
  if (length(genes) < 1 || length(groupings) < 1){
    showModal(modalDialog(title = "Box Plot Error",
                          "Have you specified at least one grouping and one gene?",
                          easyClose = TRUE,
                          footer = NULL))
    stop()
  }
  
  # pull gene counts and left_join with colData
  pdata <- assay((rse), input$slot)[genes, ,drop = FALSE] %>%
    data.frame() %>% 
    rownames_to_column('Gene') %>% 
    pivot_longer(-Gene, values_to = 'counts', names_to = 'sample_unique_id') %>%
    left_join(colData((rse)) %>%
                data.frame() %>% 
                rownames_to_column('sample_unique_id') %>% 
                mutate(rowid = row_number()),
              by = 'sample_unique_id')
  if (length(input$table_rows_selected)){
    pfdata <- pdata %>% filter(rowid %in% input$table_rows_selected)
  } else {
    pfdata <- pdata
  }
  # optional (but set as default) log2 scaling
  if (input$expression_scale){
    pfdata$counts <- log2(pfdata$counts + 1)
    ylab_text <- paste0("log2(", input$slot, ")")
  } else {
    ylab_text <- input$slot
  }
  output <- list()
  # if only one grouping given then just column as is to retain potential factor
  if (length(groupings) == 1){
    pfdata$geyser_group <- pfdata[,groupings] %>% pull(1)
  } else {
    pfdata <- pfdata %>%
      # make custom column with user selected groupings of columns
      unite("geyser_group", all_of(groupings), remove = FALSE, sep = " | ")
  }
  
  if (input$color_by != ''){
    pfdata$geyser_color_by <- pfdata[,input$color_by] %>% pull(1)
    output$plot <- pfdata %>%
      ggplot(aes(x=geyser_group,y=counts, color = geyser_color_by, group = geyser_group)) +
      geom_boxplot() +
      geom_beeswarm(dodge.width = 0.75) +
      coord_flip() +
      xlab(paste0(groupings, collapse = ' | ')) +
      ylab(ylab_text) +
      theme_linedraw(base_size = 16) +
      facet_wrap(~Gene, ncol = 1) +
      guides(col= guide_legend(title= input$color_by))
  } else {
    output$plot <- pfdata %>%
      ggplot(aes(x=geyser_group,y=counts, group = geyser_group)) +
      geom_boxplot() +
      geom_beeswarm() +
      coord_flip() +
      xlab(paste0(groupings, collapse = ' | ')) +
      ylab(ylab_text) +
      theme_linedraw(base_size = 16) +
      facet_wrap(~Gene, ncol = 1)
  }
  output$grouping_length <- pfdata$geyser_group %>% unique() %>% length()
  output
}
