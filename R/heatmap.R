#' @title hm_plot
#'
#' @description draws the expression heatmap
#'
#' @keywords internal
#'
#' @import SummarizedExperiment
#' @import ComplexHeatmap
#' @import dplyr
#' @import tibble
#' @importFrom tidyr pivot_longer pivot_wider
#'
#'
#' @param input From ui.R
#' @param rse_name Name of the rse object
#' @param slot which slot to pull the count data from the rse assay
#'
#' @details
#'
#' Makes the box plot for the geyser Shiny app
#'
#' @author David McGaughey
#'
#' @examples
#'
#' \dontrun{
#' exp_plot(input, rse_name, 'counts')
#' }
#'

hm_plot <- function(input, rse_name, slot){
  Gene <- rowid <- sample_unique_id <- counts <- group <- NULL
  genes <- input$genes
  groupings <- input$groupings

  if (length(genes) < 1 || length(groupings) < 1){
    showModal(modalDialog(title = "Heatmap Error",
                          "Have you specified at least one grouping and one gene?",
                          easyClose = T,
                          footer = NULL))
    stop()
  }

  # pull gene counts and left_join with colData
  pdata <- assay(get(rse_name),  input$slot)[genes, ,drop = FALSE] %>%
    data.frame() %>% 
    tibble::rownames_to_column('Gene') %>% 
    pivot_longer(-Gene, values_to = 'counts', names_to = 'sample_unique_id') %>%
    left_join(colData(get(rse_name)) %>%
                data.frame() %>% 
                tibble::rownames_to_column('sample_unique_id') %>% 
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
    lab_text <- "scale(log2(c))"
  } else {
    lab_text <- "scale(c)"
  }
  output <- list()
  pfdata <- pfdata %>%
    # make custom column with user selected groupings of columns
    tidyr::unite("group", all_of(groupings), remove = FALSE, sep = " | ")
  # make df for ComplexHeatmap
  pfdf <- pfdata %>%
    dplyr::select(Gene, sample_unique_id, counts) %>%
    pivot_wider(names_from = sample_unique_id, values_from = counts) 
  col_labels <- colnames(pfdf)[-1] # yank out col names before data.frame conversion to preserve special characters
  pfdf <- data.frame(pfdf)
  row.names(pfdf) <- pfdf$Gene
  pfdf <- pfdf[,-1]
  hm_data <- t(scale(t(pfdf[,pfdata$sample_unique_id %>% unique()])))
  # row clustering fails if there are zero count rows
  row_clustering <- input$row_clust
  if (min(rowSums(hm_data)) == 0){
    row_clustering <- FALSE
  }
  output$plot <- Heatmap(hm_data,
                         column_split = pfdata %>% filter(Gene == genes[1]) %>% pull(group),
                         column_title_rot = 90,
                         column_labels = col_labels,
                         cluster_columns = input$col_clust,
                         cluster_rows = row_clustering,
                         name = lab_text)
  output$grouping_length <- nrow(pfdata) + (nchar(pfdata$group) %>% max()) + (nchar(pfdf) %>% max())
  output
}
