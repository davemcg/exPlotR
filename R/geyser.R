#' @title geyser
#'
#' @description Run shiny app to use SummarizedExperiment object to display genomics data
#'
#' @export
#' 
#' @import bslib
#' @importFrom dplyr any_of select pull row_number
#' @importFrom ComplexHeatmap draw
#' @import htmltools
#' @import SummarizedExperiment
#' @importFrom magrittr "%>%"
#' @importFrom tibble rownames_to_column
#' @import shiny
#' 
#' @param rse SummarizedExperiment object
#' @param app_name Title name that goes on the top left of the Shiny app
#' @param primary_color The title bar color
#' @param secondary_color The plot action button color
#'
#' @details
#'
#' Shiny app uses the rowData rownames to define the genes. The colData field is made
#' fully available to make custom plot groupings.
#'
#' @author David McGaughey
#'
#' @examples
#'
#' if (interactive()){
#'   load(system.file('extdata/tiny_rse.Rdata', package = 'geyser'))
#'   geyser(tiny_rse)
#' }
#'

geyser <- function(rse, 
                   app_name = "geyser",
                   primary_color = "#3A5836",
                   secondary_color = "#d5673e") {
  ui <- page_navbar(
    title = app_name,
    theme = theme_ui(primary_color = primary_color, 
                     secondary_color = secondary_color),
    tags$style(HTML('table.dataTable tr.active td, table.dataTable tr.active 
                    {background-color: #3A5836 !important;}')),
    tags$style(HTML('table.dataTable tr.selected td, table.dataTable td.selected 
                    {background-color: pink !important;}')),
    selected = "Full Sample Metadata",
    collapsible = TRUE,
    nav_panel(
      title = "Full Sample Metadata",
      card(
        full_screen = TRUE,
        card_header("colData", class = 'bg-dark'),
        card_body(
          DT::dataTableOutput("table_full", width = "85%",fill = FALSE)
        )
      )
    ),
    nav_panel(
      title = "Plotting",
      page_fluid(
        layout_sidebar(
          height = '100%',
          sidebar = sidebar(
            title = 'Plot Parameters',
            open = TRUE,
            accordion(
              multiple = TRUE,
              accordion_panel(
                "Grouping and Features",
                selectizeInput("groupings",
                               "Sample Grouping(s):",
                               choices = NULL,
                               multiple = TRUE,
                ),
                selectizeInput("genes",
                               "Assay Feature(s): ",
                               choices = NULL,
                               multiple = TRUE),
                selectizeInput("slot",
                               "Assay Type:",
                               choices = NULL,
                               multiple = FALSE
                ),
                layout_column_wrap(width = 0.5, 
                                   checkboxInput("expression_scale", 
                                                 label = 'log2(expression)', 
                                                 value = TRUE)),
              ),
              accordion_panel("Sample Filtering",
                              card(full_screen = TRUE,
                                   DT::dataTableOutput("table",width = "105%",
                                                       fill = FALSE),
                                   actionButton('clear_colData_row_selections', 
                                                'Clear Rows'))
              )
            ),
            width = '40%'),
          navset_card_tab(
            full_screen = TRUE,
            nav_panel("Box Plot",
                      actionButton('exp_plot_button','Draw Box Plot'),
                      plotOutput("exp_plot",height = '100%')
            ),
            nav_panel("Heatmap",
                      layout_columns(fill = FALSE,
                                     checkboxInput("col_clust", 
                                                   label = "Cluster Columns", 
                                                   value = TRUE),
                                     checkboxInput("row_clust", 
                                                   label = "Cluster Rows", 
                                                   value = TRUE)),
                      actionButton('hm_plot_button','Draw Heatmap'),
                      plotOutput("hm_plot",height = '100%')
            )
          )
        )
      )
    ),
    nav_spacer(),
    nav_menu(
      title = "Info",
      align = "right",
      quick_start_ui(),
      #overview_ui(),
      nav_item(tags$a("Code Source (external link)", 
                      href = "https://github.com/davemcg/geyser", 
                      target = "_blank"))
    )
  )
  
  # this argument yanked via the R/geyser.R function
  rse_name <- deparse(substitute(rse))
  
  server <- function(input, output, session) {
    # select sample columns to group on -----
    updateSelectizeInput(session, 'groupings',
                         choices = colnames(colData(get(rse_name))) %>% 
                           sort(),
                         selected = '',
                         server = TRUE)
    # select genes ----
    updateSelectizeInput(session, 'genes',
                         choices = rownames(get(rse_name)) %>% sort(),
                         selected = '',
                         server = TRUE)
    # select assay slot (usually counts) ----
    updateSelectizeInput(session, 'slot',
                         choices = assays(get(rse_name)) %>% 
                           names() %>% sort(),
                         selected = if ('counts' %in% 
                                        names(assays(get(rse_name)))){'counts'} 
                         else {names(assays(get(rse_name)))[1]},
                         server = TRUE)
    # expression plot ----
    # R/exp_plot.R
    exp_plot_reactive <- eventReactive(input$exp_plot_button, {
      .exp_plot(input, get(rse_name))
    })
    output$exp_plot <- renderPlot({
      exp_plot_reactive()$plot},
      height = eventReactive(input$exp_plot_button,
                             {max(600, 30 * length(input$genes) * exp_plot_reactive()$grouping_length)})
    )
    # hm plot -----
    # R/heatmap.R
    hm_plot_reactive <- eventReactive(input$hm_plot_button, {
      .hm_plot(input, get(rse_name))
    })
    output$hm_plot <- renderPlot({
      draw(hm_plot_reactive()$plot)},
      height = eventReactive(input$hm_plot_button,
                             {max(400, 0.7 * hm_plot_reactive()$grouping_length)})
    )
    
    # sample data table -----
    output$table <- DT::renderDataTable(
      colData(get(rse_name)) %>%
        data.frame() %>% 
        rownames_to_column('rse_sample_id') %>% 
        select('rse_sample_id', any_of(input$groupings)) %>%
        DT::datatable(rownames= FALSE,
                      options = list(autoWidth = TRUE,
                                     pageLength = 15,
                                     dom = 'tp'),
                      filter = list(position = 'top', clear = FALSE)),
      server = TRUE
    )
    ## proxy to clear row selection -----
    ## https://yihui.shinyapps.io/DT-proxy/
    proxy <- DT::dataTableProxy('table')
    observeEvent(input$clear_colData_row_selections, {
      proxy %>% DT::selectRows(NULL)
    })
    
    # sample data table full -----
    output$table_full <- DT::renderDataTable(
      colData(get(rse_name)) %>%
        data.frame() %>% 
        rownames_to_column('rse_sample_id') %>% 
        DT::datatable(rownames= FALSE,
                      options = list(autoWidth = TRUE,
                                     pageLength = 25),
                      filter = list(position = 'top', clear = FALSE),
                      selection = 'none'
        ),
      server = TRUE
    )
    session$onSessionEnded(function() {
      stopApp()
    })
  }
  
  server_env <- environment(server)
  
  # variables for server.R
  server_env$rse <- rse
  server_env$app_name <- app_name
  server_env$primary_color <- primary_color
  server_env$secondary_coor <- secondary_color
  
  app <- shinyApp(ui, server)
  return(app)
}
