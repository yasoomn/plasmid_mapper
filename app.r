library(shiny)
library(bslib)
library(plasmapR)
library(ggplot2)
library(dplyr)
library(DT)
library(stringr)

remove_junk_features <- function(df) {
  
  junk_features.name <- c("Geneious type: ligation", "Geneious type: Fragment")
  junk_features.type <- c("source")
  
  df <- df |> filter(!note %in% junk_features.name, !type %in% junk_features.type)


  return(df)
}

fix_geneious_type <- function(df) {
  if ("note" %in% colnames(df)) {
    df$type[df$note == "Geneious type: promoter"] <- "Promoter"
    df$type[df$note == "Geneious type: terminator"] <- "terminator"
    df$type[df$note == "Geneious type: origin of replication"] <- "rep_origin"
  }
  
  return(df)
}

standardize_ori_position <- function(df) {
  # ORI should start at the position 0
  ori_index <- which(df$type == "rep_origin")
}

Geneious_palette = c("CDS" = "#FFFF00", 
  "gene" = "#00B200", 
  "misc_feature" = "#AAAAAA", 
  "Promoter" = "#BAFF00", 
  "promoter" = "#BAFF00",
  "terminator" = "#FF5300", 
  "rRNA" = "#FF00FF", 
  "tRNA" = "#00FFFF", 
  "rep_origin" = "#50B3FF",
  "5'UTR" = "#F5F5F5", 
  "3'UTR" = "#F5F5F5",
  "intron" = "#F5F5F5", 
  "regulatory" = "#0EA089")

# Define UI for app that draws a histogram ----
ui <- page_sidebar(
  # App title ----
  title = "Plasmid Viewer",
  # Sidebar panel for inputs ----
  sidebar = sidebar(
    # Input: Slider for the number of bins ----
    fileInput("genbank", "Choose a Genbank File"),
    checkboxInput("remove_junk", "Remove Junk Features", value = TRUE), 
    checkboxInput("invert_seq", "Invert sequence", value = FALSE), 
    checkboxInput("circular", "Circular view", value = FALSE), 
    checkboxInput("standardize_ori_position", "Standardize ORI position", value = FALSE)
  ),
  plotOutput("plasmidMap"),
  DTOutput("featuresTable"),
  tags$script(HTML("
  $(document).on('change', '.row-select', function() {
    var checked = [];
    $('.row-select:checked').each(function() {
      checked.push(parseInt($(this).data('row')));
    });
    Shiny.setInputValue('selected_rows', checked, {priority: 'event'});
  });
"))
)

# Define server logic required to draw a histogram ----
server <- function(input, output) {
  # cat in red bold letters

  cat(paste0("\033[1;31m", "Starting server", "\033[0m\n"))

  genbank <- reactiveVal(NULL)
  features_df <- reactiveVal(NULL)  # now a reactiveVal, not reactive()
  selected_rows <- reactiveVal(integer(0))
  seq_length <- reactiveVal(NULL)
  observeEvent(input$genbank, {
  
    gbk <- read_gb(file = input$genbank$datapath)
    print(gbk$length)
    seq_length(gbk$length)
    # remove features with NA start_end and add index column
    for (i in rev(seq_along(gbk$features))) {
      if (is.na(sum(gbk$features[[i]]$start_end))) {
        gbk$features[[i]] <- NULL
      } else {
        gbk$features[[i]]$index <- i
      }
    }
    df <- as.data.frame(gbk)
    # remove useless featured added by geneious prime
    df <- remove_junk_features(df)
    df <- fix_geneious_type(df)

    df$selected <- sprintf(
  '<input type="checkbox" class="row-select" data-row="%d" checked/>', 
  seq_len(nrow(df))
)
  features_df(df)
  selected_rows(seq_len(nrow(df)))  # select all rows by default
  })
  
  output$plasmidMap <- renderPlot({
    req(features_df())

    plot_plasmid(features_df()[selected_rows(), ], name = input$genbank$name, seq_length = seq_length()) +  # now uses features_df()
      {if(input$circular) ggplot2::coord_polar() else ggplot2::coord_cartesian()} + 
      ggplot2::scale_y_continuous(limits = NULL) +
      scale_fill_manual(values = Geneious_palette) + 
      theme()
  })
  
  output$featuresTable <- renderDT({
    req(features_df())
    features_df()},
    escape = FALSE,
    editable = list(target = "cell", disable = list(columns = which(names(features_df()) == "selected") - 1)),
    options = list(
    columnDefs = list(list(orderable = FALSE, 
    targets = which(names(features_df()) == "selected") - 1, 
    selection = "none"
    )), 
    pageLength = 50
  )
  )
  
  # Cell edits update features_df(), which reactively updates both table and plot
  observeEvent(input$featuresTable_cell_edit, {
    info <- input$featuresTable_cell_edit
    features_df(editData(features_df(), info))
  })

  observeEvent(input$selected_rows, {
    selected_rows(input$selected_rows)
    
  })

  observeEvent(input$standardize_ori_position, {
    if (input$standardize_ori_position) {
      standardize_ori_position(features_df())
    }
   
  })

}

shinyApp(ui = ui, server = server)