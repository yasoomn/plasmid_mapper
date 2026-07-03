library(shiny)
library(bslib)
library(ggplot2)
library(dplyr)
library(DT)
library(stringr)

## donwalod from github if it is not installed
if (!requireNamespace("plasmapR", quietly = TRUE)) {
  # Using my custom fork of plasmapR
  remotes::install_github("yasoomn/plasmapR")
}
library(plasmapR)

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

# pure coordinate shift, usable both forward and backward
shift_positions <- function(df, seq_length, shift) {
  wrap_pos <- function(pos) ((pos - shift - 1) %% seq_length) + 1
  df$start <- wrap_pos(df$start)
  df$end   <- wrap_pos(df$end)
  df[order(df$start), ]
}

standardize_ori_position <- function(df, seq_length, ori_row = NULL) {
  ori_candidates <- which(df$type == "rep_origin")
  if (length(ori_candidates) == 0) {
    warning("No rep_origin feature found; returning data unchanged")
    return(list(df = df, shift = 0))
  }
  if (is.null(ori_row)) ori_row <- ori_candidates[1]
  
  shift <- df$start[ori_row] - 1
  list(df = shift_positions(df, seq_length, shift), shift = shift)
}

Geneious_palette = c("CDS" = "#FFFF00", 
  "gene" = "#00B200", 
  "misc_feature" = "#AAAAAA", 
  "Promoter" = "#BAFF00", 
  "promoter" = "#BAFF00",
  "terminator" = "#FF5300", 
  "rRNA" = "#FF0000", 
  "tRNA" = "#FF00A0", 
  "rep_origin" = "#50B3FF",
  "5'UTR" = "#F5F5F5", 
  "3'UTR" = "#F5F5F5",
  "intron" = "#F5F5F5", 
  "regulatory" = "#0EA089")

my_theme <- bs_theme(
  version = 5,
  bootswatch = "flatly",
  primary = "#279b3e"
)

# Define UI for app that draws a histogram ----
ui <- page_sidebar(
  
  title = "Plasmid Viewer",
  theme = my_theme,
  tags$head(
  tags$style(HTML("
    table.dataTable tbody tr.selected td,
    table.dataTable tbody td.selected {
      box-shadow: inset 0 0 0 9999px #7b8a8b40 !important;
      border-left: 3px solid #89ab94 !important;
      color: #000 !important;
    }
    :root {
      --dt-row-selected: transparent !important;
    }
    table.dataTable tbody tr:hover, table.dataTable tbody tr:hover td {
      background-color: rgba(255, 192, 203, 0.15) !important;
    }
  "))
),
  # Sidebar panel for inputs ----
  sidebar = sidebar(
    # Input: Slider for the number of bins ----
    fileInput("genbank", "Choose a Genbank File"),
    #checkboxInput("remove_junk", "Remove Junk Features", value = TRUE), 
    #checkboxInput("invert_seq", "Invert sequence", value = FALSE), 
    checkboxInput("circular", "Circular view", value = FALSE), 
    checkboxInput("standardize_ori_position", "ORI at the start", value = TRUE),
    # TODO add input for font size and wrap threshold
    # TODO add input for color palette
    downloadButton("downloadMap", "Download Map")
  ),
  plotOutput("plasmidMap"),
  DTOutput("featuresTable")

)


server <- function(input, output) {


  cat(paste0("\033[1;31m", "Starting server", "\033[0m\n"))

  genbank <- reactiveVal(NULL)
  features_df <- reactiveVal(NULL)  # now a reactiveVal, not reactive()
  selected_rows <- reactiveVal(integer(0))
  seq_length <- reactiveVal(NULL)
  render_trigger <- reactiveVal(0)
  ori_shift <- reactiveVal(TRUE)
  
  observeEvent(input$genbank, {
  gbk <- read_gb(file = input$genbank$datapath)
  seq_length(gbk$length)
  for (i in rev(seq_along(gbk$features))) {
    if (is.na(sum(gbk$features[[i]]$start_end))) {
      gbk$features[[i]] <- NULL
    } else {
      gbk$features[[i]]$index <- i
    }
  }
  df <- as.data.frame(gbk)
  df <- remove_junk_features(df)
  df <- fix_geneious_type(df)
  result <- standardize_ori_position(df, seq_length())
  ori_shift(result$shift)
  df <- result$df
  features_df(df)

  render_trigger(isolate(render_trigger()) + 1)  # only bump on new file
})
  
  output$plasmidMap <- renderPlot({
    req(features_df())

    plot_plasmid(features_df()[input$featuresTable_rows_selected, ], name = input$genbank$name, seq_length = seq_length()) +  # now uses features_df()
      {if(input$circular) ggplot2::coord_polar() else ggplot2::coord_cartesian()} + 
      ggplot2::scale_y_continuous(limits = NULL) +
      scale_fill_manual(values = Geneious_palette) + 
      theme()
  })
  
  output$featuresTable <- renderDT({
  render_trigger()                # <- the ONLY reactive dependency for a full rebuild
  df <- isolate(features_df())    # read data without depending on it
  req(df)

  datatable(
    df,
    class = "hover",
    selection = list(mode = "multiple",
                      selected = which(!(df$type %in% c("primer_bind")))),
    editable = list(target = "cell"),
    rownames = TRUE,
    options = list(
      columnDefs = list(list(visible = FALSE, targets = c(1))),
      pageLength = 50,
      dom = 'ft',
      server = TRUE
    )
  ) |> formatStyle(
  'type',
  backgroundColor = styleEqual(levels = names(Geneious_palette), values = unname(Geneious_palette))
)
})

proxy <- dataTableProxy("featuresTable")

observeEvent(input$featuresTable_cell_edit, {
  info <- input$featuresTable_cell_edit
  df <- features_df()
  df <- editData(df, info, rownames = TRUE)
  features_df(df)

  replaceData(proxy, df, resetPaging = FALSE, 
    rownames = TRUE, ,
    clearSelection = "none")
})

  observeEvent(input$standardize_ori_position, {
  req(features_df(), seq_length())
  df <- features_df()
  
  if (input$standardize_ori_position) {
    result <- standardize_ori_position(df, seq_length())
    features_df(result$df)
    ori_shift(result$shift)
  } else {
    # undo by applying the exact inverse shift
    df_restored <- shift_positions(df, seq_length(), -ori_shift())
    features_df(df_restored)
    ori_shift(0)
  }
  
  render_trigger(isolate(render_trigger()) + 1)
})

output$downloadMap <- downloadHandler(
  filename = function() {
    paste0(tools::file_path_sans_ext(input$genbank$name), "_plasmid_map.png")
  },
  content = function(file) {
    req(features_df())
    ggsave(file, plot = output$plasmidMap, device = "png", width = 10, height = 6)
  }
)

}

shinyApp(ui = ui, server = server)