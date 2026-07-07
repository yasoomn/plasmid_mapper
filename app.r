library(shiny)
library(bslib)
library(ggplot2)
library(dplyr)
library(reactable)
library(reactable.extras)
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


color_palettes <- list(
Geneious_palette = c(
	"CDS" = "#FFFF00", 
	"gene" = "#00B200", 
	"misc_feature" = "#AAAAAA", 
	"Promoter" = "#BAFF00", 
	"promoter" = "#BAFF00",
	"terminator" = "#FF5300",
	"Terminator" = "#FF5300",
	"rRNA" = "#FF0000", 
	"tRNA" = "#FF00A0", 
	"rep_origin" = "#50B3FF",
	"5'UTR" = "#F5F5F5", 
	"3'UTR" = "#F5F5F5",
	"intron" = "#F5F5F5", 
	"regulatory" = "#0EA089",
	"misc_recomb" = "#FF0000"), 
	Another_palette = c(
	"CDS" = "#4C72B0",           # blue
	"gene" = "#55A868",          # green
	"misc_feature" = "#8C8C8C",  # neutral gray
	"Promoter" = "#8172B2",      # purple
	"promoter" = "#8172B2",
	"terminator" = "#C44E52",    # red
	"Terminator" = "#C44E52",
	"rRNA" = "#DD8452",          # orange
	"tRNA" = "#CC79A7",          # pink/magenta
	"rep_origin" = "#64B5CD",    # light blue/teal
	"5'UTR" = "#EAEAEA",         # light gray
	"3'UTR" = "#EAEAEA",
	"intron" = "#EAEAEA",
	"regulatory" = "#009E73",    # teal-green
	"misc_recomb" = "#B22222"    # firebrick red
)
)



my_theme <- bs_theme(
	version = 5,
	bootswatch = "flatly",
	primary = "#279b3e"
)


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
		textInput("custom_seq_name", "Custom Sequence Name", value = ""),
		textInput("custom_seq_length", "Custom Sequence Length", value = ""),
		# TODO add checkbox for showing/hiding sequence name
		# TODO add input for font size and wrap threshold
		# TODO add input for color palette
		# TODO switch table to reactable for better performance and more features
		# TODO find a way of highlighting important features. Maybe by decreasing the opacity of unselected features. Or by adding a border to selected features.
		radioButtons("filetype_download", "Download Image", choices = c("png", "svg"), selected = "svg"),
		downloadButton("downloadMap", "Download")
	),
	plotOutput("plasmidMap"),
	 navset_card_underline(
		nav_panel("Features Table", 
			h3("Features Table"),
			"Select features from the table to display in the plasmid map. Double click to edit the text.",
			DTOutput("featuresTable")),
			nav_panel("Colors", 
				h3("Color Palette"),
				selectInput("color_palette", "Select Color Palette", choices = names(color_palettes), selected = "Geneious_palette"),
				reactableOutput("colorTable")
				
	)
 

),
reactable.extras::reactable_extras_dependency()
)

server <- function(input, output) {

	cat(paste0("\033[1;31m", "Starting server", "\033[0m\n"))

	genbank <- reactiveVal(NULL)
	features_df <- reactiveVal(NULL)  # now a reactiveVal, not reactive()
	seq_length <- reactiveVal(NULL)
	render_trigger <- reactiveVal(0)
	ori_shift <- reactiveVal(TRUE)
	colorPalettes <- reactiveVal(color_palettes)

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
	
	plasmidMap <- reactive({
			req(features_df())

		plot_plasmid(features_df()[input$featuresTable_rows_selected, ], 
		name = ifelse(input$custom_seq_name == "", input$genbank$name, input$custom_seq_name), 
		seq_length = ifelse(input$custom_seq_length == "", seq_length(), as.numeric(input$custom_seq_length))) +  
			{if(input$circular) ggplot2::coord_polar() else ggplot2::coord_cartesian()} + 
			ggplot2::scale_y_continuous(limits = NULL) +
			ggplot2::scale_fill_manual(values = colorPalettes()[[input$color_palette]]) + 
			theme()
	})

	output$plasmidMap <- renderPlot({
		plasmidMap()
	})
	
	output$featuresTable <- renderDT({
	render_trigger()                # <- the ONLY reactive dependency for a full rebuild
	df <- isolate(features_df())    # read data without depending on it
	req(df)

	datatable(
		df,
		class = "hover",
		selection = list(mode = "multiple",
											selected = which(!str_detect(df$type, "primer_bind") & !str_detect(df$note, "sequence: "))),
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
	backgroundColor = styleEqual(levels = names(colorPalettes()[[input$color_palette]]), values = unname(colorPalettes()[[input$color_palette]]))
)
})

proxy <- dataTableProxy("featuresTable")

observeEvent(input$featuresTable_cell_edit, {
	info <- input$featuresTable_cell_edit
	df <- features_df()
	df <- editData(df, info, rownames = TRUE)
	features_df(df)

	replaceData(proxy, df, resetPaging = FALSE, 
		rownames = TRUE,
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
		paste0(tools::file_path_sans_ext(input$genbank$name), ".", input$filetype_download)
	},
	content = function(file) {
		if (input$filetype_download == "svg") {
			svg(file, width = 12, height = 6)
		} else {
			png(file, width = 1200, height = 600)
		}
		print(plasmidMap())
		dev.off()
	}
)
output$colorTable <- renderReactable({
	data.frame(
		Type = names(colorPalettes()[[input$color_palette]]),
		Color = unname(colorPalettes()[[input$color_palette]])
	) |> reactable(
		columns = list(
			Color = reactable::colDef(cell = reactable.extras::text_extra("customColorPalette", class = "table-text"), style = function(value) list(background = value, color = ifelse(value == "#FFFF00", "black", "white")))
			
		))
})

observeEvent(input$customColorPalette, {
	# update the colorPalettes reactiveVal with the new colors
	req(input$customColorPalette$row, input$customColorPalette$value)
	if (nchar(input$customColorPalette$value) == 7) {
	cp <- colorPalettes()
	cp[[input$color_palette]][input$customColorPalette$row] <- input$customColorPalette$value
	colorPalettes(cp)
	}
	
})



}

shinyApp(ui = ui, server = server)