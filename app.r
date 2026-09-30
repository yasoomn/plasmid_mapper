library(shiny)
library(bslib)
library(ggplot2)
library(dplyr)
library(reactable)
library(reactable.extras)
library(DT)
library(stringr)
library(shinyjqui)

## donwalod from github if it is not installed
if (!requireNamespace("plasmapR", quietly = TRUE)) {
	# Using my fork of plasmapR
	remotes::install_github("yasoomn/plasmapR")
}
library(plasmapR)

remove_junk_features <- function(df) {
	
	junk_features.name <- c("Geneious type: ligation", "Geneious type: Fragment")
	junk_features.type <- c("source", "Fragment", "ligation")
	
	df <- df |> filter(!note %in% junk_features.name, !type %in% junk_features.type)

	return(df)
}

fix_geneious_type <- function(df) {
	print(df$name)
	if ("note" %in% colnames(df)) {
		df$name[is.na(df$name)] <- " "
		df$type[df$note == "Geneious type: promoter"] <- "Promoter"
		df$type[df$note == "Geneious type: terminator"] <- "terminator"
		df$type[df$note == "Geneious type: origin of replication"] <- "rep_origin"
		df$type[df$note == "Geneious type: 5'_UTR"] <- "5'UTR"
		df$type[df$note == "Geneious type: 3'_UTR"] <- "3'UTR"
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
	"cds" = "#FFFF00", 
	"gene" = "#00B200", 
	"misc_feature" = "#AAAAAA", 
	"promoter" = "#BAFF00",
	"terminator" = "#FF5300",
	"rrna" = "#FF0000", 
	"trna" = "#FF00A0", 
	"rep_origin" = "#50B3FF",
	"5'utr" = "#F5F5F5", 
	"3'utr" = "#F5F5F5",
	"intron" = "#F5F5F5", 
	"regulatory" = "#0EA089",
	"misc_recomb" = "#FF0000", 
	"mrna" = "#BE3232", 
	"polya_signal" = "#00CCCC"), 
	Another_palette = c(
	"cds" = "#4C72B0",           # blue
	"gene" = "#55A868",          # green
	"misc_feature" = "#8C8C8C",  # neutral gray
	"promoter" = "#8172B2",      # purple
	"terminator" = "#C44E52",    # red
	"rrna" = "#DD8452",          # orange
	"trna" = "#CC79A7",          # pink/magenta
	"rep_origin" = "#64B5CD",    # light blue/teal
	"5'utr" = "#EAEAEA",         # light gray
	"3'utr" = "#EAEAEA",
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
	# Sidebar panel for inputs ----
	sidebar = sidebar(
		# Input: Slider for the number of bins ----
		fileInput("genbank", "Choose a Genbank File"),
		#checkboxInput("remove_junk", "Remove Junk Features", value = TRUE), 
		#checkboxInput("invert_seq", "Invert sequence", value = FALSE), 
		checkboxInput("circular", "Circular view", value = FALSE), 
		checkboxInput("standardize_ori_position", "ORI at the start", value = TRUE),
		textInput("custom_seq_name", "Custom Sequence Name", value = ""),
		numericInput("custom_seq_length", "Custom Sequence Length", value = NA, min = 1),
		#numericInput("base_size", "Base Font Size", value = 14, min = 1),
		
		# TODO add suport for truncated features
		# TODO add option to show feature size in the lable, togglable by individual feature
		# TODO add input for font size and wrap threshold
		# TODO find a way of highlighting important features. Maybe by decreasing the opacity of unselected features. Or by adding a border to selected features.
		radioButtons("filetype_download", "Download Image", choices = c("png", "svg"), selected = "svg"),
		downloadButton("downloadMap", "Download")
	),
	div(
		id = "plot_container", 
		plotOutput("plasmidMap")
		),
	 navset_card_underline(
		nav_panel("Features Table", 
			"Select features from the table to display in the plasmid map.",
			reactableOutput("featuresTable")),
			nav_panel("Colors", 
				h3("Color Palette"),
				selectInput("color_palette", "Select Color Palette", choices = names(color_palettes), selected = "Geneious_palette"),
				reactableOutput("colorTable")
				
	),
	height = "1200"
 
),
  tags$head(
    tags$style(HTML("
      #plot_container {
        resize: both;
        overflow: hidden;
        border: 1px solid #ccc;
        width: 100%;
        height: 400px;
        min-width: 200px;
        min-height: 150px;
        padding: 0;
      }
      /* Let the browser stretch the existing image instantly while dragging */
      #plot_container img {
        width: 100% !important;
        height: 100% !important;
        object-fit: contain;
      }
    "))
  ),  
  tags$script(HTML("
    $(function() {
      var el = document.getElementById('plot_container');
      var timeout = null;
 
      var ro = new ResizeObserver(function(entries) {
        for (var entry of entries) {
          var w = entry.contentRect.width;
          var h = entry.contentRect.height;
 
          // Debounce: only notify Shiny after resizing pauses for 250ms
          clearTimeout(timeout);
          timeout = setTimeout(function() {
            Shiny.setInputValue('plot_dims', {w: w, h: h}, {priority: 'event'});
          }, 250);
        }
      });
 
      ro.observe(el);
    });
  ")),

reactable.extras::reactable_extras_dependency()
)

server <- function(input, output) {

	cat(paste0("\033[1;31m", "Starting server", "\033[0m\n"))

### Reactive values
	genbank <- reactiveVal(NULL)
	features_df <- reactiveVal(NULL) 
	seq_length <- reactiveVal(NULL)
	render_trigger <- reactiveVal(0)
	ori_shift <- reactiveVal(TRUE)
	colorPalettes <- reactiveVal(color_palettes)

	features_table_selected <- reactive(getReactableState("featuresTable", "selected"))
	feature_table_name_edit <- debounce(reactive({
		input[[paste0("feature_table_name_", render_trigger())]]
	}), millis = 500)
	feature_table_type_edit <- debounce(reactive({
		input[[paste0("feature_table_type_", render_trigger())]]
	}), millis = 500)

	feature_table_show_label_edit <- debounce(reactive({
		input[[paste0("feature_table_show_label_", render_trigger())]]
	}), millis = 200)
	feature_table_show_size_edit <- debounce(reactive({
		input[[paste0("feature_table_show_size_", render_trigger())]]
	}), millis = 200)
	feature_table_direction_edit <- debounce(reactive({
		input[[paste0("feature_table_direction_", render_trigger())]]
	}), millis = 200)
	dims <- reactive({
    if (is.null(input$plot_dims)) list(w = 500, h = 400) else input$plot_dims
  	})


	### Triggered by file upload
	observeEvent(input$genbank, {
	parsed <- tryCatch({
		gbk <- read_gb(file = input$genbank$datapath)
		if (!is.numeric(gbk$length) || length(gbk$length) != 1 ||
			is.null(gbk$features) || !is.list(gbk$features)) {
			stop("The uploaded file is not a valid GenBank file.")
		}

		## Remove features that have invalid start or end positions
		for (i in rev(seq_along(gbk$features))) {
			if (is.na(sum(gbk$features[[i]]$start_end))) {
				gbk$features[[i]] <- NULL
			} else {
				gbk$features[[i]]$index <- i
			}
		}
		df <- as.data.frame(gbk)
		## Remove/fix features that are not useful for visualization
		df <- remove_junk_features(df)
		df <- fix_geneious_type(df)
		## Create new column to control whether to show the label for each feature. Default is to show all labels.
		df$show_label <- 1
		df$show_size <- 0
		df$size <- df$end - df$start + 1

		result <- standardize_ori_position(df, gbk$length)
		list(gbk = gbk, df = result$df, shift = result$shift)
	}, error = function(error) {
		showNotification(
			"The uploaded file could not be read as a GenBank file.",
			type = "error"
		)
		NULL
	})

	if (is.null(parsed)) return()
## Store the parsed GenBank data and features in reactive values
	genbank(parsed$gbk)
	seq_length(parsed$gbk$length)
	ori_shift(parsed$shift)
	features_df(parsed$df)
	render_trigger(isolate(render_trigger()) + 1)  # only bump on new file

})
	
	## Create the plasmid map plot based on the selected features and user inputs
	plasmidMap <- reactive({
		req(features_df(), seq_length())
		print(head(features_df()))
		df <- features_df() 
		df$type <- tolower(df$type)
		df$name <- ifelse(df$show_size == 1, paste0(df$name, " (", df$size, " bp)"), df$name)
		plot_plasmid(df[features_table_selected(), ], 
		name =  ifelse(input$custom_seq_name == "", input$genbank$name, input$custom_seq_name),
		seq_length = ifelse(is.na(input$custom_seq_length), seq_length(), as.numeric(input$custom_seq_length))) +  
			{if(input$circular) ggplot2::coord_polar() else ggplot2::coord_cartesian()} + 
			ggplot2::scale_y_continuous() +
			ggplot2::scale_fill_manual(values = colorPalettes()[[input$color_palette]])
			
	})

	output$plasmidMap <- renderPlot({
		plasmidMap()
	},
	width  = function() dims()$w,
	height = function() dims()$h
)
	
	## Render the features table
	output$featuresTable <- renderReactable({
	render_trigger()
	df <- isolate(features_df())    # read data without depending on it
	req(df)

	name_input_id <- paste0("feature_table_name_", render_trigger())
	type_input_id <- paste0("feature_table_type_", render_trigger())
	show_label_input_id <- paste0("feature_table_show_label_", render_trigger())
	show_size_input_id <- paste0("feature_table_show_size_", render_trigger())
	direction_input_id <- paste0("feature_table_direction_", render_trigger())

	reactable(
		df, 
		defaultPageSize = 100,
		striped = TRUE,
		compact = TRUE,
		searchable = TRUE,
		selection = "multiple",
		highlight = TRUE,
		defaultSelected = which(!str_detect(df$type, "primer_bind") & !str_detect(df$note, "sequence: ") & !str_detect(df$type, "PCR_primer") & !str_detect(df$type, "gene") & !str_detect(df$type, "primer")),
		columns = list(
			index = colDef(show = FALSE), 
			name = colDef(cell = reactable.extras::text_extra(name_input_id, class = "table-text")),
			type = colDef(cell = reactable.extras::text_extra(type_input_id, class = "table-text")), 
			direction = colDef(cell = reactable.extras::dropdown_extra(direction_input_id, choices = c(1, 0, -1), class = "table-dropdown")),
			show_label = colDef(cell = reactable.extras::checkbox_extra(show_label_input_id, class = "checkbox-extra")),
			size = colDef(cell = function(value) {
				if (is.na(value)) return("")
				format(value, big.mark = ",")
			}),
			show_size = colDef(cell = reactable.extras::checkbox_extra(show_size_input_id, class = "checkbox-extra")), 
			start = colDef(show = FALSE), 
			end = colDef(show = FALSE)
			)
	
		)
	})

## Update the features_df reactiveVal when the user edits the name or type in the features table
observeEvent(feature_table_name_edit(), {
	edit <- feature_table_name_edit()
	req(edit$row, edit$value)

	# fix the row index to match the original data frame
	info <- data.frame(
		row = edit$row,
		col = 2,
		value = edit$value
	)
	df <- features_df()
	df <- editData(df, info, rownames = TRUE)
	# update the features_df reactiveVal with the new name
	features_df(df)


}, ignoreNULL = TRUE)

## Update the features_df reactiveVal when the user edits the type in the features table
observeEvent(feature_table_type_edit(), {
	edit <- feature_table_type_edit()
	req(edit$row, edit$value)

	# fix the row index to match the original data frame
	info <- data.frame(
		row = edit$row,
		col = 3,
		value = edit$value
	)
	df <- features_df()
	df <- editData(df, info, rownames = TRUE)
	# update the features_df reactiveVal with the new name
	features_df(df)


}, ignoreNULL = TRUE)

observeEvent(feature_table_show_label_edit(), {
  edit <- feature_table_show_label_edit()
	# `FALSE` is a valid checkbox value; `req(edit$value)` would drop the
	# event when the checkbox is unchecked.
	req(!is.null(edit$row), !is.null(edit$value))
  # fix the row index to match the original data frame
  info <- data.frame(
    row = edit$row,
    col = 8,
    value = edit$value
  )
  df <- features_df()
  df <- editData(df, info, rownames = TRUE)
  # update the features_df reactiveVal with the new name
  features_df(df)
}, ignoreNULL = TRUE)

observeEvent(feature_table_show_size_edit(), {
  edit <- feature_table_show_size_edit()
  req(!is.null(edit$row), !is.null(edit$value))
  # fix the row index to match the original data frame
  info <- data.frame(
	row = edit$row,
	col = 9,
	value = edit$value
  )
  df <- features_df()
  df <- editData(df, info, rownames = TRUE)
  # update the features_df reactiveVal with the new name
  features_df(df)
}, ignoreNULL = TRUE)

observeEvent(feature_table_direction_edit(), {
  edit <- feature_table_direction_edit()
  req(!is.null(edit$row), !is.null(edit$value))
  # fix the row index to match the original data frame
  info <- data.frame(
	row = edit$row,
	col = 7,
	value = edit$value
  )
  df <- features_df()
  df <- editData(df, info, rownames = TRUE)
  # update the features_df reactiveVal with the new name
  features_df(df)
}, ignoreNULL = TRUE)

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

})

output$downloadMap <- downloadHandler(
	filename = function() {
		paste0(tools::file_path_sans_ext(input$genbank$name), ".", input$filetype_download)
	},
	content = function(file) {
		if (input$filetype_download == "svg") {
			svg(file, width = dims()$w/dims()$h * 6, height = 6)
		} else {
			png(file, width = dims()$w, height = dims()$h)
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
		pagination = FALSE,
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

observeEvent(features_table_selected(), {
	# when the user selects or deselects rows, refresh the plasmid ma
	# print(features_table_selected())

})

}
shinyApp(ui = ui, server = server)
