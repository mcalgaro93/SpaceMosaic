#' Run Interactive Spatial Plotter
#'
#' Launches a local Shiny application for interactive visualization of SpaceMosaic results.
#' This app allows users to overlay cells and patches, and fill patches using
#' either meta-analysis Z-scores or arbitrary patch-level annotations.
#'
#' @param spe A \code{SpatialExperiment} object containing cell coordinates and metadata.
#' @param patches Vector of patch assignments for each cell (output of \code{getPatches}).
#' @param metats Optional matrix of meta-analysis Z-scores (genes x patches).
#' @param meaningful_vars Optional vector of column names from \code{colData(spe)} to be used for coloring and highlighting. If \code{NULL}, all columns are used.
#' @param patch_data Optional data frame with one row per patch and arbitrary
#'   numeric or categorical annotations. It must contain a `patch` column, or
#'   use patch identifiers as row names.
#'
#' @export
runInteractivePlotter <- function(spe, patches, metats = NULL,
                                  meaningful_vars = NULL,
                                  patch_data = NULL) {
  if (!requireNamespace("shiny", quietly = TRUE)) {
    stop("Package 'shiny' is required for this function. Please install it.")
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for this function. Please install it.")
  }
  if (!requireNamespace("SpatialExperiment", quietly = TRUE)) {
    stop("Package 'SpatialExperiment' is required for this function. Please install it.")
  }

  # Extract coordinates and metadata
  xy <- SpatialExperiment::spatialCoords(spe)
  if (is.null(rownames(xy))) {
    rownames(xy) <- paste0("cell_", seq_len(nrow(xy)))
  }

  meta <- as.data.frame(SummarizedExperiment::colData(spe))

  # Determine available variables for coloring and highlighting
  avail_vars <- if (!is.null(meaningful_vars)) {
    intersect(meaningful_vars, colnames(meta))
  } else {
    colnames(meta)
  }

  if (length(avail_vars) == 0) {
    stop("No meaningful variables found in the object's colData.")
  }

  # Helper to align named vectors to xy rownames
  align_to_xy <- function(vec, xy_mat, name_label) {
    if (is.null(vec)) return(rep(NA, nrow(xy_mat)))
    if (length(vec) == nrow(xy_mat)) return(vec)
    if (!is.null(names(vec)) && !is.null(rownames(xy_mat))) {
      res <- rep(NA, nrow(xy_mat))
      names(res) <- rownames(xy_mat)
      common <- intersect(names(vec), rownames(xy_mat))
      res[common] <- vec[common]
      return(res)
    }
    stop(sprintf("Input %s length (%d) does not match xy rows (%d) and cannot be aligned via names.", 
                 name_label, length(vec), nrow(xy_mat)))
  }

  # Align patches to xy
  patches_aligned <- align_to_xy(patches, xy, "patches")

  # Prepare polygon data and attach arbitrary patch-level annotations.
  poly_df <- getPatchPolys(xy, patches_aligned, patch_data = patch_data)
  patch_vars <- setdiff(names(poly_df), c("x", "y", "patch"))

  if (!is.null(metats)) {
    metats <- as.matrix(metats)
    if (is.null(rownames(metats)) || is.null(colnames(metats))) {
      stop("metats must have gene row names and patch column names.")
    }
  }
  if (is.null(metats) && length(patch_vars) == 0L) {
    stop("Provide metats, patch_data with at least one annotation, or both.")
  }

  patch_fill_choices <- c(
    if (!is.null(metats)) c("Meta-analysis Z-score" = ".zscore"),
    if (length(patch_vars) > 0L) {
      stats::setNames(patch_vars, paste0("Patch annotation: ", patch_vars))
    }
  )

  # Ensure cell data is in a tidy format
  cell_df <- data.frame(
    x = xy[, 1],
    y = xy[, 2],
    patch = patches_aligned,
    stringsAsFactors = FALSE
  )
  # Merge with colData
  cell_df <- cbind(cell_df, meta)

  # Available layers for ordering
  available_layers <- c("Background Cells", "Patches", "Highlighted Cells")

  continuous_palette_choices <- c(
    "Viridis" = "Viridis",
    "Cividis" = "Cividis",
    "Plasma" = "Plasma",
    "Inferno" = "Inferno",
    "Yellow-Orange-Red" = "YlOrRd",
    "Blue-Red" = "Blue-Red 3"
  )
  qualitative_palette_choices <- c(
    "Automatic" = "Automatic",
    "Okabe-Ito" = "Okabe-Ito",
    "Dark 3" = "Dark 3",
    "Set 2" = "Set 2",
    "Set 3" = "Set 3",
    "Harmonic" = "Harmonic",
    "Dynamic HCL" = "Dynamic"
  )

  qualitative_colors <- function(n, palette = "Automatic", dark = FALSE) {
    if (n == 0) return(character(0))

    if (palette == "Automatic") {
      if (n <= 8) {
        # Color-vision-deficiency-friendly Okabe-Ito palette.
        colors <- c(
          "#E69F00", "#56B4E9", "#009E73", "#F0E442",
          "#0072B2", "#D55E00", "#CC79A7", "#000000"
        )[seq_len(n)]
        if (dark) colors[colors == "#000000"] <- "#F2F2F2"
        return(colors)
      }

      # Set 3 remains easy to scan for medium-sized sets. For larger sets,
      # Dynamic distributes hues around the full HCL color wheel and always
      # returns exactly the requested number of colors.
      palette <- if (n <= 12) "Set 3" else "Dynamic"
    }

    if (palette == "Okabe-Ito") {
      base <- c(
        "#E69F00", "#56B4E9", "#009E73", "#F0E442",
        "#0072B2", "#D55E00", "#CC79A7", "#000000"
      )
      if (n <= length(base)) {
        colors <- base[seq_len(n)]
        if (dark) colors[colors == "#000000"] <- "#F2F2F2"
        return(colors)
      }

      # Okabe-Ito contains only eight colors; extend it without recycling
      # colors when the selected variable unexpectedly has more categories.
      colors <- c(base, grDevices::hcl.colors(n - length(base), "Dynamic"))
      if (dark) colors[colors == "#000000"] <- "#F2F2F2"
      return(colors)
    }

    grDevices::hcl.colors(n, palette = palette)
  }

  ui <- shiny::fluidPage(
    shiny::tags$head(
      shiny::tags$script(shiny::HTML("
        Shiny.addCustomMessageHandler('spaceMosaicDarkMode', function(enabled) {
          document.body.classList.toggle('space-mosaic-dark', enabled);
        });
      ")),
      shiny::tags$style(shiny::HTML("
          body,
          .well,
          .form-control,
          .selectize-input,
          .selectize-dropdown,
          #legend_panel {
            transition: background-color 180ms ease, color 180ms ease,
                        border-color 180ms ease;
          }

          body.space-mosaic-dark {
            background-color: #111417;
            color: #e8ecef;
          }

          body.space-mosaic-dark .well {
            background-color: #20252a;
            border-color: #3d454c;
          }

          body.space-mosaic-dark .help-block {
            color: #aeb8c0;
          }

          body.space-mosaic-dark .form-control,
          body.space-mosaic-dark .selectize-input,
          body.space-mosaic-dark .selectize-dropdown {
            color: #e8ecef;
            background-color: #171b1f;
            border-color: #4b555e;
          }

          body.space-mosaic-dark .selectize-dropdown .active {
            color: #ffffff;
            background-color: #34424d;
          }

          body.space-mosaic-dark #legend_panel {
            background-color: #161a1d;
            border-color: #3d454c;
          }

          #spatial_plot_container,
          #spatialPlot,
          #legend_ui,
          #legendPlot {
            background-color: #ffffff;
            border: 0 !important;
            box-shadow: none !important;
          }

          #spatialPlot img,
          #legendPlot img {
            display: block;
            background-color: transparent !important;
            border: 0 !important;
            box-shadow: none !important;
          }

          body.space-mosaic-dark #spatial_plot_container,
          body.space-mosaic-dark #spatialPlot,
          body.space-mosaic-dark #legend_ui,
          body.space-mosaic-dark #legendPlot {
            background-color: #161a1d;
          }

          .export-buttons {
            display: flex;
            flex-wrap: wrap;
            gap: 8px;
            margin-bottom: 8px;
          }

          .dark-mode-control .checkbox label {
            display: flex;
            align-items: center;
          }

          .dark-mode-control input[type='checkbox'] {
            appearance: none;
            -webkit-appearance: none;
            width: 42px;
            height: 22px;
            margin: 0 9px 0 0;
            border: 0;
            border-radius: 11px;
            background-color: #aab2b8;
            cursor: pointer;
            position: relative;
            transition: background-color 160ms ease;
          }

          .dark-mode-control input[type='checkbox']::before {
            content: '';
            position: absolute;
            width: 18px;
            height: 18px;
            left: 2px;
            top: 2px;
            border-radius: 50%;
            background: #ffffff;
            box-shadow: 0 1px 3px rgba(0, 0, 0, 0.35);
            transition: transform 160ms ease;
          }

          .dark-mode-control input[type='checkbox']:checked {
            background-color: #337ab7;
          }

          .dark-mode-control input[type='checkbox']:checked::before {
            transform: translateX(20px);
          }

          #main_panel_container {
            display: flex;
            align-items: stretch;
            gap: 16px;
            width: 100%;
            min-width: 0;
          }

          #spatial_plot_container {
            flex: 1 1 auto;
            min-width: 0;
          }

          #spatialPlot {
            width: 100% !important;
            height: calc(100vh - 150px) !important;
            min-height: 400px;
          }

          #legend_panel {
            flex: 0 0 280px;
            max-height: calc(100vh - 150px);
            overflow-y: auto;
            overflow-x: auto;
            padding: 8px 10px;
            border-left: 1px solid #e5e5e5;
          }

          @media (max-width: 991px) {
            #main_panel_container {
              flex-direction: column;
            }

            #spatialPlot {
              height: 68vh !important;
            }

            #legend_panel {
              flex-basis: auto;
              width: 100%;
              max-height: none;
              overflow-y: visible;
              border-left: 0;
              border-top: 1px solid #e5e5e5;
            }
          }
        "))
    ),
    shiny::titlePanel("SpaceMosaic Interactive Plotter"),
    shiny::sidebarLayout(
      shiny::sidebarPanel(
        shiny::wellPanel(
          shiny::h4("Layer Management"),
          shiny::p("Define the order of visual layers from bottom to top to control occlusion."),
          shiny::selectInput("layer_order", "Layer Order (Bottom to Top)",
                      choices = available_layers, 
                      selected = available_layers, 
                      multiple = TRUE)
        ),

        shiny::wellPanel(
          shiny::h4("Appearance"),
          shiny::div(
            class = "dark-mode-control",
            shiny::checkboxInput("dark_mode", "Dark Mode", value = FALSE)
          ),
          shiny::p("Use a dark background for the interface, spatial map, and legend."),
          shiny::tags$hr(),
          shiny::div(
            class = "export-buttons",
            shiny::downloadButton(
              "download_plot_code", "Export R code",
              class = "btn-primary"
            ),
            shiny::downloadButton(
              "download_plot_png", "Export PNG",
              class = "btn-success"
            ),
            shiny::downloadButton(
              "download_plot_pdf", "Export PDF",
              class = "btn-danger"
            )
          ),
          shiny::helpText(
            paste0(
              "Download the current view as reproducible ggplot2 code, ",
              "a high-resolution PNG, or a vector PDF."
            )
          )
        ),

        shiny::wellPanel(
          shiny::h4("Cell Coloring"),
          shiny::p("Choose a metadata variable to color all cells and adjust their opacity."),
          shiny::selectInput("color_var", "Color Cells By",
                      choices = avail_vars, 
                      selected = avail_vars[1]),
          shiny::uiOutput("color_palette_ui"),
          shiny::sliderInput("cell_alpha", "Cell Opacity", 0, 1, 0.5, 0.1)
        ),

        shiny::wellPanel(
          shiny::h4("Cell Highlighting"),
          shiny::p("Isolate specific cell categories and customize their highlight color."),
          shiny::selectInput("highlight_var", "Highlight Variable",
                      choices = avail_vars, 
                      selected = avail_vars[1]),
          shiny::uiOutput("highlight_categories_ui"),
          shiny::selectInput("highlight_color", "Highlight Color",
                      choices = c("Black" = "black", "Red" = "red", "Blue" = "blue", 
                                  "Green" = "green", "Yellow" = "yellow", "Magenta" = "magenta"), 
                      selected = "black"),
          shiny::checkboxInput("grey_background", "Grey out other cells", value = TRUE)
        ),

        shiny::wellPanel(
          shiny::h4("Patch Overlay"),
          shiny::p("Color patches by a gene Z-score or any supplied patch annotation."),
          shiny::checkboxInput("show_patches", "Show Patches", value = TRUE),
          shiny::selectInput("patch_fill_var", "Color Patches By",
                      choices = patch_fill_choices,
                      selected = unname(patch_fill_choices[1])),
          shiny::uiOutput("patch_fill_controls_ui"),
          shiny::sliderInput("patch_alpha", "Patch Opacity", 0, 1, 0.5, 0.1)
        )
      ),
      shiny::mainPanel(
        shiny::div(id = "main_panel_container",
            shiny::div(id = "spatial_plot_container",
                shiny::plotOutput("spatialPlot")
            ),
            shiny::div(id = "legend_panel",
                shiny::uiOutput("legend_ui")
            )
        )
      )
    )
  )

  server <- function(input, output, session) {
    shiny::observe({
      session$sendCustomMessage(
        "spaceMosaicDarkMode",
        isTRUE(input$dark_mode)
      )
    })

    output$color_palette_ui <- shiny::renderUI({
      values <- meta[[input$color_var]]

      if (is.numeric(values)) {
        return(shiny::tagList(
          shiny::selectInput(
            "color_palette", "Continuous Palette",
            choices = continuous_palette_choices,
            selected = "Viridis"
          ),
          shiny::helpText("The palette is applied continuously across the observed range.")
        ))
      }

      n_categories <- length(unique(values[!is.na(values)]))
      shiny::tagList(
        shiny::selectInput(
          "color_palette", "Qualitative Palette",
          choices = qualitative_palette_choices,
          selected = "Automatic"
        ),
        shiny::helpText(sprintf(
          paste0(
            "%d observed categories. Automatic uses Okabe-Ito up to 8, ",
            "Set 3 up to 12, and a dynamically generated HCL palette beyond 12."
          ),
          n_categories
        ))
      )
    })

    # Dynamic UI for highlight categories based on the selected variable
    output$highlight_categories_ui <- shiny::renderUI({
      var <- input$highlight_var
      cats <- sort(unique(meta[[var]]))
      shiny::selectInput("highlight_cats", "Categories to Highlight",
                  choices = cats, 
                  multiple = TRUE)
    })

    output$patch_fill_controls_ui <- shiny::renderUI({
      fill_var <- input$patch_fill_var
      if (is.null(fill_var)) return(NULL)

      if (identical(fill_var, ".zscore")) {
        return(shiny::selectInput(
          "gene_select", "Gene for Patch Fill",
          choices = rownames(metats),
          selected = rownames(metats)[1]
        ))
      }

      values <- poly_df[[fill_var]]
      if (is.numeric(values)) {
        return(shiny::selectInput(
          "patch_palette", "Patch Palette",
          choices = continuous_palette_choices,
          selected = "Viridis"
        ))
      }

      shiny::selectInput(
        "patch_palette", "Patch Palette",
        choices = qualitative_palette_choices,
        selected = "Automatic"
      )
    })

    # Reactive expression to build the base ggplot object
    base_plot <- shiny::reactive({
      is_dark <- isTRUE(input$dark_mode)
      background <- if (is_dark) "#161A1D" else "white"
      foreground <- if (is_dark) "#E8ECEF" else "#222222"
      grid_major <- if (is_dark) "#394047" else "#D9D9D9"
      grid_minor <- if (is_dark) "#292F34" else "#EEEEEE"

      p <- ggplot2::ggplot() +
        ggplot2::theme_minimal() +
        ggplot2::coord_fixed() +
        ggplot2::theme(
          plot.background = ggplot2::element_rect(fill = background, color = NA),
          panel.background = ggplot2::element_rect(fill = background, color = NA),
          legend.background = ggplot2::element_rect(fill = background, color = NA),
          legend.box.background = ggplot2::element_rect(fill = background, color = NA),
          legend.key = ggplot2::element_rect(fill = background, color = NA),
          text = ggplot2::element_text(color = foreground),
          axis.text = ggplot2::element_text(color = foreground),
          axis.title = ggplot2::element_text(color = foreground),
          plot.title = ggplot2::element_text(color = foreground),
          plot.subtitle = ggplot2::element_text(color = foreground),
          legend.title = ggplot2::element_text(color = foreground),
          legend.text = ggplot2::element_text(color = foreground),
          panel.grid.major = ggplot2::element_line(color = grid_major),
          panel.grid.minor = ggplot2::element_line(color = grid_minor)
        )

      order <- input$layer_order
      if (is.null(order)) order <- available_layers

      for (layer in order) {
        if (layer == "Background Cells") {
          color_var <- input$color_var
          if (length(input$highlight_cats) > 0 && input$grey_background) {
            cell_df$bg_color <- ifelse(cell_df[[input$highlight_var]] %in% input$highlight_cats, 
                                       "highlight", "grey80")
            p <- p + ggplot2::geom_point(data = cell_df, ggplot2::aes(x = x, y = y, color = bg_color),
                                size = 0.1, alpha = input$cell_alpha) +
                  ggplot2::scale_color_manual(values = c("highlight" = "transparent", "grey80" = "grey80"), guide = "none")
          } else {
            p <- p + ggplot2::geom_point(data = cell_df, ggplot2::aes(x = x, y = y, color = .data[[color_var]]),
                                size = 0.1, alpha = input$cell_alpha)

            values <- meta[[color_var]]
            selected_palette <- input$color_palette
            if (is.numeric(values)) {
              if (is.null(selected_palette) ||
                  !selected_palette %in% continuous_palette_choices) {
                selected_palette <- "Viridis"
              }
              p <- p + ggplot2::scale_color_gradientn(
                colors = grDevices::hcl.colors(256, selected_palette),
                name = color_var,
                na.value = "grey80"
              )
            } else {
              if (is.factor(values)) {
                categories <- levels(droplevels(values))
              } else {
                categories <- sort(unique(as.character(values[!is.na(values)])))
              }
              if (is.null(selected_palette) ||
                  !selected_palette %in% qualitative_palette_choices) {
                selected_palette <- "Automatic"
              }
              colors <- qualitative_colors(
                length(categories), selected_palette, dark = is_dark
              )
              names(colors) <- categories
              p <- p + ggplot2::scale_color_manual(
                values = colors,
                name = color_var,
                na.value = "grey80",
                drop = TRUE,
                guide = ggplot2::guide_legend(
                  override.aes = list(size = 3, alpha = 1)
                )
              )
            }
          }
        } else if (layer == "Patches") {
          if (input$show_patches) {
            plot_poly_df <- poly_df
            fill_var <- input$patch_fill_var
            if (is.null(fill_var)) fill_var <- unname(patch_fill_choices[1])

            if (identical(fill_var, ".zscore")) {
              gene <- input$gene_select
              if (is.null(gene)) gene <- rownames(metats)[1]
              patch_columns <- match(plot_poly_df$patch, colnames(metats))
              plot_poly_df$.patch_value <- metats[gene, patch_columns]
              fill_scale <- ggplot2::scale_fill_gradient2(
                low = "blue", mid = "white", high = "red", midpoint = 0,
                name = paste0(gene, " Z-score"), na.value = "grey80"
              )
            } else {
              plot_poly_df$.patch_value <- plot_poly_df[[fill_var]]
              selected_palette <- input$patch_palette

              if (is.numeric(plot_poly_df$.patch_value)) {
                if (is.null(selected_palette) ||
                    !selected_palette %in% continuous_palette_choices) {
                  selected_palette <- "Viridis"
                }
                fill_scale <- ggplot2::scale_fill_gradientn(
                  colors = grDevices::hcl.colors(256, selected_palette),
                  name = fill_var,
                  na.value = "grey80"
                )
              } else {
                plot_poly_df$.patch_value <- as.character(plot_poly_df$.patch_value)
                categories <- sort(unique(
                  plot_poly_df$.patch_value[!is.na(plot_poly_df$.patch_value)]
                ))
                if (is.null(selected_palette) ||
                    !selected_palette %in% qualitative_palette_choices) {
                  selected_palette <- "Automatic"
                }
                colors <- qualitative_colors(
                  length(categories), selected_palette, dark = is_dark
                )
                names(colors) <- categories
                fill_scale <- ggplot2::scale_fill_manual(
                  values = colors,
                  name = fill_var,
                  na.value = "grey80",
                  drop = TRUE
                )
              }
            }

            p <- p +
              ggplot2::geom_polygon(
                data = plot_poly_df,
                ggplot2::aes(x = x, y = y, group = patch, fill = .patch_value),
                color = "white", linewidth = 0.2,
                alpha = input$patch_alpha
              ) +
              fill_scale
          }
        } else if (layer == "Highlighted Cells") {
          if (!is.null(input$highlight_cats) && length(input$highlight_cats) > 0) {
            hi_df <- cell_df[cell_df[[input$highlight_var]] %in% input$highlight_cats, ]
            highlight_color <- input$highlight_color
            if (is_dark && identical(highlight_color, "black")) {
              highlight_color <- "white"
            }
            p <- p + ggplot2::geom_point(data = hi_df, ggplot2::aes(x = x, y = y),
                                color = highlight_color, size = 0.2, alpha = 1)
          }
        }
      }

      active_fill <- input$patch_fill_var
      if (is.null(active_fill)) active_fill <- unname(patch_fill_choices[1])
      patch_focus <- if (identical(active_fill, ".zscore")) {
        gene <- input$gene_select
        if (is.null(gene)) gene <- rownames(metats)[1]
        paste("Z-score:", gene)
      } else {
        paste("Patch annotation:", active_fill)
      }
      p + ggplot2::labs(title = "SpaceMosaic Interactive Spatial Map",
               subtitle = patch_focus,
               x = "X", y = "Y")
    })

    # Generate a readable static ggplot2 script from the current controls.
    # Data are deliberately not embedded: the downloaded script expects the
    # same canonical input objects used by this app (spe, patches, and, when
    # applicable, metats and patch_data) to exist in the R session.
    export_plot_code <- shiny::reactive({
      literal <- function(x) {
        paste(utils::capture.output(dput(x)), collapse = "\n")
      }

      is_dark <- isTRUE(input$dark_mode)
      background <- if (is_dark) "#161A1D" else "white"
      foreground <- if (is_dark) "#E8ECEF" else "#222222"
      grid_major <- if (is_dark) "#394047" else "#D9D9D9"
      grid_minor <- if (is_dark) "#292F34" else "#EEEEEE"
      order <- input$layer_order
      if (is.null(order)) order <- available_layers

      active_fill <- input$patch_fill_var
      if (is.null(active_fill)) active_fill <- unname(patch_fill_choices[1])
      gene <- input$gene_select
      if (identical(active_fill, ".zscore") && is.null(gene)) {
        gene <- rownames(metats)[1]
      }

      code <- c(
        "# SpaceMosaic plot exported from runInteractivePlotter()",
        "# Required objects: spe and patches; metats and/or patch_data when used below.",
        "",
        "library(ggplot2)",
        "library(SpatialExperiment)",
        "",
        "if (!exists(\"metats\")) metats <- NULL",
        "if (!exists(\"patch_data\")) patch_data <- NULL",
        "",
        "xy <- SpatialExperiment::spatialCoords(spe)",
        "if (is.null(rownames(xy))) rownames(xy) <- paste0(\"cell_\", seq_len(nrow(xy)))",
        "meta <- as.data.frame(SummarizedExperiment::colData(spe))",
        "",
        "align_to_xy <- function(vec, xy_mat) {",
        "  if (length(vec) == nrow(xy_mat)) return(vec)",
        "  if (!is.null(names(vec)) && !is.null(rownames(xy_mat))) {",
        "    out <- rep(NA, nrow(xy_mat))",
        "    names(out) <- rownames(xy_mat)",
        "    common <- intersect(names(vec), rownames(xy_mat))",
        "    out[common] <- vec[common]",
        "    return(out)",
        "  }",
        "  stop(\"patches cannot be aligned to the spatial coordinates.\")",
        "}",
        "",
        "patches_aligned <- align_to_xy(patches, xy)",
        "poly_df <- SpaceMosaic::getPatchPolys(xy, patches_aligned, patch_data = patch_data)",
        "cell_df <- cbind(",
        "  data.frame(x = xy[, 1], y = xy[, 2], patch = patches_aligned),",
        "  meta",
        ")",
        "",
        paste0("background <- ", literal(background)),
        paste0("foreground <- ", literal(foreground)),
        paste0("grid_major <- ", literal(grid_major)),
        paste0("grid_minor <- ", literal(grid_minor)),
        "",
        "p <- ggplot2::ggplot() +",
        "  ggplot2::theme_minimal() +",
        "  ggplot2::coord_fixed() +",
        "  ggplot2::theme(",
        "    plot.background = ggplot2::element_rect(fill = background, color = NA),",
        "    panel.background = ggplot2::element_rect(fill = background, color = NA),",
        "    legend.background = ggplot2::element_rect(fill = background, color = NA),",
        "    legend.box.background = ggplot2::element_rect(fill = background, color = NA),",
        "    legend.key = ggplot2::element_rect(fill = background, color = NA),",
        "    text = ggplot2::element_text(color = foreground),",
        "    axis.text = ggplot2::element_text(color = foreground),",
        "    axis.title = ggplot2::element_text(color = foreground),",
        "    plot.title = ggplot2::element_text(color = foreground),",
        "    plot.subtitle = ggplot2::element_text(color = foreground),",
        "    legend.title = ggplot2::element_text(color = foreground),",
        "    legend.text = ggplot2::element_text(color = foreground),",
        "    panel.grid.major = ggplot2::element_line(color = grid_major),",
        "    panel.grid.minor = ggplot2::element_line(color = grid_minor)",
        "  )"
      )

      for (layer in order) {
        if (layer == "Background Cells") {
          color_var <- input$color_var
          has_highlights <- !is.null(input$highlight_cats) &&
                            length(input$highlight_cats) > 0

          if (has_highlights && input$grey_background) {
            code <- c(
              code, "",
              paste0("highlight_var <- ", literal(input$highlight_var)),
              paste0("highlight_cats <- ", literal(input$highlight_cats)),
              "cell_df$.background_group <- ifelse(",
              "  cell_df[[highlight_var]] %in% highlight_cats,",
              "  \"highlight\", \"grey80\"",
              ")",
              paste0(
                "p <- p + ggplot2::geom_point(data = cell_df, ggplot2::aes(x, y, color = .background_group), ",
                "size = 0.1, alpha = ", literal(input$cell_alpha), ") +"
              ),
              "  ggplot2::scale_color_manual(",
              "    values = c(highlight = \"transparent\", grey80 = \"grey80\"),",
              "    guide = \"none\"",
              "  )"
            )
          } else {
            values <- meta[[color_var]]
            selected_palette <- input$color_palette
            code <- c(
              code, "",
              paste0("color_var <- ", literal(color_var)),
              paste0(
                "p <- p + ggplot2::geom_point(data = cell_df, ggplot2::aes(x, y, color = .data[[color_var]]), ",
                "size = 0.1, alpha = ", literal(input$cell_alpha), ")"
              )
            )

            if (is.numeric(values)) {
              if (is.null(selected_palette) ||
                  !selected_palette %in% continuous_palette_choices) {
                selected_palette <- "Viridis"
              }
              code <- c(
                code,
                paste0("cell_palette <- ", literal(selected_palette)),
                "p <- p + ggplot2::scale_color_gradientn(",
                "  colors = grDevices::hcl.colors(256, cell_palette),",
                "  name = color_var, na.value = \"grey80\"",
                ")"
              )
            } else {
              categories <- if (is.factor(values)) {
                levels(droplevels(values))
              } else {
                sort(unique(as.character(values[!is.na(values)])))
              }
              if (is.null(selected_palette) ||
                  !selected_palette %in% qualitative_palette_choices) {
                selected_palette <- "Automatic"
              }
              colors <- qualitative_colors(
                length(categories), selected_palette, dark = is_dark
              )
              names(colors) <- categories
              code <- c(
                code,
                paste0("cell_colors <- ", literal(colors)),
                "p <- p + ggplot2::scale_color_manual(",
                "  values = cell_colors, name = color_var, na.value = \"grey80\",",
                "  drop = TRUE,",
                "  guide = ggplot2::guide_legend(override.aes = list(size = 3, alpha = 1))",
                ")"
              )
            }
          }
        } else if (layer == "Patches" && isTRUE(input$show_patches)) {
          code <- c(code, "", paste0("fill_var <- ", literal(active_fill)))

          if (identical(active_fill, ".zscore")) {
            code <- c(
              code,
              paste0("gene <- ", literal(gene)),
              "patch_columns <- match(poly_df$patch, colnames(metats))",
              "poly_df$.patch_value <- metats[gene, patch_columns]",
              "fill_scale <- ggplot2::scale_fill_gradient2(",
              "  low = \"blue\", mid = \"white\", high = \"red\", midpoint = 0,",
              "  name = paste0(gene, \" Z-score\"), na.value = \"grey80\"",
              ")"
            )
          } else {
            patch_values <- poly_df[[active_fill]]
            selected_palette <- input$patch_palette
            code <- c(code, "poly_df$.patch_value <- poly_df[[fill_var]]")

            if (is.numeric(patch_values)) {
              if (is.null(selected_palette) ||
                  !selected_palette %in% continuous_palette_choices) {
                selected_palette <- "Viridis"
              }
              code <- c(
                code,
                paste0("patch_palette <- ", literal(selected_palette)),
                "fill_scale <- ggplot2::scale_fill_gradientn(",
                "  colors = grDevices::hcl.colors(256, patch_palette),",
                "  name = fill_var, na.value = \"grey80\"",
                ")"
              )
            } else {
              categories <- sort(unique(as.character(
                patch_values[!is.na(patch_values)]
              )))
              if (is.null(selected_palette) ||
                  !selected_palette %in% qualitative_palette_choices) {
                selected_palette <- "Automatic"
              }
              colors <- qualitative_colors(
                length(categories), selected_palette, dark = is_dark
              )
              names(colors) <- categories
              code <- c(
                code,
                "poly_df$.patch_value <- as.character(poly_df$.patch_value)",
                paste0("patch_colors <- ", literal(colors)),
                "fill_scale <- ggplot2::scale_fill_manual(",
                "  values = patch_colors, name = fill_var,",
                "  na.value = \"grey80\", drop = TRUE",
                ")"
              )
            }
          }

          code <- c(
            code,
            paste0(
              "p <- p + ggplot2::geom_polygon(data = poly_df, ",
              "ggplot2::aes(x, y, group = patch, fill = .patch_value), "
            ),
            paste0(
              "  color = \"white\", linewidth = 0.2, alpha = ",
              literal(input$patch_alpha), ") + fill_scale"
            )
          )
        } else if (layer == "Highlighted Cells") {
          has_highlights <- !is.null(input$highlight_cats) &&
                            length(input$highlight_cats) > 0
          if (has_highlights) {
            highlight_color <- input$highlight_color
            if (is_dark && identical(highlight_color, "black")) {
              highlight_color <- "white"
            }
            code <- c(
              code, "",
              paste0("highlight_var <- ", literal(input$highlight_var)),
              paste0("highlight_cats <- ", literal(input$highlight_cats)),
              "highlight_df <- cell_df[cell_df[[highlight_var]] %in% highlight_cats, ]",
              paste0(
                "p <- p + ggplot2::geom_point(data = highlight_df, ggplot2::aes(x, y), color = ",
                literal(highlight_color), ", size = 0.2, alpha = 1)"
              )
            )
          }
        }
      }

      patch_focus <- if (identical(active_fill, ".zscore")) {
        paste("Z-score:", gene)
      } else {
        paste("Patch annotation:", active_fill)
      }
      code <- c(
        code, "",
        "p <- p + ggplot2::labs(",
        "  title = \"SpaceMosaic Interactive Spatial Map\",",
        paste0("  subtitle = ", literal(patch_focus), ","),
        "  x = \"X\", y = \"Y\"",
        ")",
        "",
        "print(p)"
      )

      paste(code, collapse = "\n")
    })

    output$download_plot_code <- shiny::downloadHandler(
      filename = function() {
        paste0("spacemosaic_plot_", format(Sys.Date(), "%Y%m%d"), ".R")
      },
      content = function(file) {
        writeLines(export_plot_code(), file, useBytes = TRUE)
      }
    )

    plot_with_legend <- shiny::reactive({
      base_plot() +
        ggplot2::theme(
          legend.position = "right",
          legend.box = "vertical",
          legend.direction = "vertical",
          legend.title = ggplot2::element_text(size = 10),
          legend.text = ggplot2::element_text(size = 8),
          legend.key.height = grid::unit(16, "pt"),
          legend.margin = ggplot2::margin(6, 6, 6, 6)
        )
    })

    # Build the legend once and reuse the same grob both for measuring and
    # drawing. Measuring the rendered grob accounts for titles, font metrics,
    # multiple guides, keys, and spacing without relying on label-count
    # approximations.
    legend_grob <- shiny::reactive({
      g <- ggplot2::ggplotGrob(plot_with_legend())
      leg_idx <- which(grepl("^guide-box", g$layout$name))

      # ggplot2 >= 3.5 stores placeholders for unused legend positions.
      # Older versions instead use a single grob named "guide-box".
      leg_idx <- leg_idx[!vapply(
        g$grobs[leg_idx],
        inherits,
        logical(1),
        what = "zeroGrob"
      )]

      if (length(leg_idx) == 0) return(NULL)
      g$grobs[[leg_idx[1]]]
    })

    legend_spec <- shiny::reactive({
      legend <- legend_grob()
      if (is.null(legend)) {
        return(list(visible = FALSE, width = 0, height = 0))
      }

      # shiny::renderPlot() interprets dimensions as pixels at the requested
      # resolution. Convert the grob's actual physical dimensions to pixels
      # and retain a small safety margin for device/font rounding.
      resolution <- 150
      padding <- 32
      width <- grid::convertWidth(
        sum(legend$widths), "inches", valueOnly = TRUE
      )
      height <- grid::convertHeight(
        sum(legend$heights), "inches", valueOnly = TRUE
      )

      list(
        visible = TRUE,
        width = max(260, ceiling(width * resolution) + padding),
        height = max(120, ceiling(height * resolution) + padding)
      )
    })

    export_dimensions <- shiny::reactive({
      legend <- legend_grob()
      if (is.null(legend)) {
        return(list(width = 9, height = 7))
      }

      legend_width <- grid::convertWidth(
        sum(legend$widths), "inches", valueOnly = TRUE
      )
      legend_height <- grid::convertHeight(
        sum(legend$heights), "inches", valueOnly = TRUE
      )

      # Reserve a roughly square area for the spatial map, then expand the
      # device when a wide or tall legend requires more room.
      list(
        width = max(9, 7.5 + legend_width + 0.5),
        height = max(7, legend_height + 1)
      )
    })

    save_current_plot <- function(file, device) {
      dimensions <- export_dimensions()
      background <- if (isTRUE(input$dark_mode)) "#161A1D" else "white"

      ggplot2::ggsave(
        filename = file,
        plot = plot_with_legend(),
        device = device,
        width = dimensions$width,
        height = dimensions$height,
        units = "in",
        dpi = if (identical(device, "png")) 300 else 72,
        bg = background,
        limitsize = FALSE
      )
    }

    output$download_plot_png <- shiny::downloadHandler(
      filename = function() {
        paste0("spacemosaic_plot_", format(Sys.Date(), "%Y%m%d"), ".png")
      },
      content = function(file) {
        save_current_plot(file, "png")
      },
      contentType = "image/png"
    )

    output$download_plot_pdf <- shiny::downloadHandler(
      filename = function() {
        paste0("spacemosaic_plot_", format(Sys.Date(), "%Y%m%d"), ".pdf")
      },
      content = function(file) {
        save_current_plot(file, "pdf")
      },
      contentType = "application/pdf"
    )

    output$spatialPlot <- shiny::renderPlot({
      # Remove legend from the main map
      base_plot() +
        ggplot2::theme(
          legend.position = "none",
          plot.margin = ggplot2::margin(8, 12, 12, 8)
        )
    }, res = 150, bg = "transparent")

    output$legend_ui <- shiny::renderUI({
      spec <- legend_spec()
      if (!spec$visible) return(NULL)

      shiny::plotOutput(
        "legendPlot",
        width = paste0(spec$width, "px"),
        height = paste0(spec$height, "px")
      )
    })

    output$legendPlot <- shiny::renderPlot({
      legend <- legend_grob()
      if (is.null(legend)) return(NULL)

      grid::grid.newpage()
      grid::grid.draw(legend)
    }, width = function() legend_spec()$width,
       height = function() legend_spec()$height,
       res = 150,
       bg = "transparent")
  }

  shiny::shinyApp(ui = ui, server = server)
}
