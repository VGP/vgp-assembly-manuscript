# Genome Assembly Bubbles
# -----------------------
# by Diego De Panis
# ERGA Sequencing and Assembly Committee

# This R code will produce a plot like the one we did for the ERGA Pilot paper (Fig.4)
# doi.org/10.1038/s44185-024-00054-6

# The plot can be used to show relationships between multiple variables simultaneously:
# - Bubble Position (x, y coord)
# - Bubble Size
# - Bubble Colour (can be categorical or continuous)
# - Bubble Shape (filled/unfilled for binary categories)
#
# Required:
# - x_var: Numeric, x-axis. Can be log-transformed
# - y_var: Numeric, y-axis. Can be log-transformed
# - size_var: Numeric, bubble size
#
# Optional:
# - shape_var: Binary/dichotomic, for filled/unfilled bubbles
# - colour_var: Can be:
#     * Factor/character for categorical (assigns palette based on number of levels)
#     * Numeric for continuous colour gradient
# - label_var: Variable to use for labeling points when they meet conditions
#
# Also:
# - Default colours, but custom palettes can be specified
# - Reference lines can be added to both axes
# - Logarithmic transformation can be set
# - Legend titles and breaks can be custom
# - Point labels can be added based on custom conditions

library(tidyverse)
library(scales)
library(RColorBrewer)

GenomeAssemblyBubbles <- function(data,
                                 x_var,
                                 y_var,
                                 size_var,
                                 shape_var = NULL,
                                 colour_var = NULL,
                                 label_var = NULL,
                                 x_log = FALSE,
                                 y_log = FALSE,
                                 size_range = c(4, 16),
                                 size_breaks = waiver(),
                                 size_name = "Size",
                                 colour_palette = NULL,
                                 h_line = NULL,
                                 v_line = NULL,
                                 label_threshold_fn = NULL) {
  
  # Default palettes
  default_cat_palettes <- list(
    "2" = c("#1f77b4", "#ff7f0e"),
    "3" = c("#1f77b4", "#ff7f0e", "#2ca02c"),
    "4" = c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728"),
    "5" = c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd"),
    "6+" = brewer.pal(8, "Set2")
  )
  
  default_continuous_palette <- "viridis"
  
  # Base mapping
  plot_mapping <- aes(
    x = if(x_log) log10(.data[[x_var]]) else .data[[x_var]],
    y = if(y_log) log10(.data[[y_var]]) else .data[[y_var]],
    size = .data[[size_var]],
    colour = if(!is.null(colour_var)) .data[[colour_var]] else NULL
  )
  
  # Add shape and fill if shape_var provided
  if (!is.null(shape_var)) {
    plot_mapping$shape <- as.name(shape_var)
    plot_mapping$fill <- as.name(colour_var)
  }
  
  # Base plot
  p <- ggplot(data, plot_mapping)
  
  # Add points with different parameters for filled/unfilled
  if (!is.null(shape_var)) {
    # Split data for different point types
    shape_col <- data[[shape_var]]
    levels_shape <- levels(shape_col)
    
    # For unfilled points
    unfilled_data <- data[shape_col == levels_shape[1], ]
    filled_data <- data[shape_col == levels_shape[2], ]
    
    p <- p + 
      geom_point(data = unfilled_data, stroke = 1.1) +
      geom_point(data = filled_data, alpha = 0.8, stroke = 1.1)
  } else {
    p <- p + geom_point(alpha = 0.8)
  }
  
  #Style
  p <- p + 
    theme_minimal() +
    theme(axis.text = element_text(size = 14),
          axis.title = element_text(size = 15),
          legend.position = "right",
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 10))
  
  # Set scales
  if (!is.null(shape_var)) {
    p <- p + scale_shape_manual(
      values = c(1, 21),  # unfilled circle and filled circle with border
      guide = guide_legend(override.aes = list(
        fill = c("white", "black"),  # Make second point appear filled in legend
        size = 3
      ))
    )
  }
  
  # Set size scale
  p <- p + scale_size_continuous(
    range = size_range,
    breaks = size_breaks,
    name = size_name
  )
  
  if (!is.null(colour_var)) {
    # Check if colour variable is discrete or cont
    if (is.factor(data[[colour_var]]) || is.character(data[[colour_var]])) {
      # Get number of unique values
      unique_values <- unique(data[[colour_var]])
      n_categories <- length(unique_values)
      
      # Select default categorical palette based on categories
      if (is.null(colour_palette)) {
        if (n_categories <= 5) {
          colour_palette <- default_cat_palettes[[as.character(n_categories)]]
        } else {
          colour_palette <- default_cat_palettes[["6+"]][1:min(n_categories, 8)]
        }
      }
      
      # Apply categorical colours
      p <- p + scale_colour_manual(
        values = colour_palette,
        breaks = unique_values,
        name = colour_var
      )
      if (!is.null(shape_var)) {
        p <- p + scale_fill_manual(
          values = colour_palette,
          breaks = unique_values,
          name = colour_var
        )
      }
    } else {
      # Apply continuous colours
      if (is.null(colour_palette)) {
        # Use default cont palette
        p <- p + scale_colour_viridis_c(name = colour_var)
        if (!is.null(shape_var)) {
          p <- p + scale_fill_viridis_c(name = colour_var)
        }
      } else {
        # Use custom cont palette (if provided)
        p <- p + scale_colour_gradientn(colours = colour_palette, name = colour_var)
        if (!is.null(shape_var)) {
          p <- p + scale_fill_gradientn(colours = colour_palette, name = colour_var)
        }
      }
    }
  }
  
  # Add reference lines if set
  if (!is.null(h_line)) {
    p <- p + geom_hline(yintercept = h_line, linetype = "dashed", colour = "black")
  }
  if (!is.null(v_line)) {
    p <- p + geom_vline(xintercept = v_line, linetype = "dashed", colour = "black")
  }
  
  # Add labels if set
  if (!is.null(label_var) && !is.null(label_threshold_fn)) {
    label_data <- data[label_threshold_fn(data), ]
    p <- p + geom_text(
      data = label_data,
      aes(label = .data[[label_var]]),
      size = 3.75,
      show.legend = FALSE,
      nudge_y = 0.02,
      vjust = -2.5,
      colour = "black"
    )
  }
  
  # Format axis labels
  if (x_log) {
    p <- p + xlab(paste0("log10(", x_var, ")"))
  } else {
    p <- p + xlab(x_var)
  }
  if (y_log) {
    p <- p + ylab(paste0("log10(", y_var, ")"))
  } else {
    p <- p + ylab(y_var)
  }
  
  return(p)
}


# EXAMPLE using mtcars dataset
mtcars_test <- mtcars
mtcars_test$am <- as.factor(ifelse(mtcars_test$am == 0, "Automatic", "Manual"))
mtcars_test$car_name <- rownames(mtcars_test)


mtcars_test$cyl <- as.factor(mtcars_test$cyl)

showcase_plot1 <- GenomeAssemblyBubbles(
  data = mtcars_test,
  x_var = "hp",
  y_var = "qsec",
  size_var = "disp",
  shape_var = "am",
  colour_var = "cyl",
  label_var = "car_name",
  x_log = TRUE,
  y_log = FALSE,
  size_name = "Engine disp",
  size_breaks = seq(100, 500, by = 100),
  h_line = 16,
  v_line = log10(200),
  label_threshold_fn = function(data) {data$qsec < 15} # use {TRUE} for print all
)
print(showcase_plot1)


mtcars_test$cyl <- as.numeric(as.character(mtcars_test$cyl))
showcase_plot2 <- GenomeAssemblyBubbles(
  data = mtcars_test,
  x_var = "wt",
  y_var = "mpg",
  size_var = "disp",
  colour_var = "cyl",
  shape_var = "am",
  size_name = "disp",
  size_breaks = seq(100, 500, by = 100)
)
print(showcase_plot2)
