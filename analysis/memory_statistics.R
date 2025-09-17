# Memory Consumption Comparison: Multiple Task Analysis
# Adapted to work with your existing memory plotting infrastructure

library(ggplot2)
library(dplyr)
library(tidyr)  # Added this for pivot_wider
library(readr)
library(scales)
library(patchwork)  # Changed from gridExtra to patchwork
library(cowplot)
library(stringr)
library(grid)

# Your existing cell extraction function (keeping it identical)
extract_cell_count <- function(dataset_name) {
  cell_mapping <- c(
    "tahoe_10K_cells" = 10000,
    "tahoe_50K_cells" = 50000,
    "tahoe_100K_cells" = 100000,
    "tahoe_250K_cells" = 250000,
    "tahoe_500K_cells" = 500000,
    "tahoe_750K_cells" = 750000,
    "tahoe_1M_cells" = 1000000,
    "tahoe_2M_cells" = 2000000,
    "tahoe_5M_cells" = 5000000,
    "tahoe_7M_cells" = 7500000,
    "tahoe_10M_cells" = 10000000,
    "tahoe_12M_cells" = 12500000,
    "tahoe_15M_cells" = 15000000,
    "tahoe_streaming_17M_cells" = 17500000,
    "tahoe_streaming_20M_cells" = 20000000,
    "tahoe_22M_cells" = 22500000,
    "tahoe_25M_cells" = 25000000,
    "tahoe_27M_cells" = 27500000,
    "tahoe_30M_cells" = 30000000
  )

  pattern <- "tahoe_(?:streaming_)?\\d+(?:\\.\\d+)?[KM]_cells"
  match <- str_extract(dataset_name, pattern)

  if (!is.na(match) && match %in% names(cell_mapping)) {
    return(cell_mapping[[match]])
  }

  return(NA)
}

# Enhanced memory data processing function
process_memory_data_enhanced <- function(file_path, language_name, task_name) {
  if (!file.exists(file_path)) {
    warning(sprintf("File not found: %s", file_path))
    return(NULL)
  }

  tryCatch({
    data <- read_csv(file_path, show_col_types = FALSE)

    processed_data <- data %>%
      mutate(
        cells = sapply(dataset, extract_cell_count),
        memory_mb = max_memory_kb / 1024,
        memory_gb = max_memory_kb / (1024 * 1024),
        memory_per_cell_kb = max_memory_kb / cells,
        memory_per_cell_mb = memory_mb / cells,
        language = language_name,
        task = task_name
      ) %>%
      filter(!is.na(cells)) %>%
      arrange(cells)

    return(processed_data)
  }, error = function(e) {
    warning(sprintf("Error reading %s: %s", file_path, e$message))
    return(NULL)
  })
}

# Main comparison function for memory usage
compare_memory_consumption <- function(
  rust_file_paths,
  python_file_paths,
  method_names = NULL,
  title = "Memory Consumption Analysis: Rust vs Python",
  subtitle = "Peak memory usage vs dataset size",
  rust_color = "#CE422B",
  python_color = "#3776AB",
  filter_load_step = TRUE,  # Match plotting.R (not used in memory but keep for compatibility)
  base_text_size = 10,
  x_scale_power = 0.65,
  # New configuration parameters (matching plotting.R exactly)
  plot_layout = list(ncol = 1, nrow = 7),  # Increased to 7 rows to accommodate all methods
  plot_width = 4,          # Reduced width for vertical stacking
  plot_height = 21,        # Increased height for 7 stacked plots (3 per plot)
  show_individual_legends = FALSE,
  show_global_legend = FALSE,  # Disable legend (matching plotting.R)
  legend_position = "bottom",
  x_axis_label_angle = 45,
  show_speedup_annotations = FALSE,  # Disable annotations (matching plotting.R)
  show_statistics_box = FALSE,  # Disable statistics box for cleaner horizontal layout
  annotation_size = 2.8,   # Slightly larger annotations
  line_size = 1.1,         # Slightly thicker lines
  point_size = 2.2,        # Slightly larger points
  failure_region_alpha = 0.15,
  grid_minor = TRUE,       # Enable minor grid lines for better readability
  custom_breaks = NULL,    # Allow custom x-axis breaks
  y_axis_on_all_plots = TRUE,  # Show y-axis on all plots for vertical layout
  x_axis_on_all_plots = TRUE,  # Show x-axis on all plots
  # Memory-specific parameters (keeping for backwards compatibility)
  point_alpha = 0.9,
  smooth_span = 0.75,
  smooth_alpha = 0.2,
  filter_outliers = TRUE,
  outlier_threshold_gb = 0.1,
  dpi = 300
) {

  # Set method names
  if (is.null(method_names)) {
    method_names <- paste("Method", 1:length(rust_file_paths))
  }

  # Read and process all data
  all_data_list <- list()
  successful_tasks <- 0

  cat("Processing memory data files...\n")

  for (i in seq_along(rust_file_paths)) {
    task_name <- method_names[i]
    cat(sprintf("Processing task %d/%d: %s\n", i, length(rust_file_paths), task_name))

    # Process Rust data
    rust_data <- process_memory_data_enhanced(rust_file_paths[i], "Rust", task_name)
    if (!is.null(rust_data)) {
      all_data_list[[paste0("rust_", i)]] <- rust_data
    }

    # Process Python data
    python_data <- process_memory_data_enhanced(python_file_paths[i], "Python", task_name)
    if (!is.null(python_data)) {
      all_data_list[[paste0("python_", i)]] <- python_data
    }

    # Check if we have data for this task
    if (!is.null(rust_data) || !is.null(python_data)) {
      successful_tasks <- successful_tasks + 1
    }
  }

  if (length(all_data_list) == 0) {
    stop("No valid data files found. Please check your file paths.")
  }

  # Combine all data
  all_data <- bind_rows(all_data_list)

  # Filter outliers if requested
  if (filter_outliers) {
    outlier_filter <- !(all_data$cells == 100000 & all_data$memory_gb < outlier_threshold_gb)
    all_data <- all_data[outlier_filter, ]

    n_filtered <- sum(!outlier_filter)
    if (n_filtered > 0) {
      cat(sprintf("Filtered %d outlier data points\n", n_filtered))
    }
  }

  # Memory efficiency analysis
  memory_efficiency <- all_data %>%
    group_by(task, language) %>%
    summarise(
      max_cells = max(cells),
      memory_at_max = memory_gb[which.max(cells)],
      avg_memory_per_cell_kb = mean(memory_per_cell_kb),
      min_memory_gb = min(memory_gb),
      max_memory_gb = max(memory_gb),
      n_datapoints = n(),
      .groups = "drop"
    )

  # Failure analysis using base R approach (no tidyr dependency)
  temp_data <- memory_efficiency %>%
    select(task, language, max_cells, memory_at_max)

  # Split by language
  rust_data <- temp_data %>%
    filter(language == "Rust") %>%
    select(task, max_cells_Rust = max_cells, memory_at_max_Rust = memory_at_max)

  python_data <- temp_data %>%
    filter(language == "Python") %>%
    select(task, max_cells_Python = max_cells, memory_at_max_Python = memory_at_max)

  # Merge and calculate ratios
  failure_analysis <- merge(rust_data, python_data, by = "task", all = TRUE) %>%
    mutate(
      # Handle potential NA values
      max_cells_Rust = ifelse(is.na(max_cells_Rust), 0, max_cells_Rust),
      max_cells_Python = ifelse(is.na(max_cells_Python), 0, max_cells_Python),
      memory_at_max_Rust = ifelse(is.na(memory_at_max_Rust), 0, memory_at_max_Rust),
      memory_at_max_Python = ifelse(is.na(memory_at_max_Python), 0, memory_at_max_Python),

      # Calculate ratios
      rust_advantage_cells = ifelse(max_cells_Python > 0, max_cells_Rust / max_cells_Python, NA),
      python_failed_first = max_cells_Python < max_cells_Rust,
      memory_efficiency_ratio = ifelse(memory_at_max_Rust > 0, memory_at_max_Python / memory_at_max_Rust, NA)
    )

  # Define colors
  colors <- c("Rust" = rust_color, "Python" = python_color)

  # Create custom power transformation (exactly like performance plot)
  power_trans <- trans_new(
    "power",
    transform = function(x) {
      x_safe <- pmax(x, 1e-10)
      x_safe^x_scale_power
    },
    inverse = function(x) {
      x_safe <- pmax(x, 1e-10)
      x_safe^(1/x_scale_power)
    },
    breaks = function(x) {
      if (!is.null(custom_breaks)) {
        return(custom_breaks)
      }

      # Default breaks (same as performance plot)
      desired_breaks <- c(
        10000, 500000, 1000000, 5000000, 10000000, 15000000, 20000000
      )

      if (length(x) > 0) {
        x <- x[!is.na(x) & x > 0]
        if (length(x) > 0) {
          range_vals <- range(x)
          desired_breaks <- desired_breaks[desired_breaks >= range_vals[1] * 0.8 &
                                         desired_breaks <= range_vals[2] * 1.2]
        }
      }

      return(desired_breaks)
    },
    domain = c(1e-10, Inf)
  )

  # Create individual plots for each task
  plot_list <- list()

  for (i in seq_along(unique(all_data$task))) {
    task_name <- unique(all_data$task)[i]
    task_data <- all_data %>% filter(task == task_name)
    task_failure <- failure_analysis %>% filter(task == task_name)
    task_efficiency <- memory_efficiency %>% filter(task == task_name)

    cat(sprintf("Creating plot for: %s\n", task_name))

    # Determine which plots should show y-axis labels (match performance plot logic)
    show_y_axis <- if (y_axis_on_all_plots) {
      TRUE
    } else {
      # Show y-axis on leftmost plots of each row
      (i - 1) %% plot_layout$ncol == 0
    }

    # Create the memory consumption plot (matching performance plot style)
    p <- ggplot(task_data, aes(x = cells, y = memory_gb, color = language)) +
      # Lines and points (matching performance plot)
      geom_line(linewidth = line_size, alpha = point_alpha) +
      geom_point(size = point_size, alpha = point_alpha) +

      # Add failure point visualization where Rust outperforms Python
      {if(show_statistics_box && nrow(task_failure) > 0 &&
          !is.na(task_failure$python_failed_first[1]) &&
          task_failure$python_failed_first[1]) {

        python_max_cells <- task_failure$max_cells_Python[1]
        rust_max_cells <- task_failure$max_cells_Rust[1]

        # Only add visualization if there's a meaningful difference
        if (rust_max_cells > python_max_cells && python_max_cells > 0) {

          # Get y-axis range for shading (match performance plot method)
          y_min <- min(task_data$memory_gb, na.rm = TRUE) * 0.5
          y_max <- max(task_data$memory_gb, na.rm = TRUE) * 2

          list(
            # Vertical line at Python failure point
            geom_vline(
              xintercept = python_max_cells,
              color = python_color,
              linetype = "dashed",
              linewidth = 1,
              alpha = 0.7
            ),
            # Green shaded region where Rust continues but Python failed
            annotate(
              "rect",
              xmin = python_max_cells,
              xmax = rust_max_cells,
              ymin = y_min,
              ymax = y_max,
              fill = "darkgreen",
              alpha = failure_region_alpha
            )
          )
        }
      }} +

      # Scales (matching performance plot formatting exactly)
      scale_x_continuous(
        trans = power_trans,  # Add power transformation
        labels = function(x) {
          sapply(x, function(val) {
            if (is.na(val)) return(NA)

            rounded_val <- round(val, 0)

            if (rounded_val >= 1e6) {
              if (rounded_val %% 1e6 == 0) {
                paste0(rounded_val/1e6, "M")
              } else {
                paste0(round(rounded_val/1e6, 1), "M")
              }
            } else if (rounded_val >= 1e3) {
              thousands <- round(rounded_val/1e3, 0)
              if (thousands >= 1000) {
                paste0(thousands/1000, "M")
              } else {
                paste0(thousands, "K")
              }
            } else {
              as.character(rounded_val)
            }
          })
        },
        breaks = if(is.null(custom_breaks)) c(10000, 500000, 1000000, 5000000, 10000000, 15000000, 20000000, 25000000, 30000000) else custom_breaks
      ) +
      scale_y_log10(labels = function(x) {
        # Custom formatter to avoid mixing 1000 and 1000.0
        formatted <- scales::comma(x, accuracy = 1)
        # Remove unnecessary .0 from whole numbers
        formatted <- gsub("\\.0$", "", formatted)
        return(formatted)
      }) +
      scale_color_manual(
        values = c("Rust" = rust_color, "Python" = python_color),
        labels = c("Rust" = "Rust", "Python" = "Python")
      ) +

      # Labels (matching performance plot style)
      labs(
        title = task_name,
        x = if(x_axis_on_all_plots) "Dataset Size (rows)" else NULL,
        y = if(show_y_axis) "Memory (GB)" else NULL,  # Keep memory-specific label
        color = NULL
      ) +

      # Theme (exactly matching performance plot with optimized spacing)
      theme_minimal(base_size = base_text_size) +
      theme(
        plot.title = element_text(face = "bold", hjust = 0.5, size = base_text_size + 1),
        axis.text = element_text(size = base_text_size - 1),
        axis.text.x = element_text(angle = x_axis_label_angle, hjust = 1, vjust = 1),
        axis.title = element_text(size = base_text_size),
        axis.title.x = if(x_axis_on_all_plots) element_text() else element_blank(),
        axis.title.y = if(show_y_axis) element_text() else element_blank(),
        legend.position = if(show_individual_legends) "bottom" else "none",
        legend.background = element_rect(fill = alpha("white", 0.9)),
        legend.text = element_text(size = base_text_size - 1),
        panel.grid.minor = if(grid_minor) element_line() else element_blank(),  # Match performance plot
        plot.margin = margin(t = 15, r = 15, b = 15, l = 15, unit = "pt")  # Equal margins all around
      ) +
      coord_cartesian(clip = "off")  # Prevent clipping issues

    # Add memory efficiency annotation if requested
    if (show_speedup_annotations && nrow(task_efficiency) >= 2) {
      rust_eff <- task_efficiency %>% filter(language == "Rust")
      python_eff <- task_efficiency %>% filter(language == "Python")

      if (nrow(rust_eff) > 0 && nrow(python_eff) > 0) {
        efficiency_ratio <- python_eff$avg_memory_per_cell_kb[1] / rust_eff$avg_memory_per_cell_kb[1]
        max_cell_advantage <- rust_eff$max_cells[1] / python_eff$max_cells[1]

        # Position annotation
        x_pos <- max(task_data$cells) * 0.3
        y_pos <- max(task_data$memory_gb) * 0.9

        annotation_text <- sprintf("Rust: %.1fx more efficient\n%.1fx larger datasets",
                                 efficiency_ratio, max_cell_advantage)

        p <- p + annotate(
          "label",
          x = x_pos,
          y = y_pos,
          label = annotation_text,
          size = annotation_size,  # Use the matching parameter name
          color = "gray30",
          fill = alpha("white", 0.95),
          label.padding = unit(0.4, "lines"),
          label.r = unit(0.15, "lines"),
          fontface = "bold"
        )
      }
    }

    plot_list[[length(plot_list) + 1]] <- p

    # Print brief summary for this task
    summary_stats <- task_data %>%
      group_by(language) %>%
      summarise(
        n = n(),
        min_gb = round(min(memory_gb), 2),
        max_gb = round(max(memory_gb), 2),
        avg_efficiency_kb_per_cell = round(mean(memory_per_cell_kb), 1),
        .groups = 'drop'
      )
    print(summary_stats)
  }

  # Create combined grid using patchwork (exactly like performance plots)
  if (length(plot_list) > 0) {
    cat(sprintf("\nCreating grid with %d plots...\n", length(plot_list)))

    # Combine plots using patchwork with explicit equal heights for vertical stacking
    if (plot_layout$ncol == 1 && length(plot_list) > 1) {
      # For vertical stacking, create explicit design matrix
      design_matrix <- paste(1:length(plot_list), collapse = "\n")
      combined_plot <- wrap_plots(plot_list, design = design_matrix, heights = rep(1, length(plot_list)))
    } else {
      combined_plot <- wrap_plots(plot_list, ncol = plot_layout$ncol, nrow = plot_layout$nrow)
    }

    # Add global legend if requested (matching performance plot method)
    if (show_global_legend) {
      # Create a dummy plot for the legend
      legend_data <- data.frame(
        x = c(1, 2),
        y = c(1, 2),
        language = c("Rust", "Python")
      )

      legend_plot <- ggplot(legend_data, aes(x = x, y = y, color = language)) +
        geom_point() +
        scale_color_manual(
          values = c("Rust" = rust_color, "Python" = python_color),
          labels = c("Rust" = "Rust", "Python" = "Python"),
          name = ""
        ) +
        theme_void() +
        theme(
          legend.position = legend_position,
          legend.text = element_text(size = base_text_size),
          legend.margin = margin(t = 10)
        )

      # Extract legend
      legend <- get_legend(legend_plot)

      # Combine plots with legend (matching performance plot logic exactly)
      if (legend_position == "bottom") {
        final_plot <- combined_plot / legend + plot_layout(heights = c(10, 1))  # Match plotting.R exactly
      } else if (legend_position == "right") {
        final_plot <- combined_plot | legend + plot_layout(widths = c(10, 1))   # Match plotting.R exactly
      } else {
        final_plot <- combined_plot
      }
    } else {
      final_plot <- combined_plot
    }

    # Add titles (exactly matching performance plot with optimized spacing)
    final_plot <- final_plot +
      plot_annotation(
        title = title,
        subtitle = subtitle,
        caption = "Lines show memory usage vs dataset size. Dashed line = maximum size achieved. Green area = Rust's scalability advantage.",
        theme = theme(
          plot.title = element_text(size = base_text_size + 4, face = "bold", margin = margin(b = 5)),
          plot.subtitle = element_text(size = base_text_size + 2, margin = margin(b = 10)),
          plot.caption = element_text(size = base_text_size - 1, color = "gray50", hjust = 0, margin = margin(t = 5))
        )
      )

    # Print comprehensive analysis
    cat("\nMemory Consumption Analysis:\n")
    cat(paste(rep("=", 70), collapse = ""), "\n")

    for (i in 1:nrow(failure_analysis)) {
      task_name <- failure_analysis$task[i]
      cat(sprintf("\n%s:\n", task_name))

      # Format cell counts
      python_max <- failure_analysis$max_cells_Python[i]
      rust_max <- failure_analysis$max_cells_Rust[i]

      python_formatted <- case_when(
        python_max >= 1e6 ~ paste0(round(python_max/1e6, 1), "M"),
        python_max >= 1e3 ~ paste0(round(python_max/1e3, 0), "K"),
        TRUE ~ as.character(python_max)
      )

      rust_formatted <- case_when(
        rust_max >= 1e6 ~ paste0(round(rust_max/1e6, 1), "M"),
        rust_max >= 1e3 ~ paste0(round(rust_max/1e3, 0), "K"),
        TRUE ~ as.character(rust_max)
      )

      cat(sprintf("  Dataset Capacity:\n"))
      cat(sprintf("    Python: Up to %s cells (%.1f GB)\n",
                  python_formatted, failure_analysis$memory_at_max_Python[i]))
      cat(sprintf("    Rust:   Up to %s cells (%.1f GB)\n",
                  rust_formatted, failure_analysis$memory_at_max_Rust[i]))

      if (failure_analysis$python_failed_first[i]) {
        cat(sprintf("  → Python FAILED beyond %s cells\n", python_formatted))
        cat(sprintf("  → Rust processes %.1fx more cells\n",
                    failure_analysis$rust_advantage_cells[i]))
        cat(sprintf("  → GREEN REGION shows Rust's scalability advantage\n"))
      } else {
        cat("  → Both achieved similar maximum dataset sizes\n")
      }

      if (!is.na(failure_analysis$memory_efficiency_ratio[i])) {
        cat(sprintf("  → Python uses %.1fx more memory at max capacity\n",
                    failure_analysis$memory_efficiency_ratio[i]))
      }
    }

    # Overall statistics
    avg_rust_advantage <- mean(failure_analysis$rust_advantage_cells[failure_analysis$python_failed_first], na.rm = TRUE)
    avg_memory_efficiency <- mean(failure_analysis$memory_efficiency_ratio, na.rm = TRUE)

    cat(sprintf("\nOverall Summary:\n"))
    cat(sprintf("  Tasks where Rust handles larger datasets: %d/%d\n",
                sum(failure_analysis$python_failed_first, na.rm = TRUE), nrow(failure_analysis)))
    cat(sprintf("  Average cell capacity advantage: %.1fx\n", avg_rust_advantage))
    cat(sprintf("  Average memory efficiency advantage: %.1fx\n", avg_memory_efficiency))
    cat(sprintf("  Successfully processed %d tasks\n", successful_tasks))

    # Return results (matching performance plot structure)
    return(list(
      plot = final_plot,
      analysis = failure_analysis,
      efficiency_analysis = memory_efficiency,
      config = list(
        layout = plot_layout,
        dimensions = c(width = plot_width, height = plot_height)
      )
    ))
  } else {
    stop("No plots were created. Please check your data.")
  }
}

# Helper function to extract legend (matching performance plot)
get_legend <- function(plot) {
  if (requireNamespace("cowplot", quietly = TRUE)) {
    return(cowplot::get_legend(plot))
  } else {
    warning("cowplot package not available. Install it for legend extraction.")
    return(NULL)
  }
}

# # Usage example with your memory data files (matching performance plot style)
# rust_memory_files <- c(
#   "/local/06-24_single-rust-project/single-test-comparison/statistics/rust/filtering_memory.csv",
#   "/local/06-24_single-rust-project/single-test-comparison/statistics/rust/hvg_memory.csv",
#   "/local/06-24_single-rust-project/single-test-comparison/statistics/rust/deg_memory.csv",
#   "/local/06-24_single-rust-project/single-test-comparison/statistics/rust/pca_memory.csv",
#   "/local/06-24_single-rust-project/single-test-comparison/statistics/rust/normalization_memory.csv",
#   "/local/06-24_single-rust-project/single-test-comparison/statistics/rust/nn_memory_t.csv"
# )

# python_memory_files <- c(
#   "/local/06-24_single-rust-project/single-test-comparison/statistics/python/filtering_memory.csv",
#   "/local/06-24_single-rust-project/single-test-comparison/statistics/python/hvg_memory.csv",
#   "/local/06-24_single-rust-project/single-test-comparison/statistics/python/deg_memory.csv",
#   "/local/06-24_single-rust-project/single-test-comparison/statistics/python/pca_memory.csv",
#   "/local/06-24_single-rust-project/single-test-comparison/statistics/python/normalization_memory.csv",
#   "/local/06-24_single-rust-project/single-test-comparison/statistics/python/nn_memory.csv"
# )

# # Create memory consumption comparison (matching performance plot call exactly)
# memory_result <- compare_memory_consumption(
#   rust_file_paths = rust_memory_files,
#   python_file_paths = python_memory_files,
#   method_names = c("Filtering", "HVG", "DEG", "PCA", "Normalization", "KNN"),
#   plot_layout = list(ncol = 1, nrow = 7),  # Increased to 7 rows to accommodate all methods
#   show_speedup_annotations = FALSE,  # Disable speedup annotations entirely
#   show_statistics_box = TRUE,       # Disable statistics box for cleaner horizontal layout
#   plot_width = 4,                    # Reduced width for vertical stacking
#   plot_height = 21,                  # Increased height for 7 stacked plots (3 per plot)
#   y_axis_on_all_plots = TRUE,        # Show y-axis on all plots for vertical layout
#   x_axis_label_angle = 45,           # Angled labels for better readability
#   base_text_size = 10,               # Slightly larger text for better readability
#   annotation_size = 2.8,             # Slightly larger annotations
#   point_size = 2.2,                  # Slightly larger points
#   line_size = 1.1,                   # Slightly thicker lines
#   grid_minor = TRUE                  # Enable minor grid lines for better readability
# )

# Save the plot (matching performance plot save)
# if (!is.null(memory_result)) {
#   ggsave(
#     "memory_consumption_comparison.png",
#     plot = memory_result$plot,
#     width = memory_result$config$dimensions["width"],
#     height = memory_result$config$dimensions["height"],
#     dpi = 300,
#     bg = "white"
#   )

#   cat("\nSaved memory consumption comparison plot!\n")
# }
