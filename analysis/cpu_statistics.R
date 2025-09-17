compare_performance_failure_points <- function(
  rust_file_paths,
  python_file_paths,
  method_names = NULL,
  title = "Performance & Failure Analysis: Rust vs Python",
  subtitle = "Showing maximum successful dataset size and performance",
  rust_color = "#CE422B",
  python_color = "#3776AB",
  filter_load_step = TRUE,
  base_text_size = 10,
  x_scale_power = 0.65,
  # New configuration parameters
  plot_layout = list(ncol = 4, nrow = 1),  # 4 in a row by default
  plot_width = 16,  # Increased width for 4 columns
  plot_height = 4,  # Reduced height for single row
  show_individual_legends = FALSE,
  show_global_legend = TRUE,
  legend_position = "bottom",
  x_axis_label_angle = 45,
  show_speedup_annotations = TRUE,
  show_statistics_box = TRUE,  # New flag to control statistics box
  annotation_size = 3,
  line_size = 1.2,
  point_size = 2.5,
  failure_region_alpha = 0.15,
  grid_minor = FALSE,
  custom_breaks = NULL,  # Allow custom x-axis breaks
  y_axis_on_all_plots = FALSE,  # Only show y-axis on leftmost plots
  x_axis_on_all_plots = TRUE   # Show x-axis on all plots
) {
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(scales)

  # Set method names
  if (is.null(method_names)) {
    method_names <- paste("Method", 1:length(rust_file_paths))
  }

  # Read all data
  all_data_list <- list()
  for (i in seq_along(rust_file_paths)) {
    if (file.exists(rust_file_paths[i])) {
      rust_data <- read.csv(rust_file_paths[i])
      rust_data$version <- "rust"
      rust_data$method <- method_names[i]
      all_data_list[[paste0("rust_", i)]] <- rust_data
    }

    if (file.exists(python_file_paths[i])) {
      python_data <- read.csv(python_file_paths[i])
      python_data$version <- "python"
      python_data$method <- method_names[i]
      all_data_list[[paste0("python_", i)]] <- python_data
    }
  }

  # Harmonize and combine
  all_columns <- unique(unlist(lapply(all_data_list, colnames)))
  all_data_list <- lapply(all_data_list, function(df) {
    missing_cols <- setdiff(all_columns, colnames(df))
    for (col in missing_cols) df[[col]] <- NA
    return(df[all_columns])
  })
  all_data <- do.call(rbind, all_data_list)

  # Filter
  if (filter_load_step && "test_step" %in% colnames(all_data)) {
    all_data <- all_data %>% filter(test_step != "load")
  }

  # Calculate summary statistics
  summary_data <- all_data %>%
    group_by(method, rows, version) %>%
    summarize(
      mean_time = mean(time_ms, na.rm = TRUE),
      se_time = sd(time_ms, na.rm = TRUE) / sqrt(sum(!is.na(time_ms))),
      n_runs = sum(!is.na(time_ms)),
      .groups = "drop"
    )

  # Find failure points - maximum successful size for each version
  failure_analysis <- summary_data %>%
    group_by(method, version) %>%
    summarize(
      max_successful_rows = max(rows),
      time_at_max = mean_time[which.max(rows)],
      n_sizes_tested = n_distinct(rows),
      .groups = "drop"
    ) %>%
    pivot_wider(
      names_from = version,
      values_from = c(max_successful_rows, time_at_max, n_sizes_tested)
    ) %>%
    mutate(
      # Calculate how much more data Rust can handle
      rust_advantage = max_successful_rows_rust / max_successful_rows_python,
      # Did Python fail before Rust?
      python_failed_first = max_successful_rows_python < max_successful_rows_rust,
      # Format for display
      rust_max_formatted = case_when(
        max_successful_rows_rust >= 1e6 ~ paste0(round(max_successful_rows_rust/1e6, 1), "M"),
        max_successful_rows_rust >= 1e3 & round(max_successful_rows_rust/1e3, 0) >= 1000 ~
          paste0(round(max_successful_rows_rust/1e3, 0)/1000, "M"),
        max_successful_rows_rust >= 1e3 ~ paste0(round(max_successful_rows_rust/1e3, 0), "K"),
        TRUE ~ as.character(max_successful_rows_rust)
      ),
      python_max_formatted = case_when(
        max_successful_rows_python >= 1e6 ~ paste0(round(max_successful_rows_python/1e6, 1), "M"),
        max_successful_rows_python >= 1e3 & round(max_successful_rows_python/1e3, 0) >= 1000 ~
          paste0(round(max_successful_rows_python/1e3, 0)/1000, "M"),
        max_successful_rows_python >= 1e3 ~ paste0(round(max_successful_rows_python/1e3, 0), "K"),
        TRUE ~ as.character(max_successful_rows_python)
      )
    )

  # Calculate average speedup across all common sizes
  avg_speedup_all <- summary_data %>%
    pivot_wider(
      names_from = version,
      values_from = mean_time,
      id_cols = c(method, rows)
    ) %>%
    filter(!is.na(rust) & !is.na(python)) %>%
    group_by(method) %>%
    summarize(
      avg_speedup = mean(python / rust),
      n_common_sizes = n(),
      .groups = "drop"
    )

  # Create custom power transformation
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

      # Default breaks
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

  # Create plots
  plot_list <- list()

  for (i in seq_along(unique(summary_data$method))) {
    m <- unique(summary_data$method)[i]
    method_data <- summary_data %>% filter(method == m)
    method_failure <- failure_analysis %>% filter(method == m)
    method_avg_speedup <- avg_speedup_all %>% filter(method == m)

    # Determine which plots should show y-axis labels
    show_y_axis <- if (y_axis_on_all_plots) {
      TRUE
    } else {
      # Show y-axis on leftmost plots of each row
      (i - 1) %% plot_layout$ncol == 0
    }

    # Create line plot showing performance and failure points
    p <- ggplot(method_data, aes(x = rows, y = mean_time, color = version)) +
      # Lines and points
      geom_line(linewidth = line_size, alpha = 0.9) +
      geom_point(size = point_size, alpha = 0.9) +

      # Mark failure points with vertical lines
      geom_vline(
        xintercept = method_failure$max_successful_rows_python[1],
        color = python_color,
        linetype = "dashed",
        linewidth = 1,
        alpha = 0.7
      ) +

      # Add shaded region where Python failed
      {if(method_failure$python_failed_first[1]) {
        y_min <- min(method_data$mean_time, na.rm = TRUE) * 0.5
        y_max <- max(method_data$mean_time, na.rm = TRUE) * 2

        annotate(
          "rect",
          xmin = method_failure$max_successful_rows_python[1],
          xmax = method_failure$max_successful_rows_rust[1],
          ymin = y_min,
          ymax = y_max,
          fill = "darkgreen",
          alpha = failure_region_alpha
        )
      }} +

      # Scales
      scale_x_continuous(
        trans = power_trans,
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
      scale_y_log10(labels = scales::comma) +
      scale_color_manual(
        values = c("rust" = rust_color, "python" = python_color),
        labels = c("rust" = "Rust", "python" = "Python")
      ) +

      # Labels
      labs(
        title = m,
        x = if(x_axis_on_all_plots) "Dataset Size (rows)" else NULL,
        y = if(show_y_axis) "Time (ms)" else NULL,
        color = NULL
      ) +

      # Theme
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
        panel.grid.minor = if(grid_minor) element_line() else element_blank(),
        plot.margin = margin(t = 15, r = 15, b = 15, l = 15, unit = "pt")  # Equal margins all around
      ) +
      coord_cartesian(clip = "off")  # Prevent clipping issues

    # Add speedup annotations if enabled
    if (show_speedup_annotations && show_statistics_box && nrow(method_avg_speedup) > 0) {
      rows_data <- method_data$rows
      rows_data <- rows_data[!is.na(rows_data) & rows_data > 0]

      if (length(rows_data) > 0) {
        rows_range <- range(rows_data)
        time_range <- range(method_data$mean_time, na.rm = TRUE)

        # Transform to power scale for positioning
        x_pos_transformed <- rows_range[1]^x_scale_power +
          (rows_range[2]^x_scale_power - rows_range[1]^x_scale_power) * 0.4
        x_pos <- x_pos_transformed^(1/x_scale_power)

        y_pos <- 10^(log10(time_range[1]) +
                     (log10(time_range[2]) - log10(time_range[1])) * 0.9)

        # Format the cell increase
        cells_increase <- method_failure$rust_advantage[1]
        cells_text <- if(!is.na(cells_increase) && cells_increase > 1) {
          sprintf("%.1fx more cells", cells_increase)
        } else {
          "Same capacity"
        }

        p <- p +
          annotate(
            "label",
            x = x_pos,
            y = y_pos,
            label = sprintf("Avg: %.1fx faster\n%s",
                           method_avg_speedup$avg_speedup[1],
                           cells_text),
            size = annotation_size,
            color = "gray30",
            fill = alpha("white", 0.95),
            label.padding = unit(0.4, "lines"),
            label.r = unit(0.15, "lines"),
            fontface = "bold"
          )
      }
    }

    plot_list[[m]] <- p
  }

  # Combine plots with explicit equal heights
  if (plot_layout$ncol == 1 && length(plot_list) > 1) {
    # For vertical stacking, create explicit design matrix
    design_matrix <- paste(1:length(plot_list), collapse = "\n")
    combined_plot <- wrap_plots(plot_list, design = design_matrix, heights = rep(1, length(plot_list)))
  } else {
    combined_plot <- wrap_plots(plot_list, ncol = plot_layout$ncol, nrow = plot_layout$nrow)
  }

  # Add global legend if requested
  if (show_global_legend) {
    # Create a dummy plot for the legend
    legend_data <- data.frame(
      x = c(1, 2),
      y = c(1, 2),
      version = c("rust", "python")
    )

    legend_plot <- ggplot(legend_data, aes(x = x, y = y, color = version)) +
      geom_point() +
      scale_color_manual(
        values = c("rust" = rust_color, "python" = python_color),
        labels = c("rust" = "Rust", "python" = "Python"),
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

    # Combine plots with legend
    if (legend_position == "bottom") {
      final_plot <- combined_plot / legend + plot_layout(heights = c(10, 1))
    } else if (legend_position == "right") {
      final_plot <- combined_plot | legend + plot_layout(widths = c(10, 1))
    } else {
      final_plot <- combined_plot
    }
  } else {
    final_plot <- combined_plot
  }

  # Add titles
  final_plot <- final_plot +
    plot_annotation(
      title = title,
      subtitle = subtitle,
      caption = "Lines show execution time vs dataset size. Dashed line = maximum size achieved. Green area = Rust's scalability advantage.",
      theme = theme(
        plot.title = element_text(size = base_text_size + 4, face = "bold"),
        plot.subtitle = element_text(size = base_text_size + 2),
        plot.caption = element_text(size = base_text_size - 1, color = "gray50", hjust = 0)
      )
    )

  # Print failure analysis summary
  cat("\nFailure Point Analysis:\n")
  cat(paste(rep("=", 70), collapse = ""), "\n")

  for (i in 1:nrow(failure_analysis)) {
    cat(sprintf("\n%s:\n", failure_analysis$method[i]))
    cat(sprintf("  Python: Successfully processed up to %s rows\n",
                failure_analysis$python_max_formatted[i]))
    cat(sprintf("  Rust:   Successfully processed up to %s rows\n",
                failure_analysis$rust_max_formatted[i]))

    if (failure_analysis$python_failed_first[i]) {
      cat(sprintf("  → Python FAILED beyond %s rows\n",
                  failure_analysis$python_max_formatted[i]))
      cat(sprintf("  → Rust handled %.1fx more data before failing\n",
                  failure_analysis$rust_advantage[i]))
    } else {
      cat("  → Both handled the same maximum dataset size\n")
    }

    # Add speedup info
    avg_speedup_info <- avg_speedup_all %>% filter(method == failure_analysis$method[i])
    if (nrow(avg_speedup_info) > 0) {
      cat(sprintf("  → Average speedup across %d sizes: Rust is %.1fx faster\n",
                  avg_speedup_info$n_common_sizes[1],
                  avg_speedup_info$avg_speedup[1]))
    }
  }

  # Calculate overall statistics
  cat(sprintf("\nOverall: %d/%d methods showed Python failing before Rust\n",
              sum(failure_analysis$python_failed_first),
              nrow(failure_analysis)))
  cat(sprintf("Average data handling advantage for Rust: %.1fx\n",
              mean(failure_analysis$rust_advantage[failure_analysis$python_failed_first])))

  return(list(
    plot = final_plot,
    analysis = failure_analysis,
    config = list(
      layout = plot_layout,
      dimensions = c(width = plot_width, height = plot_height)
    )
  ))
}

# Helper function to extract legend (you might need cowplot for this)
get_legend <- function(plot) {
  if (requireNamespace("cowplot", quietly = TRUE)) {
    return(cowplot::get_legend(plot))
  } else {
    warning("cowplot package not available. Install it for legend extraction.")
    return(NULL)
  }
}

# Geometric mean helper function
geometric_mean <- function(x) {
  exp(mean(log(x)))
}

# Usage examples:

# Example 1: 4 plots in a row (default)
# rust_files <- c(
#   "/local/06-24_single-rust-project/single-test-comparison/statistics/rust/filtering_performance.csv",
#   "/local/06-24_single-rust-project/single-test-comparison/statistics/rust/hvg_performance.csv",
#   "/local/06-24_single-rust-project/single-test-comparison/statistics/rust/deg_performance.csv",
#   "/local/06-24_single-rust-project/single-test-comparison/statistics/rust/pca_performance.csv",
#   "/local/06-24_single-rust-project/single-test-comparison/statistics/rust/normalization_performance.csv",
#   "/local/06-24_single-rust-project/single-test-comparison/statistics/rust/nn_performance_t.csv"
# )

# python_files <- c(
#   "/local/06-24_single-rust-project/single-test-comparison/statistics/python/filtering_performance.csv",
#   "/local/06-24_single-rust-project/single-test-comparison/statistics/python/hvg_performance.csv",
#   "/local/06-24_single-rust-project/single-test-comparison/statistics/python/deg_performance.csv",
#   "/local/06-24_single-rust-project/single-test-comparison/statistics/python/pca_performance.csv",
#   "/local/06-24_single-rust-project/single-test-comparison/statistics/python/normalization_performance.csv",
#   "/local/06-24_single-rust-project/single-test-comparison/statistics/python/nn_performance.csv"
# )

# result_clean <- compare_performance_failure_points(
#   rust_file_paths = rust_files,
#   python_file_paths = python_files,
#   method_names = c("Filtering", "HVG", "DEG", "PCA", "Normalization", "KNN"),
#   plot_layout = list(ncol = 1, nrow = 7),  # Increased to 7 rows to accommodate all methods
#   show_speedup_annotations = FALSE,  # Disable speedup annotations entirely
#   show_statistics_box = FALSE,       # Disable statistics box for cleaner horizontal layout
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

# ggsave(
#   "performance_horizontal_layout.png",  # Updated filename to reflect layout
#   plot = result_clean$plot,
#   width = result_clean$config$dimensions["width"],
#   height = result_clean$config$dimensions["height"],
#   dpi = 300,
#   bg = "white"
# )
