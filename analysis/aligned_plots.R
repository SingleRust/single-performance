# Align Memory and Performance Statistics
# This script compares and aligns statistics from both memory and performance analyses

library(dplyr)
library(tidyr)
library(readr)

# Source the plotting functions
source("cpu_statistics.R")
source("memory_statistics.R")

setwd("/local/06-24_single-rust-project/single-performance/analysis-results")

# Function to align and compare memory and performance statistics
align_memory_performance_statistics <- function(
  rust_performance_files,
  python_performance_files,
  rust_memory_files,
  python_memory_files,
  method_names = NULL,
  output_file = "aligned_statistics.csv",
  detailed_output = TRUE,
  generate_plots = TRUE,
  plot_output_dir = "aligned_plots"
) {

  if (is.null(method_names)) {
    method_names <- paste("Method", 1:length(rust_performance_files))
  }

  cat("Generating performance statistics...\n")

  # Generate performance statistics with explicit parameters for consistency
  performance_result <- compare_performance_failure_points(
    rust_file_paths = rust_performance_files,
    python_file_paths = python_performance_files,
    method_names = method_names,
    plot_layout = list(ncol = 1, nrow = length(method_names)),
    show_speedup_annotations = FALSE,
    show_statistics_box = TRUE,
    plot_width = 4,
    plot_height = 21,
    y_axis_on_all_plots = TRUE,
    show_global_legend = FALSE,
    base_text_size = 10,
    x_scale_power = 0.65,
    rust_color = "#CE422B",
    python_color = "#3776AB"
  )

  cat("Performance analysis completed.\n")
  cat(sprintf("Performance plots generated: %d\n", length(performance_result$plot$layers) > 0))

  cat("Generating memory statistics...\n")

  # Generate memory statistics with identical parameters
  memory_result <- compare_memory_consumption(
    rust_file_paths = rust_memory_files,
    python_file_paths = python_memory_files,
    method_names = method_names,
    plot_layout = list(ncol = 1, nrow = length(method_names)),
    show_speedup_annotations = FALSE,
    show_statistics_box = TRUE,
    plot_width = 4,
    plot_height = 21,
    y_axis_on_all_plots = TRUE,
    show_global_legend = FALSE,
    base_text_size = 10,
    x_scale_power = 0.65,
    rust_color = "#CE422B",
    python_color = "#3776AB"
  )

  cat("Memory analysis completed.\n")
  cat(sprintf("Memory plots generated: %d\n", length(memory_result$plot$layers) > 0))

  # Extract analysis data
  perf_analysis <- performance_result$analysis
  mem_analysis <- memory_result$analysis
  mem_efficiency <- memory_result$efficiency_analysis

  # Align the statistics by method
  aligned_stats <- perf_analysis %>%
    select(
      method,
      perf_max_rows_python = max_successful_rows_python,
      perf_max_rows_rust = max_successful_rows_rust,
      perf_rust_advantage = rust_advantage,
      perf_python_failed_first = python_failed_first,
      perf_python_max_formatted = python_max_formatted,
      perf_rust_max_formatted = rust_max_formatted
    ) %>%
    left_join(
      mem_analysis %>%
        select(
          method = task,
          mem_max_cells_python = max_cells_Python,
          mem_max_cells_rust = max_cells_Rust,
          mem_rust_advantage = rust_advantage_cells,
          mem_python_failed_first = python_failed_first,
          mem_efficiency_ratio = memory_efficiency_ratio,
          mem_at_max_python = memory_at_max_Python,
          mem_at_max_rust = memory_at_max_Rust
        ),
      by = "method"
    ) %>%
    left_join(
      mem_efficiency %>%
        filter(language == "Rust") %>%
        select(
          method = task,
          rust_avg_mem_per_cell = avg_memory_per_cell_kb,
          rust_max_mem = max_memory_gb,
          rust_min_mem = min_memory_gb
        ),
      by = "method"
    ) %>%
    left_join(
      mem_efficiency %>%
        filter(language == "Python") %>%
        select(
          method = task,
          python_avg_mem_per_cell = avg_memory_per_cell_kb,
          python_max_mem = max_memory_gb,
          python_min_mem = min_memory_gb
        ),
      by = "method"
    ) %>%
    mutate(
      # Calculate combined metrics
      consistent_failure_pattern = (!is.na(perf_python_failed_first) & !is.na(mem_python_failed_first)) &
                                   (perf_python_failed_first == mem_python_failed_first),
      perf_mem_correlation = case_when(
        is.na(perf_python_failed_first) | is.na(mem_python_failed_first) ~ "Insufficient data",
        perf_python_failed_first & mem_python_failed_first ~ "Both show Rust advantage",
        !perf_python_failed_first & !mem_python_failed_first ~ "Both show similar capacity",
        perf_python_failed_first & !mem_python_failed_first ~ "Performance advantage only",
        !perf_python_failed_first & mem_python_failed_first ~ "Memory advantage only",
        TRUE ~ "Mixed results"
      ),
      # Format memory efficiency ratio
      mem_efficiency_advantage = ifelse(!is.na(mem_efficiency_ratio),
                                       sprintf("%.1fx", mem_efficiency_ratio),
                                       "N/A"),
      # Memory per cell efficiency
      mem_per_cell_advantage = ifelse(!is.na(python_avg_mem_per_cell) & !is.na(rust_avg_mem_per_cell),
                                     python_avg_mem_per_cell / rust_avg_mem_per_cell,
                                     NA)
    ) %>%
    arrange(method)

  # Create summary statistics
  summary_stats <- list(
    total_methods = nrow(aligned_stats),
    performance_rust_advantages = sum(aligned_stats$perf_python_failed_first, na.rm = TRUE),
    memory_rust_advantages = sum(aligned_stats$mem_python_failed_first, na.rm = TRUE),
    consistent_patterns = sum(aligned_stats$consistent_failure_pattern, na.rm = TRUE),
    avg_perf_advantage = mean(aligned_stats$perf_rust_advantage[aligned_stats$perf_python_failed_first], na.rm = TRUE),
    avg_mem_advantage = mean(aligned_stats$mem_rust_advantage[aligned_stats$mem_python_failed_first], na.rm = TRUE),
    avg_mem_efficiency = mean(aligned_stats$mem_efficiency_ratio, na.rm = TRUE),
    avg_mem_per_cell_efficiency = mean(aligned_stats$mem_per_cell_advantage, na.rm = TRUE)
  )

  # Print detailed analysis
  cat("\n" , paste(rep("=", 80), collapse = ""), "\n")
  cat("ALIGNED MEMORY & PERFORMANCE STATISTICS\n")
  cat(paste(rep("=", 80), collapse = ""), "\n\n")

  cat("SUMMARY:\n")
  cat(sprintf("  Total methods analyzed: %d\n", summary_stats$total_methods))
  cat(sprintf("  Methods with Rust performance advantage: %d/%d\n",
              summary_stats$performance_rust_advantages, summary_stats$total_methods))
  cat(sprintf("  Methods with Rust memory advantage: %d/%d\n",
              summary_stats$memory_rust_advantages, summary_stats$total_methods))
  cat(sprintf("  Consistent failure patterns: %d/%d\n",
              summary_stats$consistent_patterns, summary_stats$total_methods))
  cat(sprintf("  Average performance advantage: %.1fx\n", summary_stats$avg_perf_advantage))
  cat(sprintf("  Average memory capacity advantage: %.1fx\n", summary_stats$avg_mem_advantage))
  cat(sprintf("  Average memory efficiency advantage: %.1fx\n", summary_stats$avg_mem_efficiency))
  cat(sprintf("  Average memory per cell efficiency: %.1fx\n", summary_stats$avg_mem_per_cell_efficiency))

  cat("\nDETAILED METHOD ANALYSIS:\n")
  cat(paste(rep("-", 80), collapse = ""), "\n")

  for (i in 1:nrow(aligned_stats)) {
    method_data <- aligned_stats[i, ]
    cat(sprintf("\n%s:\n", method_data$method))

    # Performance analysis
    cat("  Performance:\n")
    if (!is.na(method_data$perf_python_failed_first) && method_data$perf_python_failed_first) {
      cat(sprintf("    Python failed at: %s rows\n", method_data$perf_python_max_formatted))
      cat(sprintf("    Rust succeeded up to: %s rows\n", method_data$perf_rust_max_formatted))
      cat(sprintf("    Rust advantage: %.1fx more data\n", method_data$perf_rust_advantage))
    } else {
      cat("    Both achieved similar maximum dataset sizes\n")
    }

    # Memory analysis
    cat("  Memory:\n")
    if (!is.na(method_data$mem_python_failed_first) && method_data$mem_python_failed_first) {
      mem_py_formatted <- case_when(
        method_data$mem_max_cells_python >= 1e6 ~ paste0(round(method_data$mem_max_cells_python/1e6, 1), "M"),
        method_data$mem_max_cells_python >= 1e3 ~ paste0(round(method_data$mem_max_cells_python/1e3, 0), "K"),
        TRUE ~ as.character(method_data$mem_max_cells_python)
      )
      mem_rust_formatted <- case_when(
        method_data$mem_max_cells_rust >= 1e6 ~ paste0(round(method_data$mem_max_cells_rust/1e6, 1), "M"),
        method_data$mem_max_cells_rust >= 1e3 ~ paste0(round(method_data$mem_max_cells_rust/1e3, 0), "K"),
        TRUE ~ as.character(method_data$mem_max_cells_rust)
      )

      cat(sprintf("    Python failed at: %s cells (%.1f GB)\n",
                  mem_py_formatted, method_data$mem_at_max_python))
      cat(sprintf("    Rust succeeded up to: %s cells (%.1f GB)\n",
                  mem_rust_formatted, method_data$mem_at_max_rust))
      cat(sprintf("    Rust advantage: %.1fx more cells\n", method_data$mem_rust_advantage))
      cat(sprintf("    Memory efficiency: %s\n", method_data$mem_efficiency_advantage))
    } else {
      cat("    Both achieved similar maximum dataset sizes\n")
    }

    # Memory per cell efficiency
    if (!is.na(method_data$mem_per_cell_advantage)) {
      cat(sprintf("    Memory per cell: Rust is %.1fx more efficient\n", method_data$mem_per_cell_advantage))
    }

    # Correlation
    cat(sprintf("  Pattern: %s\n", method_data$perf_mem_correlation))
    cat(sprintf("  Consistent: %s\n", ifelse(method_data$consistent_failure_pattern, "Yes", "No")))
  }

  # Save detailed results
  if (detailed_output) {
    write_csv(aligned_stats, output_file)
    cat(sprintf("\nDetailed statistics saved to: %s\n", output_file))
  }

  # Generate and save plots
  plots_generated <- list()
  if (generate_plots) {
    cat("\nGenerating and saving plots...\n")

    # Create output directory if it doesn't exist
    if (!dir.exists(plot_output_dir)) {
      dir.create(plot_output_dir, recursive = TRUE)
    }

    # Save performance plot
    perf_plot_file <- file.path(plot_output_dir, "performance_comparison.pdf")
    ggsave(perf_plot_file, performance_result$plot,
           width = 4, height = 21, units = "in")
    cat(sprintf("  Performance plot saved to: %s\n", perf_plot_file))
    plots_generated$performance <- perf_plot_file

    # Save memory plot
    mem_plot_file <- file.path(plot_output_dir, "memory_comparison.pdf")
    ggsave(mem_plot_file, memory_result$plot,
           width = 4, height = 21, units = "in")
    cat(sprintf("  Memory plot saved to: %s\n", mem_plot_file))
    plots_generated$memory <- mem_plot_file

    # Create combined plot with 3x2 layout (3 rows, 2 columns)
    # Each row shows one pair: performance on left, memory on right
    library(patchwork)

    # Extract individual plots from both results
    # We need to create individual plots for each method to arrange them properly

    # Create a new combined layout: 3 rows x 2 columns
    # Row 1: DEG performance | DEG memory
    # Row 2: Filtering performance | Filtering memory
    # Row 3: HVG performance | HVG memory
    # Row 4: KNN performance | KNN memory
    # Row 5: Normalization performance | Normalization memory
    # Row 6: PCA performance | PCA memory

    # Since we have 6 methods, we'll create a 6x2 layout (not 3x2 as originally requested)
    # or we can create separate combined plots if you prefer 3x2

    # Option 1: 6x2 layout (all methods)
    combined_plot_6x2 <- performance_result$plot | memory_result$plot

    # Option 2: Create two separate 3x2 plots
    # For this, we need to generate separate plots for subsets of methods

    # Save the 6x2 version
    combined_plot_file <- file.path(plot_output_dir, "combined_performance_memory.pdf")
    ggsave(combined_plot_file, combined_plot_6x2,
           width = 8, height = 21, units = "in")
    cat(sprintf("  Combined plot (6x2) saved to: %s\n", combined_plot_file))
    plots_generated$combined <- combined_plot_file

    # Create a more compact 3x2 layout by grouping methods
    # First 3 methods
    tryCatch({
      # Generate plots for first 3 methods
      perf_result_1 <- compare_performance_failure_points(
        rust_file_paths = rust_performance_files[1:3],
        python_file_paths = python_performance_files[1:3],
        method_names = method_names[1:3],
        plot_layout = list(ncol = 1, nrow = 3),
        show_speedup_annotations = FALSE,
        show_statistics_box = TRUE,
        plot_width = 4,
        plot_height = 9,  # 3 methods * 3 inches each
        y_axis_on_all_plots = TRUE,
        show_global_legend = FALSE,
        base_text_size = 10,
        x_scale_power = 0.65,
        rust_color = "#CE422B",
        python_color = "#3776AB"
      )

      mem_result_1 <- compare_memory_consumption(
        rust_file_paths = rust_memory_files[1:3],
        python_file_paths = python_memory_files[1:3],
        method_names = method_names[1:3],
        plot_layout = list(ncol = 1, nrow = 3),
        show_speedup_annotations = FALSE,
        show_statistics_box = TRUE,
        plot_width = 4,
        plot_height = 9,  # 3 methods * 3 inches each
        y_axis_on_all_plots = TRUE,
        show_global_legend = FALSE,
        base_text_size = 10,
        x_scale_power = 0.65,
        rust_color = "#CE422B",
        python_color = "#3776AB"
      )

      combined_3x2_part1 <- perf_result_1$plot | mem_result_1$plot
      combined_3x2_file1 <- file.path(plot_output_dir, "combined_3x2_part1.pdf")
      ggsave(combined_3x2_file1, combined_3x2_part1,
             width = 8, height = 9, units = "in")
      cat(sprintf("  Combined plot (3x2 part 1) saved to: %s\n", combined_3x2_file1))
      plots_generated$combined_3x2_part1 <- combined_3x2_file1

      # Second 3 methods
      perf_result_2 <- compare_performance_failure_points(
        rust_file_paths = rust_performance_files[4:6],
        python_file_paths = python_performance_files[4:6],
        method_names = method_names[4:6],
        plot_layout = list(ncol = 1, nrow = 3),
        show_speedup_annotations = FALSE,
        show_statistics_box = TRUE,
        plot_width = 4,
        plot_height = 9,  # 3 methods * 3 inches each
        y_axis_on_all_plots = TRUE,
        show_global_legend = FALSE,
        base_text_size = 10,
        x_scale_power = 0.65,
        rust_color = "#CE422B",
        python_color = "#3776AB"
      )

      mem_result_2 <- compare_memory_consumption(
        rust_file_paths = rust_memory_files[4:6],
        python_file_paths = python_memory_files[4:6],
        method_names = method_names[4:6],
        plot_layout = list(ncol = 1, nrow = 3),
        show_speedup_annotations = FALSE,
        show_statistics_box = TRUE,
        plot_width = 4,
        plot_height = 9,  # 3 methods * 3 inches each
        y_axis_on_all_plots = TRUE,
        show_global_legend = FALSE,
        base_text_size = 10,
        x_scale_power = 0.65,
        rust_color = "#CE422B",
        python_color = "#3776AB"
      )

      combined_3x2_part2 <- perf_result_2$plot | mem_result_2$plot
      combined_3x2_file2 <- file.path(plot_output_dir, "combined_3x2_part2.pdf")
      ggsave(combined_3x2_file2, combined_3x2_part2,
             width = 8, height = 9, units = "in")
      cat(sprintf("  Combined plot (3x2 part 2) saved to: %s\n", combined_3x2_file2))
      plots_generated$combined_3x2_part2 <- combined_3x2_file2

    }, error = function(e) {
      cat("  Note: Could not create 3x2 layouts, using standard layout only\n")
    })

    # Also save as PNG for easier viewing
    perf_png_file <- file.path(plot_output_dir, "performance_comparison.png")
    ggsave(perf_png_file, performance_result$plot,
           width = 4, height = 21, units = "in", dpi = 300)
    plots_generated$performance_png <- perf_png_file

    mem_png_file <- file.path(plot_output_dir, "memory_comparison.png")
    ggsave(mem_png_file, memory_result$plot,
           width = 4, height = 21, units = "in", dpi = 300)
    plots_generated$memory_png <- mem_png_file

    combined_png_file <- file.path(plot_output_dir, "combined_performance_memory.png")
    ggsave(combined_png_file, combined_plot_6x2,
           width = 8, height = 21, units = "in", dpi = 300)
    plots_generated$combined_png <- combined_png_file

    # Also save 3x2 PNG versions if they were created
    if ("combined_3x2_part1" %in% names(plots_generated)) {
      combined_3x2_png1 <- file.path(plot_output_dir, "combined_3x2_part1.png")
      ggsave(combined_3x2_png1, combined_3x2_part1,
             width = 8, height = 9, units = "in", dpi = 300)
      plots_generated$combined_3x2_part1_png <- combined_3x2_png1

      combined_3x2_png2 <- file.path(plot_output_dir, "combined_3x2_part2.png")
      ggsave(combined_3x2_png2, combined_3x2_part2,
             width = 8, height = 9, units = "in", dpi = 300)
      plots_generated$combined_3x2_part2_png <- combined_3x2_png2
    }

    cat(sprintf("  PNG versions also saved in: %s\n", plot_output_dir))
  }

  # Return results
  return(list(
    aligned_statistics = aligned_stats,
    summary = summary_stats,
    performance_result = performance_result,
    memory_result = memory_result,
    plots_generated = plots_generated
  ))
}

# Usage example
cat("Loading file paths...\n")

# Define file paths (matching your existing setup)
rust_performance_files <- c(
  "/local/06-24_single-rust-project/single-performance/timings/rust/deg_performance.csv",
  "/local/06-24_single-rust-project/single-performance/timings/rust/filtering_performance.csv",
  "/local/06-24_single-rust-project/single-performance/timings/rust/hvg_performance.csv",
  "/local/06-24_single-rust-project/single-performance/timings/rust/knn_performance.csv",
  "/local/06-24_single-rust-project/single-performance/timings/rust/normalization_performance.csv",
  "/local/06-24_single-rust-project/single-performance/timings/rust/pca_performance.csv"
)

python_performance_files <- c(
  "/local/06-24_single-rust-project/single-performance/timings/python/deg_performance.csv",
  "/local/06-24_single-rust-project/single-performance/timings/python/filtering_performance.csv",
  "/local/06-24_single-rust-project/single-performance/timings/python/hvg_performance.csv",
  "/local/06-24_single-rust-project/single-performance/timings/python/knn_performance.csv",
  "/local/06-24_single-rust-project/single-performance/timings/python/normalization_performance.csv",
  "/local/06-24_single-rust-project/single-performance/timings/python/pca_performance.csv"
)

rust_memory_files <- c(
  "/local/06-24_single-rust-project/single-performance/timings/rust/deg_memory.csv",
  "/local/06-24_single-rust-project/single-performance/timings/rust/filtering_memory.csv",
  "/local/06-24_single-rust-project/single-performance/timings/rust/hvg_memory.csv",
  "/local/06-24_single-rust-project/single-performance/timings/rust/knn_memory.csv",
  "/local/06-24_single-rust-project/single-performance/timings/rust/normalization_memory.csv",
  "/local/06-24_single-rust-project/single-performance/timings/rust/pca_memory.csv"
)

python_memory_files <- c(
  "/local/06-24_single-rust-project/single-performance/timings/python/deg_memory.csv",
  "/local/06-24_single-rust-project/single-performance/timings/python/filtering_memory.csv",
  "/local/06-24_single-rust-project/single-performance/timings/python/hvg_memory.csv",
  "/local/06-24_single-rust-project/single-performance/timings/python/knn_memory.csv",
  "/local/06-24_single-rust-project/single-performance/timings/python/normalization_memory.csv",
  "/local/06-24_single-rust-project/single-performance/timings/python/pca_memory.csv"
)

# Run the alignment analysis
cat("Running alignment analysis...\n")
aligned_results <- align_memory_performance_statistics(
  rust_performance_files = rust_performance_files,
  python_performance_files = python_performance_files,
  rust_memory_files = rust_memory_files,
  python_memory_files = python_memory_files,
  method_names = c("DEG", "Filtering", "HVG", "KNN", "Normalization", "PCA"),
  output_file = "aligned_memory_performance_statistics.csv",
  detailed_output = TRUE,
  generate_plots = TRUE,
  plot_output_dir = "aligned_plots"
)

cat("\nAlignment analysis complete!\n")
