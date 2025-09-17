# Single Performance Benchmarking Repository

This repository contains the complete benchmarking suite used to evaluate and compare the performance of Python and Rust implementations for various computational biology and data analysis tasks.

## Purpose

We have made this repository publicly available to provide full transparency into our benchmarking methodology and results. This allows users to:

- **Understand our methodology**: See exactly how we conducted our performance comparisons
- **Reproduce our results**: Run the same benchmarks on your own systems
- **Verify our findings**: Independently validate our performance claims
- **Extend our work**: Use our framework to benchmark additional methods or implementations

## Repository Structure

### Implementation Directories
- **`py-methods/`** - Python implementations of various algorithms
- **`rs-methods/`** - Rust implementations of the same algorithms

### Benchmarking Components
- **`runners/`** - Shell scripts to execute benchmarks for both Python and Rust
- **`timings/`** - Raw performance and memory usage data from benchmark runs
- **`analysis/`** - R scripts for statistical analysis and visualization
- **`analysis-results/`** - Generated plots and statistical summaries

### Benchmarked Algorithms

The repository includes performance comparisons for:

- **DEG (Differential Expression Gene) Analysis** - Gene expression analysis
- **Filtering** - Data filtering operations
- **HVG (Highly Variable Genes)** - Gene variability analysis
- **KNN (K-Nearest Neighbors)** - Nearest neighbor computations
- **Normalization** - Data normalization methods
- **PCA (Principal Component Analysis)** - Dimensionality reduction

## Running Benchmarks

### Prerequisites
- Python 3.x with required scientific computing libraries (Scanpy 1.11.1)
  - The complete Python environment specification can be found in `environment-py.yml`
  - To recreate the environment: `conda env create -f environment-py.yml`
- Rust toolchain (latest stable)
- R for analysis scripts

### Test Data
Due to file size restrictions, we are currently unable to publish the count files and datasets used in our testing within this repository. However, if you are interested in obtaining these data files to reproduce our exact benchmarks, please feel free to contact us directly and we will be happy to share them with you.

### Execution
1. Navigate to the appropriate runner directory (`runners/python/` or `runners/rust/`)
2. Execute the desired benchmark script (e.g., `./pca_performance.sh`)
3. Results will be saved to the `timings/` directory

### Analysis
Run the R scripts in the `analysis/` directory to generate statistical summaries and visualizations from the timing data.

## Transparency and Reproducibility

This repository represents our commitment to open science and reproducible research. All code, data, and analysis methods are provided as-is, allowing the community to:

- Scrutinize our methods
- Identify potential improvements
- Build upon our work
- Contribute back to the community

## Contributing

We welcome contributions, suggestions, and feedback. Feel free to open issues or submit pull requests to help improve our benchmarking methodology.

---

*This benchmarking suite was created to provide objective, reproducible performance comparisons between Python and Rust implementations of common computational biology algorithms.*