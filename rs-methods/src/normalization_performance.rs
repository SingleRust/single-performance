use anyhow::Result;
use clap::Parser;
use rayon::ThreadPoolBuilder;
use single_rust::io;
use single_rust::shared::Precision;
use single_utilities::types::Direction;
use std::fs::{File, OpenOptions};
use std::io::Write;
use std::path::PathBuf;
use std::time::Instant;
use single_rust::memory::processing::{log1p_expression, normalize_expression};

#[derive(Parser)]
#[command(
    author = "Ian Ferenc Diks",
    version,
    about = "Benchmark for H5AD normalization and log1p transformation"
)]
struct Args {
    /// Path to H5AD file
    #[arg(short = 'i', long = "input")]
    h5ad_path: PathBuf,
    /// Path for output CSV file
    #[arg(short = 'o', long = "output")]
    csv_path: PathBuf,
    /// Number of benchmark iterations to run
    #[arg(short = 'n', long = "iterations", default_value = "5")]
    iterations: usize,
    /// Normalization target (counts per cell)
    #[arg(short = 't', long = "target", default_value = "10000")]
    normalization_target: u32,
}

fn log_benchmark_result(
    file: &mut File,
    test_case: &str,
    step: &str,
    duration_ms: u128,
    success: bool,
    dimensions: Option<(usize, usize)>,
) -> Result<()> {
    let dims_str = if let Some((obs, vars)) = dimensions {
        format!("{},{}", obs, vars)
    } else {
        String::from(",")
    };

    writeln!(
        file,
        "{},{},{},{},{}",
        test_case, step, duration_ms, success, dims_str
    )?;
    Ok(())
}

fn run_normalization_benchmark(
    h5ad_path: &PathBuf,
    csv_path: &PathBuf,
    iterations: usize,
    normalization_target: u32,
) -> Result<()> {
    println!("Starting normalization + log1p benchmark");

    let thread_pool = ThreadPoolBuilder::new().num_threads(128).build()?;

    // Create or open the output CSV file
    let mut file = OpenOptions::new()
        .create(true)
        .append(true)
        .open(csv_path)?;

    // Add header if file is new
    if file.metadata()?.len() == 0 {
        writeln!(
            file,
            "test_case,test_step,time_ms,success,rows,cols"
        )?;
    }

    // First, do a quick load to get dimensions
    println!("Getting dataset dimensions from {}", h5ad_path.display());
    let dimensions = match io::read_h5ad_fast_memory(h5ad_path) {
        Ok(adata) => {
            let dims = (adata.n_obs(), adata.n_vars());
            println!("Dataset dimensions: {} cells x {} genes", dims.0, dims.1);
            drop(adata); // not necessary here
            Some(dims)
        }
        Err(e) => {
            println!("Failed to read dataset for dimensions: {}", e);
            None
        }
    };

    println!("Running normalization + log1p benchmarks for {} iterations", iterations);

    // Run benchmarks with fresh load each iteration
    for i in 1..=iterations {
        println!("  Iteration {}/{}", i, iterations);

        let load_start = Instant::now();
        let adata_result = io::read_h5ad_fast_memory(h5ad_path);
        let load_time = load_start.elapsed().as_millis();
        let load_success = adata_result.is_ok();

        log_benchmark_result(
            &mut file,
            &format!("normalization_log1p_iter{}", i),
            "load",
            load_time,
            load_success,
            dimensions,
        )?;

        if let Ok(adata) = adata_result {
            let x = adata.x();

            let combined_start = Instant::now();

            let norm_result = thread_pool.install(|| {
                normalize_expression(
                &x,
                normalization_target,
                &Direction::ROW,
                Some(Precision::Single)
            )
            });

            let log1p_result = if norm_result.is_ok() {
                thread_pool.install(|| {log1p_expression(&x, Some(Precision::Single))})
            } else {
                Err(anyhow::anyhow!("Normalization failed"))
            };

            let combined_time = combined_start.elapsed().as_millis();
            let combined_success = log1p_result.is_ok();

            log_benchmark_result(
                &mut file,
                &format!("normalization_log1p_iter{}", i),
                "process",
                combined_time,
                combined_success,
                dimensions,
            )?;

            drop(adata);
        } else if let Err(e) = adata_result {
            println!("    Failed to load dataset in iteration {}: {}", i, e);
        }
    }

    println!("Benchmark completed successfully");
    Ok(())
}

fn main() -> Result<()> {
    let args = Args::parse();
    println!("Starting normalization + log1p benchmark with:");
    println!("  Input file: {}", args.h5ad_path.display());
    println!("  Output file: {}", args.csv_path.display());
    println!("  Iterations: {}", args.iterations);
    println!("  Normalization target: {}", args.normalization_target);

    run_normalization_benchmark(
        &args.h5ad_path,
        &args.csv_path,
        args.iterations,
        args.normalization_target
    )?;

    println!("Benchmark completed");
    Ok(())
}