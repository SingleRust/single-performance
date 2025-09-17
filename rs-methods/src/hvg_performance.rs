use anndata_memory::DeepClone;
use anyhow::Result;
use clap::Parser;
use rayon::ThreadPoolBuilder;
use single_rust::io;
use single_rust::memory::processing::compute_highly_variable_genes;
use single_rust::shared::Precision;
use single_utilities::types::Direction;
use std::fs::{File, OpenOptions};
use std::io::Write;
use std::path::PathBuf;
use std::time::Instant;

#[derive(Parser)]
#[command(
    author = "Ian Ferenc Diks",
    version,
    about = "Benchmark for H5AD operations"
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
}

fn log_benchmark_result(
    file: &mut File,
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
        "hvg,{},{},{},{}",
        step, duration_ms, success, dims_str
    )?;
    Ok(())
}

fn run_hvg_benchmark(h5ad_path: &PathBuf, csv_path: &PathBuf, iterations: usize) -> Result<()> {
    println!("Starting HVG benchmark");
    let thread_pool = ThreadPoolBuilder::new().num_threads(128).build()?;

    // Create or open the output CSV file
    let mut file = OpenOptions::new()
        .create(true)
        .append(true)
        .open(csv_path)?;

    // Add header if file is new
    if file.metadata()?.len() == 0 {
        writeln!(file, "test_case,test_step,time_ms,success,rows,cols")?;
    }

    println!("Loading dataset from {}", h5ad_path.display());
    let load_start = Instant::now();
    let adata_result = io::read_h5ad_fast_memory(h5ad_path);
    let load_time = load_start.elapsed().as_millis();
    let load_success = adata_result.is_ok();

    if let Ok(adata) = adata_result {
        log_benchmark_result(
            &mut file,
            "load",
            load_time,
            load_success,
            Some((adata.n_obs(), adata.n_vars())),
        )?;

        let x = adata.x();
        single_rust::memory::processing::normalize_expression(
            &x,
            1e4 as u32,
            &Direction::ROW,
            Some(Precision::Single),
        )?;

        println!(
            "Dataset loaded: {} cells x {} genes",
            adata.n_obs(),
            adata.n_vars()
        );
        println!("Running HVG computation for {} iterations", iterations);

        let var_df = adata.var().deep_clone();

        for i in 1..=iterations {
            println!("  Iteration {}/{}", i, iterations);

            let hvg_start = Instant::now();
            let hvg_result = thread_pool.install(|| compute_highly_variable_genes(&adata, None));
            let hvg_time = hvg_start.elapsed().as_millis();
            let hvg_success = hvg_result.is_ok();

            log_benchmark_result(
                &mut file,
                "compute_hvg",
                hvg_time,
                hvg_success,
                Some((adata.n_obs(), adata.n_vars())),
            )?;
            adata.var().set_data(var_df.get_data())?;
        }

        println!("HVG benchmark completed successfully");
    } else if let Err(e) = adata_result {
        // Log the failed loading step
        log_benchmark_result(&mut file, "load", load_time, load_success, None)?;
        println!("Failed to load dataset: {}", e);
    }

    Ok(())
}

fn main() -> Result<()> {
    let args = Args::parse();

    println!("Starting HVG benchmark with:");
    println!("  Input file: {}", args.h5ad_path.display());
    println!("  Output file: {}", args.csv_path.display());
    println!("  Iterations: {}", args.iterations);

    run_hvg_benchmark(&args.h5ad_path, &args.csv_path, args.iterations)?;

    println!("Benchmark completed");
    Ok(())
}
