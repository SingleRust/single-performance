use anyhow::Result;
use clap::Parser;
use rayon::ThreadPoolBuilder;
use single_rust::io;
use single_rust::memory::processing::diffexp::{CorrectionMethod, rank_gene_groups};
use single_rust::shared::Precision;
use single_statistics::testing::{TTestType, TestMethod};
use single_utilities::types::Direction;
use std::fs::{File, OpenOptions};
use std::io::Write;
use std::path::PathBuf;
use std::time::Instant;

#[derive(Parser)]
#[command(
    author = "Ian Ferenc Diks",
    version,
    about = "Benchmark for DEG (Differential Expression Gene) operations"
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
        "deg,{},{},{},{}",
        step, duration_ms, success, dims_str
    )?;
    Ok(())
}

fn run_deg_benchmark(h5ad_path: &PathBuf, csv_path: &PathBuf, iterations: usize) -> Result<()> {
    println!("Starting DEG benchmark");

    // Create or open the output CSV file
    let mut file = OpenOptions::new()
        .create(true)
        .append(true)
        .open(csv_path)?;

    // Add header if file is new
    if file.metadata()?.len() == 0 {
        writeln!(file, "test_case,test_step,time_ms,success,rows,cols")?;
    }

    let thread_pool = ThreadPoolBuilder::new().num_threads(128).build()?;

    // Load the dataset
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

        println!("Preprocessing data");

        let x = adata.x();
        single_rust::memory::processing::normalize_expression(
            &x,
            10_000,
            &Direction::ROW,
            Some(Precision::Single),
        )?;
        single_rust::memory::processing::log1p_expression(&x, Some(Precision::Single))?;
        println!(
            "Dataset loaded and preprocessed: {} cells x {} genes",
            adata.n_obs(),
            adata.n_vars()
        );
        println!("Running DEG computation for {} iterations", iterations);
        let key = "neural_vs_meso".to_string();

        for i in 1..=iterations {
            println!("  Iteration {}/{}", i, iterations);

            let deg_start = Instant::now();
            let deg_result = thread_pool.install(|| {
                rank_gene_groups(
                    &adata,
                    "cell_name",
                    Some("SW480"),
                    Some(&["NCI-H460"]),
                    Some(&key),
                    Some(TestMethod::TTest(TTestType::Welch)),
                    None,
                    CorrectionMethod::BejaminiHochberg,
                    Some(true),
                    Some(1e-9),
                )
            });

            let deg_time = deg_start.elapsed().as_millis();
            let deg_success = deg_result.is_ok();

            log_benchmark_result(
                &mut file,
                "compute_deg_for_group",
                deg_time,
                deg_success,
                Some((adata.n_obs(), adata.n_vars())),
            )?;
        }

        println!("DEG benchmark completed successfully");
    } else if let Err(e) = adata_result {
        log_benchmark_result(&mut file, "load", load_time, load_success, None)?;
        println!("Failed to load dataset: {}", e);
    }

    Ok(())
}

fn main() -> Result<()> {
    let args = Args::parse();

    println!("Starting DEG benchmark with:");
    println!("  Input file: {}", args.h5ad_path.display());
    println!("  Output file: {}", args.csv_path.display());
    println!("  Iterations: {}", args.iterations);

    run_deg_benchmark(&args.h5ad_path, &args.csv_path, args.iterations)?;

    println!("Benchmark completed");
    Ok(())
}
