use anyhow::Result;
use clap::Parser;
use polars::io::SerReader;
use polars::prelude::CsvReadOptions;
use rayon::ThreadPoolBuilder;
use single_algebra::dimred::pca::{PowerIterationNormalizer, SVDMethod};
use single_rust::io;
use single_rust::memory::processing::compute_highly_variable_genes;
use single_rust::memory::processing::dimred::pca::run_pca_sparse_masked;
use single_rust::memory::processing::dimred::FeatureSelectionMethod;
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
        "pca_truncated,{},{},{},{}",
        step, duration_ms, success, dims_str
    )?;
    Ok(())
}

fn run_pca_benchmark(h5ad_path: &PathBuf, csv_path: &PathBuf, iterations: usize) -> Result<()> {
    println!("Starting HVG benchmark");

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
        compute_highly_variable_genes(&adata, None)?;
        single_rust::memory::processing::log1p_expression(&x, Some(Precision::Single))?;

        let df = CsvReadOptions::default()
            .with_has_header(true)
            .try_into_reader_with_file_path(Some(
                "/local/06-24_single-rust-project/single_test/var_data.csv".into(),
            ))?
            .finish()?;
        let var_df_col = df.column("highly_variable")?;
        let bool_col = var_df_col.bool()?;
        let s: Vec<bool> = bool_col
            .into_iter()
            .map(|opt_val| opt_val.unwrap_or(false))
            .collect();
        println!(
            "Dataset loaded: {} cells x {} genes",
            adata.n_obs(),
            adata.n_vars()
        );
        println!("Running PCA computation for {} iterations", iterations);

        for i in 1..=iterations {
            println!("  Iteration {}/{}", i, iterations);

            let x_clone = adata.x();
            let pca_start = Instant::now();
            let pca_result = thread_pool.install(|| {
                run_pca_sparse_masked::<f32>(
                    &x_clone,
                    Some(FeatureSelectionMethod::HighlyVariableSelection(s.clone())),
                    Some(false),
                    Some(true),
                    Some(50),
                    None,
                    Some(42),
                    //Some(SVDMethod::Lanczos)
                    Some(SVDMethod::Random {
                        n_oversamples: 10,
                        n_power_iterations: 7,
                        normalizer: PowerIterationNormalizer::LU,
                    }),
                )
            });

            let pca_time = pca_start.elapsed().as_millis();
            let pca_success = pca_result.is_ok();

            log_benchmark_result(
                &mut file,
                "compute_pca_truncated",
                pca_time,
                pca_success,
                Some((adata.n_obs(), adata.n_vars())),
            )?;
        }

        println!("PCA benchmark completed successfully");
    } else if let Err(e) = adata_result {
        log_benchmark_result(&mut file, "load", load_time, load_success, None)?;
        println!("Failed to load dataset: {}", e);
    }

    Ok(())
}

fn main() -> Result<()> {
    let args = Args::parse();

    println!("Starting PCA benchmark with:");
    println!("  Input file: {}", args.h5ad_path.display());
    println!("  Output file: {}", args.csv_path.display());
    println!("  Iterations: {}", args.iterations);

    run_pca_benchmark(&args.h5ad_path, &args.csv_path, args.iterations)?;

    println!("Benchmark completed");
    Ok(())
}
