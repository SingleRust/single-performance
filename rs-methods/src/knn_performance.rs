use anyhow::Result;
use clap::Parser;
use kiddo::SquaredEuclidean;
use polars::io::SerReader;
use polars::prelude::CsvReadOptions;
use rayon::ThreadPoolBuilder;
use single_algebra::dimred::pca::{PowerIterationNormalizer, SVDMethod};
use single_clustering::neighborhood::knn_arrayd_adaptive;
use single_rust::io;
use single_rust::memory::processing::dimred::FeatureSelectionMethod;
use single_rust::shared::Precision;
use single_utilities::types::Direction;
use std::fs::{File, OpenOptions};
use std::io::Write;
use std::path::PathBuf;
use std::time::Instant;

#[derive(Parser)]
#[command(
    author = "Your Name",
    version,
    about = "Benchmark for nearest neighbor computation using KNN"
)]
struct Args {
    /// Path to NPY file containing the data
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
    test_case: &str,
    step: &str,
    duration_ms: u128,
    success: bool,
    dimensions: Option<(usize, usize)>,
) -> Result<()> {
    let dims_str = if let Some((rows, cols)) = dimensions {
        format!("{},{}", rows, cols)
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

fn run_knn_benchmark(h5ad_path: &PathBuf, csv_path: &PathBuf, iterations: usize) -> Result<()> {
    let thread_pool = ThreadPoolBuilder::new().num_threads(128).build()?;
    println!("Starting KNN benchmark");

    let mut file = OpenOptions::new()
        .create(true)
        .append(true)
        .open(csv_path)?;

    if file.metadata()?.len() == 0 {
        writeln!(file, "test_case,test_step,time_ms,success,rows,cols")?;
    }


    println!("Loading dataset from {}", h5ad_path.display());
    let load_start = Instant::now();
    let adata = io::read_h5ad_fast_memory(h5ad_path)?;
    let x = adata.x();
    single_rust::memory::processing::normalize_expression(
        &x,
        1e4 as u32,
        &Direction::ROW,
        Some(Precision::Single),
    )?;
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
    let pca_res = thread_pool.install(|| {single_rust::memory::processing::dimred::pca::run_pca_sparse_masked::<f32>(
        &x,
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
    )})?;
    let data: ndarray::ArrayBase<ndarray::OwnedRepr<f32>, ndarray::Dim<ndarray::IxDynImpl>> =
        pca_res.transformed.into_dyn();

    let load_time = load_start.elapsed().as_millis();
    let load_success = true;

    let shape = data.shape();
    let n_samples = shape[0];
    let n_features = shape[1];

    log_benchmark_result(
        &mut file,
        "knn",
        "load",
        load_time,
        load_success,
        Some((n_samples, n_features)),
    )?;

    println!(
        "Dataset loaded: {} samples x {} features",
        n_samples, n_features
    );
    println!("Running KNN benchmarks for {} iterations", iterations);

    // KNN computation benchmark
    for i in 1..=iterations {
        println!("  Iteration {}/{}", i, iterations);

        let knn_start = Instant::now();
        let knn_result = knn_arrayd_adaptive::<_, 30, SquaredEuclidean>(data.view(), 15);
        let knn_time = knn_start.elapsed().as_millis();
        let knn_success = knn_result.is_ok();

        log_benchmark_result(
            &mut file,
            "knn",
            "compute",
            knn_time,
            knn_success,
            Some((n_samples, n_features)),
        )?;

        if let Err(e) = knn_result {
            println!("    KNN computation failed: {}", e);
        }
    }

    println!("KNN benchmark completed successfully");

    Ok(())
}

fn main() -> Result<()> {
    let args = Args::parse();

    println!("Starting KNN benchmark with:");
    println!("  Input file: {}", args.h5ad_path.display());
    println!("  Output file: {}", args.csv_path.display());
    println!("  Iterations: {}", args.iterations);
    let thread_pool = ThreadPoolBuilder::new().num_threads(128).build()?;
    thread_pool.install(|| run_knn_benchmark(&args.h5ad_path, &args.csv_path, args.iterations))?;

    println!("Benchmark completed");
    Ok(())
}
