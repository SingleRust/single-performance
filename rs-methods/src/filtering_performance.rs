use anndata::AnnDataOp;
use anyhow::Result;
use clap::Parser;
use ndarray::Array1;
use single_rust::io;
use single_rust::memory::processing::filtering::{mark_filter_cells, mark_filter_genes};
use std::fs::{File, OpenOptions};
use std::io::Write;
use std::path::PathBuf;
use std::time::Instant;

#[derive(Parser)]
#[command(
    author = "Ian Ferenc Diks",
    version,
    about = "Benchmark for H5AD filtering operations"
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
    test_case: &str,
    step: &str,
    duration_ms: u128,
    success: bool,
    dimensions: Option<(usize, usize)>,
    filtered_count: Option<usize>,
) -> Result<()> {
    let dims_str = if let Some((obs, vars)) = dimensions {
        format!("{},{}", obs, vars)
    } else {
        String::from(",")
    };

    let filtered_str = if let Some(count) = filtered_count {
        count.to_string()
    } else {
        String::from("")
    };

    writeln!(
        file,
        "{},{},{},{},{},{}",
        test_case, step, duration_ms, success, dims_str, filtered_str
    )?;
    Ok(())
}

fn run_filtering_benchmark(
    h5ad_path: &PathBuf,
    csv_path: &PathBuf,
    iterations: usize,
) -> Result<()> {
    println!("Starting filtering benchmark");

    // Create or open the output CSV file
    let mut file = OpenOptions::new()
        .create(true)
        .append(true)
        .open(csv_path)?;

    // Add header if file is new
    if file.metadata()?.len() == 0 {
        writeln!(
            file,
            "test_case,test_step,time_ms,success,rows,cols,filtered_count"
        )?;
    }

    // Initial load to get dimensions
    println!("Loading dataset once to get dimensions");
    let temp_load_start = Instant::now();
    let temp_adata = io::read_h5ad(h5ad_path, io::FileScope::Read, false);
    let temp_load_time = temp_load_start.elapsed().as_millis();

    let (n_obs, n_vars) = if let Ok(adata) = &temp_adata {
        (adata.n_obs(), adata.n_vars())
    } else {
        log_benchmark_result(
            &mut file,
            "initial_dimensions_check",
            "load",
            temp_load_time,
            false,
            None,
            None,
        )?;
        println!("Failed to load dataset for initial dimension check");
        return Err(anyhow::anyhow!(
            "Failed to load dataset for initial dimension check"
        ));
    };

    println!("Dataset dimensions: {} cells x {} genes", n_obs, n_vars);

    // Drop the temporary dataset
    drop(temp_adata);

    // Combined filtering benchmark
    for i in 1..=iterations {
        println!("  Combined filtering - Iteration {}/{}", i, iterations);

        // Load fresh data for this iteration
        let load_start = Instant::now();
        let adata_result = io::read_h5ad_fast_memory(h5ad_path);
        let load_time = load_start.elapsed().as_millis();

        if let Ok(mut adata) = adata_result {
            log_benchmark_result(
                &mut file,
                "combined_filtering",
                "load",
                load_time,
                true,
                Some((n_obs, n_vars)),
                None,
            )?;

            let combined_start = Instant::now();

            // First filter cells
            let cell_filter_result = mark_filter_cells::<u32, f64>(
                &adata,
                Some(200),    // min_genes
                Some(5000),   // max_genes
                Some(1000.0), // min_counts
                None,         // max_counts
                None,         // min_fraction
                None,         // max_fraction
            );

            // Then filter genes
            let gene_filter_result = if cell_filter_result.is_ok() {
                mark_filter_genes::<u32, f64>(
                    &adata,
                    Some(3),    // min_cells
                    None,       // max_cells
                    Some(10.0), // min_counts
                    None,       // max_counts
                    None,       // min_fraction
                    Some(0.9),  // max_fraction
                )
            } else {
                Err(anyhow::anyhow!("Cell filtering failed"))
            };

            // Calculate filtered counts if both operations succeeded
            if let (Ok(cell_filter), Ok(gene_filter)) = (cell_filter_result, gene_filter_result) {
                // Convert to Array1 and perform the subsetting
                let cell_filter = Array1::from(cell_filter);
                let gene_filter = Array1::from(gene_filter);

                let obs_sel_info =
                    single_rust::shared::get_select_info_obs(Some(cell_filter.view()))?;
                let var_sel_info =
                    single_rust::shared::get_select_info_vars(Some(gene_filter.view()))?;

                // Apply both filters
                adata.subset_inplace(&[&obs_sel_info[0], &var_sel_info[1]])?;

                let combined_time = combined_start.elapsed().as_millis();
                let filtered_total = (n_obs - adata.n_obs()) + (n_vars - adata.n_vars());

                log_benchmark_result(
                    &mut file,
                    "combined_filtering",
                    "filter",
                    combined_time,
                    true,
                    Some((n_obs, n_vars)),
                    Some(filtered_total),
                )?;
            } else {
                let combined_time = combined_start.elapsed().as_millis();
                log_benchmark_result(
                    &mut file,
                    "combined_filtering",
                    "filter",
                    combined_time,
                    false,
                    Some((n_obs, n_vars)),
                    None,
                )?;
            }
        } else {
            log_benchmark_result(
                &mut file,
                "combined_filtering",
                "load",
                load_time,
                false,
                Some((n_obs, n_vars)),
                None,
            )?;
            println!(
                "Failed to load dataset for combined filtering iteration {}",
                i
            );
        }
    }

    println!("Filtering benchmarks completed successfully");

    Ok(())
}

fn main() -> Result<()> {
    let args = Args::parse();
    println!("Starting filtering benchmark with:");
    println!("  Input file: {}", args.h5ad_path.display());
    println!("  Output file: {}", args.csv_path.display());
    println!("  Iterations: {}", args.iterations);

    run_filtering_benchmark(&args.h5ad_path, &args.csv_path, args.iterations)?;

    println!("Benchmark completed");
    Ok(())
}
