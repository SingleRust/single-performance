import os
import argparse
import time
import csv
import scanpy as sc
import numpy as np
from multiprocessing import Pool, cpu_count
import pandas as pd

def log_benchmark_result(file, step, duration_ms, success, dimensions=None):
    """Log benchmark results to a CSV file."""
    dims_str = f"{dimensions[0]},{dimensions[1]}" if dimensions else ","
    file.write(f"deg,{step},{duration_ms},{success},{dims_str}\n")
    file.flush()

def run_deg_analysis(adata, iteration):
    """Run differential expression analysis for one iteration."""
    start_time = time.time()
    success = True
    
    try:
        # Equivalent parameters to the Rust implementation
        sc.tl.rank_genes_groups(
            adata,
            groupby='cell_name',
            reference='SW480',
            groups=['NCI-H460'],
            key_added=f'neural_vs_meso',
            method='t-test',  # Equivalent to TestMethod::TTest(TTestType::Welch)
            pts=True,               # Calculate percentage of cells expressing each gene
            corr_method='benjamini-hochberg',  # Equivalent to CorrectionMethod::BejaminiHochberg
        )
    except Exception as e:
        print(f"Error in iteration {iteration}: {e}")
        success = False
    
    duration_ms = int((time.time() - start_time) * 1000)
    return {
        'step': 'compute_deg_for_group',
        'duration_ms': duration_ms,
        'success': success,
        'iteration': iteration
    }

def run_deg_benchmark(h5ad_path, csv_path, iterations):
    """Run the benchmark for differential expression analysis."""
    print("Starting DEG benchmark")
    
    # Create or open output CSV file
    file_exists = os.path.exists(csv_path) and os.path.getsize(csv_path) > 0
    with open(csv_path, 'a') as file:
        if not file_exists:
            file.write("test_case,test_step,time_ms,success,rows,cols\n")
        
        # Load dataset
        print(f"Loading dataset from {h5ad_path}")
        load_start = time.time()
        success_load = True
        try:
            adata = sc.read_h5ad(h5ad_path)
            dimensions = (adata.n_obs, adata.n_vars)
        except Exception as e:
            print(f"Failed to load dataset: {e}")
            success_load = False
            dimensions = None
        
        load_time = int((time.time() - load_start) * 1000)
        log_benchmark_result(file, "load", load_time, success_load, dimensions)
        
        if success_load:
            # Preprocessing - equivalent to the normalize_expression and log1p_expression in Rust
            print("Preprocessing data")
            sc.pp.normalize_total(adata, target_sum=1e4)
            sc.pp.log1p(adata)
            
            print(f"Dataset loaded and preprocessed: {adata.n_obs} cells x {adata.n_vars} genes")
            print(f"Running DEG computation for {iterations} iterations")
            
            # Run iterations
            for i in range(1, iterations + 1):
                print(f"  Iteration {i}/{iterations}")
                result = run_deg_analysis(adata, i)
                log_benchmark_result(
                    file, 
                    result['step'], 
                    result['duration_ms'],
                    result['success'], 
                    dimensions
                )
            
            print("DEG benchmark completed successfully")

def main():
    parser = argparse.ArgumentParser(
        description="Benchmark for DEG (Differential Expression Gene) operations"
    )
    parser.add_argument('-i', '--input', required=True, help='Path to H5AD file')
    parser.add_argument('-o', '--output', required=True, help='Path for output CSV file')
    parser.add_argument('-n', '--iterations', type=int, default=5, help='Number of benchmark iterations to run')
    
    args = parser.parse_args()
    
    print("Starting DEG benchmark with:")
    print(f"  Input file: {args.input}")
    print(f"  Output file: {args.output}")
    print(f"  Iterations: {args.iterations}")
    
    run_deg_benchmark(args.input, args.output, args.iterations)
    
    print("Benchmark completed")

if __name__ == "__main__":
    # Set up parallel processing similar to the Rust code
    sc.settings.n_jobs = cpu_count()
    main()