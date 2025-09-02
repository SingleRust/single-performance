import os
import time
import pandas as pd
import scanpy as sc
import argparse
from pathlib import Path
import gc


def parse_args():
    parser = argparse.ArgumentParser(
        description="Benchmark for H5AD normalization and log1p transformation",
        prog="normalization_log1p_benchmark"
    )
    parser.add_argument(
        "-i", "--input",
        dest="h5ad_path",
        type=Path,
        required=True,
        help="Path to H5AD file"
    )
    parser.add_argument(
        "-o", "--output",
        dest="csv_path",
        type=Path,
        required=True,
        help="Path for output CSV file"
    )
    parser.add_argument(
        "-n", "--iterations",
        dest="iterations",
        type=int,
        default=5,
        help="Number of benchmark iterations to run"
    )
    parser.add_argument(
        "-t", "--target",
        dest="normalization_target",
        type=int,
        default=10000,
        help="Normalization target (counts per cell)"
    )
    return parser.parse_args()


def log_benchmark_result(file_path, test_case, step, duration_ms, success, dimensions=None):
    """Log benchmark results to CSV file."""
    
    file_exists = os.path.isfile(file_path) and os.path.getsize(file_path) > 0
    
    df = pd.DataFrame({
        "test_case": [test_case],
        "test_step": [step],
        "time_ms": [duration_ms],
        "success": [success],
        "rows": [dimensions[0] if dimensions else None],
        "cols": [dimensions[1] if dimensions else None]
    })
    
    df.to_csv(file_path, mode='a', header=not file_exists, index=False)
    

def run_normalization_benchmark(h5ad_path, csv_path, iterations, normalization_target):
    """Run the normalization + log1p benchmark."""
    
    print("Starting normalization + log1p benchmark")
    
    print(f"Loading dataset from {h5ad_path}")
    load_start = time.time()
    
    try:
        adata = sc.read_h5ad(h5ad_path, backed='r')
        load_time = int((time.time() - load_start) * 1000)
        load_success = True
        dimensions = (adata.n_obs, adata.n_vars)
        
        print(f"Dataset loaded: {adata.n_obs} cells x {adata.n_vars} genes")
        
        log_benchmark_result(csv_path, "normalization_log1p", "load", load_time, load_success, dimensions)
        
        print(f"Running normalization + log1p benchmarks for {iterations} iterations")
        del adata
        gc.collect()
        
        for i in range(1, iterations + 1):
            print(f"  Iteration {i}/{iterations}")
            
            adata = sc.read_h5ad(h5ad_path)
            
            combined_start = time.time()
            try:
                # First normalize
                sc.pp.normalize_total(adata, target_sum=normalization_target)
                
                # Then perform log1p transform
                sc.pp.log1p(adata)
                
                combined_success = True
            except Exception as e:
                print(f"Error during processing: {e}")
                combined_success = False
            
            del adata
            gc.collect()
            
            combined_time = int((time.time() - combined_start) * 1000)
            log_benchmark_result(
                csv_path, "normalization_log1p", "process", combined_time, combined_success, dimensions
            )
            
        print("Benchmark completed successfully")
        
    except Exception as e:
        load_time = int((time.time() - load_start) * 1000)
        load_success = False
        print(f"Failed to load dataset: {e}")
        log_benchmark_result(csv_path, "normalization_log1p", "load", load_time, load_success)


def main():
    args = parse_args()
    
    print("Starting normalization + log1p benchmark with:")
    print(f"  Input file: {args.h5ad_path}")
    print(f"  Output file: {args.csv_path}")
    print(f"  Iterations: {args.iterations}")
    print(f"  Normalization target: {args.normalization_target}")
    
    run_normalization_benchmark(
        args.h5ad_path, 
        args.csv_path, 
        args.iterations, 
        args.normalization_target
    )
    
    print("Benchmark completed")


if __name__ == "__main__":
    main()