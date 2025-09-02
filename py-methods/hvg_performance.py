import os
import time
import pandas as pd
import scanpy as sc
import argparse
from pathlib import Path
import gc


def parse_args():
    parser = argparse.ArgumentParser(
        description="Benchmark for H5AD operations using Scanpy",
        prog="scanpy_hvg_benchmark"
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
    return parser.parse_args()


def log_benchmark_result(file_path, step, duration_ms, success, dimensions=None):
    """Log benchmark results to CSV file."""
    
    file_exists = os.path.isfile(file_path) and os.path.getsize(file_path) > 0
    
    dims_str = f"{dimensions[0]},{dimensions[1]}" if dimensions else ","
    
    df = pd.DataFrame({
        "test_case": ["hvg"],
        "test_step": [step],
        "time_ms": [duration_ms],
        "success": [success],
        "rows": [dimensions[0] if dimensions else None],
        "cols": [dimensions[1] if dimensions else None]
    })
    
    df.to_csv(file_path, mode='a', header=not file_exists, index=False)
    

def run_hvg_benchmark(h5ad_path, csv_path, iterations):
    """Run the HVG benchmark with specified parameters."""
    
    print("Starting HVG benchmark")
    
    print(f"Loading dataset from {h5ad_path}")
    load_start = time.time()
    
    try:
        adata = sc.read_h5ad(h5ad_path)
        load_time = int((time.time() - load_start) * 1000) 
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        load_success = True
        dimensions = (adata.n_obs, adata.n_vars)
        adata_var = adata.var.copy()
        print(f"Dataset loaded: {adata.n_obs} cells x {adata.n_vars} genes")
        
        log_benchmark_result(csv_path, "load", load_time, load_success, dimensions)
        
        print(f"Running HVG computation for {iterations} iterations")
        
        for i in range(1, iterations + 1):
            print(f"  Iteration {i}/{iterations}")
            
            
            hvg_start = time.time()
            try:
                sc.pp.highly_variable_genes(
                    adata, 
                    flavor='seurat',
                    n_bins=20,
                    min_mean=0.0125,
                    max_mean=3,
                    min_disp=0.5,
                    span=0.3
                )
                hvg_success = True
            except Exception as e:
                print(f"HVG computation failed: {e}")
                hvg_success = False
            
            hvg_time = int((time.time() - hvg_start) * 1000) 
            log_benchmark_result(csv_path, "compute_hvg", hvg_time, hvg_success, dimensions)
            adata.var = adata_var.copy()
        print("HVG benchmark completed successfully")
        
    except Exception as e:
        load_time = int((time.time() - load_start) * 1000)  
        load_success = False
        print(f"Failed to load dataset: {e}")
        log_benchmark_result(csv_path, "load", load_time, load_success)
    

def main():
    args = parse_args()
    
    print("Starting HVG benchmark with:")
    print(f"  Input file: {args.h5ad_path}")
    print(f"  Output file: {args.csv_path}")
    print(f"  Iterations: {args.iterations}")
    
    run_hvg_benchmark(args.h5ad_path, args.csv_path, args.iterations)
    
    print("Benchmark completed")


if __name__ == "__main__":
    main()