import os
import time
import pandas as pd
import scanpy as sc
import numpy as np
import argparse
from pathlib import Path
import gc
import anndata as ad


def parse_args():
    parser = argparse.ArgumentParser(
        description="Benchmark for nearest neighbor computation using Scanpy",
        prog="scanpy_knn_benchmark"
    )
    parser.add_argument(
        "-i", "--input",
        dest="input_path",
        type=Path,
        required=True,
        help="Path to NPY file containing PCA-transformed data"
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


def log_benchmark_result(file_path, test_case, step, duration_ms, success, dimensions=None, info=None):
    """Log benchmark results to CSV file."""
    
    file_exists = os.path.isfile(file_path) and os.path.getsize(file_path) > 0
    
    dims_str = f"{dimensions[0]},{dimensions[1]}" if dimensions else ","
    
    df = pd.DataFrame({
        "test_case": [test_case],
        "test_step": [step],
        "time_ms": [duration_ms],
        "success": [success],
        "rows": [dimensions[0] if dimensions else None],
        "cols": [dimensions[1] if dimensions else None]
    })
    
    df.to_csv(file_path, mode='a', header=not file_exists, index=False)


def run_knn_benchmark(input_path, csv_path, iterations):
    """Run the KNN benchmark with specified parameters."""
    
    print("Starting KNN benchmark")
    
    # Load the PCA-transformed data once
    print(f"Loading dataset from {input_path}")
    load_start = time.time()
    
    try:
        adata = sc.read_h5ad(input_path)
        
        load_time = int((time.time() - load_start) * 1000)
        load_success = True
        dimensions = adata.shape
        
        print(f"Dataset loaded: {dimensions[0]} samples x {dimensions[1]} features")
        
        var_data = pd.read_csv("/local/06-24_single-rust-project/single_test/var_data.csv")
        highly_variable_mask = var_data['highly_variable'].values.astype(bool)
        adata.var['highly_variable'] = highly_variable_mask
        sc.pp.pca(adata, use_highly_variable=True, svd_solver='randomized', zero_center=False, n_comps=50, random_state=42)
        
        # Set scanpy settings
        sc.settings.n_jobs = -1  # Use all available cores
        
        log_benchmark_result(
            csv_path, 
            "knn", 
            "load", 
            load_time, 
            load_success, 
            dimensions
        )
        
        print(f"Running KNN benchmarks for {iterations} iterations")
        
        # Fixed parameters matching Rust benchmark
        k = 15
        n_pcs = 30
        
        for i in range(1, iterations + 1):
            print(f"  Iteration {i}/{iterations}")
            
            # Clear only the specific fields created by the neighbors function
            # Based on the scanpy source code analysis
            key_added = "neighbors"
            conns_key = "connectivities"
            dists_key = "distances"
            
            # Clear neighbor results from previous iteration
            if hasattr(adata, 'uns') and key_added in adata.uns:
                del adata.uns[key_added]
            
            if hasattr(adata, 'obsp'):
                adata.obsp.pop(dists_key, None)
                adata.obsp.pop(conns_key, None)
            
            # Clear X_diffmap from obsm if it exists
            if hasattr(adata, 'obsm') and "X_diffmap" in adata.obsm:
                del adata.obsm["X_diffmap"]
            
            # Force garbage collection to clear any cached objects
            gc.collect()
            
            # KNN computation
            knn_start = time.time()
            try:
                sc.pp.neighbors(
                    adata, 
                    n_neighbors=k,
                    n_pcs=n_pcs,
                    method='gauss', 
                    random_state=42
                )
                knn_success = True
                knn_time = int((time.time() - knn_start) * 1000)
                
            except Exception as e:
                print(f"    KNN computation failed: {e}")
                knn_success = False
                knn_time = int((time.time() - knn_start) * 1000)
            
            log_benchmark_result(
                csv_path,
                "knn",
                "compute",
                knn_time,
                knn_success,
                dimensions
            )
            
            # Force garbage collection after each iteration
            gc.collect()
        
        print("KNN benchmark completed successfully")
        
    except Exception as e:
        load_time = int((time.time() - load_start) * 1000)
        load_success = False
        print(f"Failed to load dataset: {e}")
        log_benchmark_result(csv_path, "knn", "load", load_time, load_success)


def main():
    args = parse_args()
    
    print("Starting KNN benchmark with:")
    print(f"  Input file: {args.input_path}")
    print(f"  Output file: {args.csv_path}")
    print(f"  Iterations: {args.iterations}")
    
    run_knn_benchmark(
        args.input_path,
        args.csv_path,
        args.iterations
    )
    
    print("Benchmark completed")


if __name__ == "__main__":
    main()