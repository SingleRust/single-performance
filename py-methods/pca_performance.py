import argparse
import time
import csv
import os
import numpy as np
import pandas as pd
import scanpy as sc
from pathlib import Path
from datetime import datetime

def log_benchmark_result(file_handle, test_case, step, duration_ms, success, dimensions=None):
    """Log benchmark results to CSV file"""
    dims_str = f"{dimensions[0]},{dimensions[1]}" if dimensions else ","
    file_handle.writerow([test_case, step, duration_ms, success, dims_str])

def run_pca_benchmark(h5ad_path, csv_path, iterations):
    """Run PCA benchmark using scanpy"""
    print("Starting PCA benchmark")
    
    # Create or open the output CSV file
    file_exists = os.path.exists(csv_path)
    with open(csv_path, 'a', newline='') as csvfile:
        writer = csv.writer(csvfile)
        
        # Add header if file is new
        if not file_exists or os.path.getsize(csv_path) == 0:
            writer.writerow(['test_case', 'test_step', 'time_ms', 'success', 'rows,cols'])
        
        # Load and prepare the dataset
        print(f"Loading dataset from {h5ad_path}")
        load_start = time.time()
        try:
            # Load h5ad file
            adata = sc.read_h5ad(h5ad_path)
            
            # Normalize to 10,000 counts per cell (equivalent to 1e4)
            sc.pp.normalize_total(adata, target_sum=1e4)
            
            # Log1p transformation
            sc.pp.log1p(adata)
            
            # Load highly variable genes from CSV
            var_data = pd.read_csv("/local/06-24_single-rust-project/single_test/var_data.csv")
            
            # Ensure the gene order matches
            if 'highly_variable' in var_data.columns:
                # Assuming the CSV has genes in the same order as adata.var
                highly_variable_mask = var_data['highly_variable'].values.astype(bool)
                adata.var['highly_variable'] = highly_variable_mask
            else:
                # If no highly variable column, compute it
                sc.pp.highly_variable_genes(adata, n_top_genes=2000)
            
            load_time = int((time.time() - load_start) * 1000)  # Convert to milliseconds
            load_success = True
            
            n_samples = adata.n_obs
            n_features = adata.n_vars
            n_hvg = adata.var['highly_variable'].sum()
            
        except Exception as e:
            print(f"Loading failed: {e}")
            load_time = int((time.time() - load_start) * 1000)
            load_success = False
            n_samples, n_features, n_hvg = None, None, None
            
        log_benchmark_result(writer, "pca", "load", load_time, load_success, 
                           (n_samples, n_features) if n_samples else None)
        
        if not load_success:
            print("Failed to load dataset, exiting")
            return
            
        print(f"Dataset loaded: {n_samples} samples x {n_features} features ({n_hvg} highly variable)")
        print(f"Running PCA benchmarks for {iterations} iterations")
        
        # PCA computation benchmark
        for i in range(1, iterations + 1):
            print(f"  Iteration {i}/{iterations}")
            
            # Clear any existing PCA results from previous iterations
            if 'X_pca' in adata.obsm:
                del adata.obsm['X_pca']
            if 'PCs' in adata.varm:
                del adata.varm['PCs']
            if 'pca' in adata.uns:
                del adata.uns['pca']
            
            pca_start = time.time()
            try:
                # Run PCA with 50 components using mask_var
                sc.pp.pca(adata, n_comps=50, svd_solver='randomized', 
                         mask_var='highly_variable', zero_center=False, random_state=42)
                
                pca_time = int((time.time() - pca_start) * 1000)
                pca_success = True
                
                # Verify PCA results exist
                if 'X_pca' in adata.obsm:
                    pca_shape = adata.obsm['X_pca'].shape
                    print(f"    PCA completed: {pca_shape[0]} samples x {pca_shape[1]} components")
                
            except Exception as e:
                print(f"    PCA computation failed: {e}")
                pca_time = int((time.time() - pca_start) * 1000)
                pca_success = False
            
            log_benchmark_result(writer, "pca", "compute", pca_time, pca_success,
                               (n_samples, n_features))
    
    print("PCA benchmark completed successfully")

def main():
    parser = argparse.ArgumentParser(
        description="Benchmark for PCA computation on single-cell data",
        epilog="Benchmark PCA computation on single-cell data using scanpy"
    )
    
    parser.add_argument('-i', '--input', dest='h5ad_path', type=Path, required=True,
                       help='Path to h5ad file containing the data')
    parser.add_argument('-o', '--output', dest='csv_path', type=Path, required=True,
                       help='Path for output CSV file')
    parser.add_argument('-n', '--iterations', type=int, default=5,
                       help='Number of benchmark iterations to run (default: 5)')
    
    args = parser.parse_args()
    
    print("Starting PCA benchmark with:")
    print(f"  Input file: {args.h5ad_path}")
    print(f"  Output file: {args.csv_path}")
    print(f"  Iterations: {args.iterations}")
    
    # Set scanpy settings for performance
    sc.settings.n_jobs = 128  # Use 128 cores for parallel operations where supported
    
    run_pca_benchmark(args.h5ad_path, args.csv_path, args.iterations)
    
    print("Benchmark completed")

if __name__ == "__main__":
    main()