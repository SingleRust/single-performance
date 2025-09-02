import os
import time
import argparse
import scanpy as sc
import pandas as pd
import numpy as np
import gc

def log_benchmark_result(file_handle, test_case, step, duration_ms, success, dimensions=None, filtered_count=None):
    """Log a benchmark result to the CSV file."""
    dims_str = f"{dimensions[0]},{dimensions[1]}" if dimensions else ","
    filtered_str = str(filtered_count) if filtered_count is not None else ""
    
    file_handle.write(f"{test_case},{step},{duration_ms},{success},{dims_str},{filtered_str}\n")

def run_filtering_benchmark(h5ad_path, csv_path, iterations):
    """Run the filtering benchmark on an h5ad file."""
    print("Starting filtering benchmark")
    
    # Create or open the output CSV file
    file_exists = os.path.isfile(csv_path)
    with open(csv_path, 'a') as file:
        # Add header if file is new
        if not file_exists:
            file.write("test_case,test_step,time_ms,success,rows,cols,filtered_count\n")
        
    
        if True:
            print(f"Running filtering benchmarks for {iterations} iterations")
            
            # Combined filtering benchmark
            for i in range(1, iterations + 1):
                print(f"  Combined filtering - Iteration {i}/{iterations}")
                
                adata_copy = sc.read_h5ad(h5ad_path)
                n_obs = adata_copy.n_obs
                n_vars = adata_copy.n_vars
                
                combined_start = time.time()
                try:
                    # Calculate QC metrics if they don't exist
                    if 'n_genes_by_counts' not in adata_copy.obs or 'n_cells_by_counts' not in adata_copy.var:
                        sc.pp.calculate_qc_metrics(adata_copy, inplace=True)
                    
                    initial_n_obs = adata_copy.n_obs
                    initial_n_vars = adata_copy.n_vars
                    
                    # Cell filtering with scanpy functions
                    sc.pp.filter_cells(adata_copy, min_genes=200)
                    sc.pp.filter_cells(adata_copy, max_genes=5000)
                    sc.pp.filter_cells(adata_copy, min_counts=1000)
                    
                    # Gene filtering with scanpy functions
                    sc.pp.filter_genes(adata_copy, min_cells=3)
                    sc.pp.filter_genes(adata_copy, min_counts=10)
                    max_cells = int(0.9 * adata_copy.n_obs)
                    sc.pp.filter_genes(adata_copy, max_cells=max_cells)
                    
                    # Calculate filtered counts
                    filtered_cells = initial_n_obs - adata_copy.n_obs
                    filtered_genes = initial_n_vars - adata_copy.n_vars
                    filtered_total = filtered_cells + filtered_genes
                    combined_success = True
                except Exception as e:
                    print(f"Combined filtering failed: {e}")
                    filtered_total = None
                    combined_success = False
                
                del adata_copy
                gc.collect()
                
                combined_time = int((time.time() - combined_start) * 1000)
                
                
                log_benchmark_result(
                    file,
                    "combined_filtering",
                    "filter",
                    combined_time,
                    combined_success,
                    dimensions=(n_obs, n_vars),
                    filtered_count=filtered_total
                )
            
            print("Filtering benchmarks completed successfully")

def main():
    # Parse command line arguments
    parser = argparse.ArgumentParser(
        description="Benchmark for H5AD filtering operations",
    )
    parser.add_argument('-i', '--input', dest='h5ad_path', required=True,
                        help='Path to H5AD file')
    parser.add_argument('-o', '--output', dest='csv_path', required=True,
                        help='Path for output CSV file')
    parser.add_argument('-n', '--iterations', dest='iterations', type=int, default=5,
                        help='Number of benchmark iterations to run')
    
    args = parser.parse_args()
    
    print("Starting filtering benchmark with:")
    print(f"  Input file: {args.h5ad_path}")
    print(f"  Output file: {args.csv_path}")
    print(f"  Iterations: {args.iterations}")
    
    run_filtering_benchmark(args.h5ad_path, args.csv_path, args.iterations)
    
    print("Benchmark completed")

if __name__ == "__main__":
    main()