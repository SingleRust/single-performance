#!/bin/bash

# Hard-coded path to the binary
BINARY_PATH="../../rs-methods/target/release/deg-performance"

# Define separate CSV files for timing and memory
TIMING_CSV="../../timings/rust/deg_performance.csv"
MEMORY_CSV="../../timings/rust/deg_memory.csv"

# Check if binary exists
if [ ! -f "$BINARY_PATH" ]; then
    echo "Error: Binary not found at $BINARY_PATH"
    exit 1
fi

if [ ! -f "$MEMORY_CSV" ]; then
    echo "dataset,n_tries,max_memory_kb" > "$MEMORY_CSV"
fi

# Define datasets and tries
declare -A datasets=(
    ["/local/06-24_single-rust-project/performance_eval/datasets/real_world/tahoe/tahoe_10K_cells.h5ad"]=12
    ["/local/06-24_single-rust-project/performance_eval/datasets/real_world/tahoe/tahoe_50K_cells.h5ad"]=12
    ["/local/06-24_single-rust-project/performance_eval/datasets/real_world/tahoe/tahoe_100K_cells.h5ad"]=12
    ["/local/06-24_single-rust-project/performance_eval/datasets/real_world/tahoe/tahoe_250K_cells.h5ad"]=12
    ["/local/06-24_single-rust-project/performance_eval/datasets/real_world/tahoe/tahoe_500K_cells.h5ad"]=8
    ["/local/06-24_single-rust-project/performance_eval/datasets/real_world/tahoe/tahoe_750K_cells.h5ad"]=8
    ["/local/06-24_single-rust-project/performance_eval/datasets/real_world/tahoe/tahoe_1M_cells.h5ad"]=3
    ["/local/06-24_single-rust-project/performance_eval/datasets/real_world/tahoe/tahoe_2M_cells.h5ad"]=3
    ["/local/06-24_single-rust-project/performance_eval/datasets/real_world/tahoe/tahoe_5M_cells.h5ad"]=3
    ["/local/06-24_single-rust-project/performance_eval/datasets/real_world/tahoe/tahoe_7M_cells.h5ad"]=3
    ["/local/06-24_single-rust-project/performance_eval/datasets/real_world/tahoe/tahoe_10M_cells.h5ad"]=2
    ["/local/06-24_single-rust-project/performance_eval/datasets/real_world/tahoe/tahoe_12M_cells.h5ad"]=2
    ["/local/06-24_single-rust-project/performance_eval/datasets/real_world/tahoe/tahoe_15M_cells.h5ad"]=1
    ["/local/06-24_single-rust-project/performance_eval/datasets/real_world/tahoe/tahoe_streaming_17M_cells.h5ad"]=1
    ["/local/06-24_single-rust-project/performance_eval/datasets/real_world/tahoe/tahoe_streaming_20M_cells.h5ad"]=1
    ["/local/06-24_single-rust-project/performance_eval/datasets/real_world/tahoe/tahoe_22M_cells.h5ad"]=1
    ["/local/06-24_single-rust-project/performance_eval/datasets/real_world/tahoe/tahoe_25M_cells.h5ad"]=1
    ["/local/06-24_single-rust-project/performance_eval/datasets/real_world/tahoe/tahoe_27M_cells.h5ad"]=1
    ["/local/06-24_single-rust-project/performance_eval/datasets/real_world/tahoe/tahoe_30M_cells.h5ad"]=1
    # Add more datasets as needed
)

# Process each dataset
for dataset in "${!datasets[@]}"; do
    n_tries=${datasets[$dataset]}
    echo "Processing dataset: $dataset with $n_tries tries"

    # Run the command and capture memory usage with /usr/bin/time
    # Note: We pass the TIMING_CSV to the binary as it records its own timing data
    start_time=$(date +%s%N)
    memory_output=$(/usr/bin/time -v "$BINARY_PATH" -i "$dataset" -o "$TIMING_CSV" -n "$n_tries" 2>&1)
    exit_code=$?
    end_time=$(date +%s%N)

    # Check if the command executed successfully
    if [ $exit_code -ne 0 ]; then
        echo "Error running the command. Exit code: $exit_code"
        echo "Output: $memory_output"
        continue
    fi

    # Calculate execution time in milliseconds (this is just for display, the program records its own timing)
    execution_time=$((($end_time - $start_time) / 1000000))

    # Extract maximum memory usage
    max_memory=$(echo "$memory_output" | grep "Maximum resident set size" | awk '{print $6}')
    if [ -z "$max_memory" ]; then
        echo "Warning: Could not extract memory usage information"
        max_memory="N/A"
    fi

    echo "  Max memory usage: $max_memory KB"
    echo "  Shell-measured execution time: $execution_time ms"

    # Append to memory CSV
    echo "$dataset,$n_tries,$max_memory" >> "$MEMORY_CSV"
done

echo "All datasets processed."
echo "Timing data written to $TIMING_CSV"
echo "Memory usage data written to $MEMORY_CSV"