#!/bin/bash

# Check if correct number of arguments are provided
if [ "$#" -ne 4 ]; then
    echo "Usage: $0 <files.list> <my_macro.C> <output_directory> <N files>"
    exit 1
fi

# Input arguments
FILES_LIST=$1
ROOT_MACRO=$2
OUTPUT_DIR=$3
N_FILES=$4

NUM_JOBS=4  # Number of parallel jobs



process_file() {
    local file=$1
    local out=$2
    #root -l -b -q "$ROOT_MACRO+(\"$file\")"
    timeout $timeout_duration root -l -b -q "$ROOT_MACRO(\"$file\", \"$out\")" || {
        echo "Error processing $file. Skipping to the next file..."
        continue  # Skip to the next file in case of an error
    }
}


# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Timeout duration in seconds
timeout_duration=120  # Adjust as needed (e.g., 2 minutes)

count=0

# Iterate over each file in files.list
while IFS= read -r ROOT_FILE; do
    if (( count > N_FILES )); then
        break  # Wait for all background jobs to finish before starting new ones
    fi
    # Get the base name of the ROOT file (without path and extension)
    BASE_NAME=$(basename "$ROOT_FILE" .root)
    
    # Define the output file name
    OUTPUT_FILE="$OUTPUT_DIR/${BASE_NAME}_output${count}.root"
    ((count+=1))    

    process_file "$ROOT_FILE" "$OUTPUT_FILE" &  # Run in the background
    # Run the ROOT macro on the current file and save to output directory
    #root -l -b -q "$ROOT_MACRO(\"$ROOT_FILE\", \"$OUTPUT_FILE\")"
    #timeout $timeout_duration root -l -b -q "$ROOT_MACRO(\"$ROOT_FILE\", \"$OUTPUT_FILE\")" || {
    #    echo "Error processing $file. Skipping to the next file..."
    #    continue  # Skip to the next file in case of an error
    #}

    if (( count % NUM_JOBS == 0 )); then
        wait  # Wait for all background jobs to finish before starting new ones
    fi


    # Check if the ROOT command was successful
    if [ $? -eq 0 ]; then
        echo "Processed $ROOT_FILE and saved to $OUTPUT_FILE"
    else
        echo "Error processing $ROOT_FILE"
    fi
done < "$FILES_LIST"

# Wait for any remaining background processes to finish
wait

