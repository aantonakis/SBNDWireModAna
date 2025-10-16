import os
import subprocess
from pathlib import Path
from math import ceil
from concurrent.futures import ThreadPoolExecutor
import sys


def hadd(output, inputs, threads=8):
    """Run ROOT hadd with multithreading."""
    cmd = ["hadd", "-f", f"-j{threads}", str(output)] + [str(f) for f in inputs]
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

def merge_in_chunks(input_dir, output_dir, final_name, chunk_size=100, threads=8):
    input_dir, output_dir = Path(input_dir), Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    files = sorted(input_dir.glob("*.root"))
    n_files = len(files)
    n_chunks = ceil(n_files / chunk_size)
    print(f"Found {n_files} files → {n_chunks} chunks of {chunk_size}")

    batch_files = []

    # Stage 1: merge each chunk
    with ThreadPoolExecutor(max_workers=64) as executor:
        futures = []
        for i in range(n_chunks):
            start, end = i * chunk_size, (i + 1) * chunk_size
            chunk = files[start:end]
            batch_out = output_dir / f"merged_batch_{i:04d}.root"
            batch_files.append(batch_out)
            print(f"[{i+1}/{n_chunks}] Merging → {batch_out}")
            futures.append(executor.submit(hadd, batch_out, chunk, threads))
        for f in futures:
            f.result()

    # Stage 2: merge all batch outputs
    final_out = output_dir / final_name
    print(f"Merging {len(batch_files)} batch files into {final_out}")
    hadd(final_out, batch_files, threads)
    print("✅ Done:", final_out)

# Example usage
if __name__ == "__main__":
    input_dir = sys.argv[1]
    output_dir=sys.argv[2]
    output_name=sys.argv[3]
    merge_in_chunks(input_dir, output_dir, output_name, chunk_size=100, threads=8)

