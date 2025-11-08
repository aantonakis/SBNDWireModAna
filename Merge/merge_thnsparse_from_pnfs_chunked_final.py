import os
import sys
import subprocess
import ROOT
import tempfile
import shutil
from concurrent.futures import ThreadPoolExecutor, as_completed
import math
import array


# --- Configuration ---
#pnfs_dir = "/pnfs/sbnd/scratch/users/aantonak/THnSparseFiles"
pnfs_dir = sys.argv[1]

pattern = ".root"
#merged_name = "merged_sparse.root"
merged_name = sys.argv[2]

pnfs_output_dir = "/pnfs/sbnd/scratch/users/aantonak/MergedResults"
n_threads = 64          # number of concurrent ifdh cp calls
chunk_size = 100       # number of files per chunk

object_name = "hHit"  # name of your THnSparseD object
#object_name = "hTrack"  # name of your THnSparseD object

dims = array.array('i', [0, 3, 4])

# --- Create local working directory ---
local_dir = tempfile.mkdtemp(prefix="thnsparse_merge_")
print(f"📁 Created temporary working directory: {local_dir}")

# --- Step 1: list PNFS files ---
print("Listing PNFS directory...")
cmd_list = ["ifdh", "ls", pnfs_dir]
pnfs_files = subprocess.check_output(cmd_list, text=True).splitlines()
pnfs_files = [f for f in pnfs_files if f.endswith(pattern)]
n_files = len(pnfs_files)
print(f"Found {n_files} ROOT files on PNFS.")

if n_files == 0:
    raise RuntimeError("No input files found!")

n_chunks = math.ceil(n_files / chunk_size)
print(f"Splitting into {n_chunks} chunks of up to {chunk_size} files each.\n")

# --- Helper: copy one file ---
def copy_file(pnfs_path):
    try:
        subprocess.check_call(["ifdh", "cp", pnfs_path, local_dir])
        return pnfs_path, True
    except subprocess.CalledProcessError:
        return pnfs_path, False

# --- Helper: merge one chunk of local files ---
def merge_local_files(files, output_name):
    merged = [None for num in range(6)]
    for fname in files:
        f = ROOT.TFile.Open(fname)
        if not f or f.IsZombie():
            print(f"⚠️ Skipping unreadable file {fname}")
            continue
        for num in range(6):
            h = f.Get(object_name+str(num))
       
            if not h:
                print(f"⚠️ No {object_name} found in {fname}")
                continue
            if merged[num] is None:
                merged[num] = h.Clone(output_name+str(num))
                #merged.SetDirectory(0)
                merged[num].Reset()
                merged[num].Add(h)
           
            else:
                merged[num].Add(h)
        f.Close()
    return merged

# --- Step 2: process in chunks ---
partial_files = []
final_merged = [None for num in range(6)]

for chunk_idx in range(n_chunks):
    start = chunk_idx * chunk_size
    end = min((chunk_idx + 1) * chunk_size, n_files)
    chunk_files = pnfs_files[start:end]
    #print(f"🔹 Processing chunk {chunk_idx+1}/{n_chunks} ({len(chunk_files)} files)")

    # Copy in parallel
    with ThreadPoolExecutor(max_workers=n_threads) as executor:
        futures = {executor.submit(copy_file, f): f for f in chunk_files}
        for i, fut in enumerate(as_completed(futures), 1):
            path, ok = fut.result()
            status = "✅" if ok else "⚠️"
            if not ok:
                print(f"  [{i}/{len(chunk_files)}] {status} {os.path.basename(path)}")

    # Merge local files in this chunk
    local_files = [os.path.join(local_dir, f) for f in os.listdir(local_dir) if f.endswith(".root")]
    #print(f"  Merging {len(local_files)} local files...")
    merged_chunk = merge_local_files(local_files, f"merged_chunk_{chunk_idx}")

    if not merged_chunk[0]:
        print(f"  ⚠️ Skipping empty chunk {chunk_idx}")
        shutil.rmtree(local_dir)
        local_dir = tempfile.mkdtemp(prefix="thnsparse_merge_")
        continue

    if chunk_idx == 0:
        for num in range(6):
            final_merged[num] = merged_chunk[num].Projection(len(dims), dims)
        #final_merged.SetDirectory(0)
    else:
        for num in range(6):
            h_temp = merged_chunk[num].Projection(len(dims), dims)
            final_merged[num].Add(h_temp)

    # Write partial output
    #partial_path = os.path.join(local_dir, f"partial_{chunk_idx}.root")
    #fout = ROOT.TFile(partial_path, "RECREATE")
    #merged_chunk.Write()
    #fout.Close()
    #partial_files.append(partial_path)
    #print(f"  ✅ Wrote partial result: {partial_path}")

    # Cleanup local files before next chunk
    for f in os.listdir(local_dir):
        if f.startswith("partial_"):
            continue
        os.remove(os.path.join(local_dir, f))
    print("  🧹 Cleaned up chunk temporary files.\n")


# --- Step 5: Write final result locally ---
local_out = os.path.join(local_dir, merged_name)
fout = ROOT.TFile(local_out, "RECREATE")
#proj.Write()
for num in range(6):
    final_merged.Write()

fout.Close()
print(f"✅ Saived final merged projection: {local_out}")

# --- Step 6: Copy result back to PNFS ---
print(f"Uploading merged file to PNFS: {pnfs_output_dir}")
subprocess.call(["ifdh", "mkdir_p", pnfs_output_dir])
subprocess.call(["ifdh", "cp", "-D", local_out, pnfs_output_dir])

print(f"✅ Merge complete and copied to PNFS: {pnfs_output_dir}/{merged_name}")

# --- Step 7: Cleanup ---
shutil.rmtree(local_dir)
print(f"🧹 Removed temporary directory: {local_dir}")

