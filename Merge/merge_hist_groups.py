import os
import subprocess
import ROOT
import tempfile
import shutil
import math
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed

def copy_file(pnfs_path, dest_dir):
    """Copy one file from PNFS to local dir using ifdh cp."""
    try:
        subprocess.check_call(["ifdh", "cp", pnfs_path, dest_dir])
        return pnfs_path, True
    except subprocess.CalledProcessError:
        return pnfs_path, False

def merge_histograms(files):
    """Merge all TH1 histograms from a list of ROOT files."""
    merged_hists = {}
    for fname in files:
        f = ROOT.TFile.Open(fname)
        if not f or f.IsZombie():
            print(f"⚠️ Skipping unreadable file {fname}")
            continue

        for key in f.GetListOfKeys():
            name = key.GetName()
            h = f.Get(name)
            if not h:
                continue

            if name not in merged_hists:
                merged_hists[name] = h.Clone()
                merged_hists[name].Reset()
            merged_hists[name].Add(h)

        f.Close()
    return merged_hists


def main():
    parser = argparse.ArgumentParser(description="Merge ROOT histograms into N output groups.")
    parser.add_argument("-p", "--pnfs_dir", required=True, help="Path to PNFS directory with ROOT files")
    parser.add_argument("-n", "--n_outputs", type=int, default=10, help="Number of merged output files to create")
    parser.add_argument("-o", "--output_dir", required=True, help="Destination PNFS directory for merged outputs")
    parser.add_argument("-t", "--threads", type=int, default=6, help="Number of parallel IFDH copies per batch")
    args = parser.parse_args()

    # Step 1. List PNFS files
    cmd_list = ["ifdh", "ls", args.pnfs_dir]
    pnfs_files = subprocess.check_output(cmd_list, text=True).splitlines()
    pnfs_files = [f for f in pnfs_files if f.endswith(".root")]
    n_files = len(pnfs_files)
    if n_files == 0:
        print("❌ No ROOT files found.")
        return
    print(f"📦 Found {n_files} ROOT files on PNFS.")

    # Step 2. Split files into groups
    n_outputs = args.n_outputs
    group_size = math.ceil(n_files / n_outputs)
    groups = [pnfs_files[i:i + group_size] for i in range(0, n_files, group_size)]
    print(f"🔹 Splitting into {len(groups)} groups (~{group_size} files each)\n")

    # Step 3. Create local working area
    local_dir = tempfile.mkdtemp(prefix="merge_groups_")
    print(f"📁 Using temporary dir: {local_dir}")

    # Step 4. Process each group
    for i_group, group in enumerate(groups):
        print(f"\n=== Merging group {i_group+1}/{len(groups)} ({len(group)} files) ===")

        # Copy files in parallel
        with ThreadPoolExecutor(max_workers=args.threads) as executor:
            futures = {executor.submit(copy_file, f, local_dir): f for f in group}
            for j, fut in enumerate(as_completed(futures), 1):
                path, ok = fut.result()
                status = "✅" if ok else "⚠️"
                print(f"  [{j}/{len(group)}] {status} {os.path.basename(path)}")

        local_files = [os.path.join(local_dir, f) for f in os.listdir(local_dir) if f.endswith(".root")]
        print(f"  🧮 Merging {len(local_files)} files in group {i_group}...")

        merged = merge_histograms(local_files)
        if not merged:
            print("  ⚠️ No histograms merged in this group, skipping.")
            continue

        out_path = os.path.join(local_dir, f"merged_{i_group}.root")
        fout = ROOT.TFile(out_path, "RECREATE")
        for name, hist in merged.items():
            hist.Write(name)
        fout.Close()

        print(f"  💾 Saved merged file: {out_path}")

        # Copy back to PNFS
        subprocess.call(["ifdh", "mkdir_p", args.output_dir])
        subprocess.call(["ifdh", "cp", "-D", out_path, args.output_dir])
        print(f"  🚀 Copied to PNFS: {args.output_dir}/merged_{i_group}.root")

        # Cleanup local files
        for f in os.listdir(local_dir):
            os.remove(os.path.join(local_dir, f))
        print("  🧹 Cleaned up local group files.")

    print("\n✅ All groups merged and uploaded.")
    shutil.rmtree(local_dir)
    print(f"🧹 Removed temporary directory: {local_dir}")

if __name__ == "__main__":
    main()

