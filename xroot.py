import sys
import os
import subprocess

input_list = sys.argv[1]
output_list = sys.argv[2]

with open(input_list) as fin, open(output_list, "w") as fout:
    for line in fin:
        filename = line.strip()
        if not filename:
            continue  # skip empty lines

        # Run pnfsToXRootD command and capture its output
        result = subprocess.run(
            ["pnfsToXRootD", filename],
            stdout=subprocess.PIPE,
            text=True,
            check=True
        )

        new_name = result.stdout.strip()
        fout.write(new_name + "\n")
        print(f"{filename} -> {new_name}")



print("done")

