import glob
import subprocess
import os
from concurrent.futures import ThreadPoolExecutor

file_list = glob.glob("../alignment/output/out_MIMOSA_*root")
file_list.sort()

max_workers = 4

with ThreadPoolExecutor(max_workers=max_workers) as executor:
    # submit jobs inline without a separate function
    futures = []
    for file_path in file_list:
        base_name = os.path.basename(file_path)
        postfix = os.path.splitext(base_name)[0]
        print(f"Running make_plots.sh on {file_path} with postfix {postfix}")
        futures.append(
            executor.submit(
                subprocess.run, ["bash", "./make_plots.sh", file_path, postfix]
            )
        )
    # wait for all futures to finish (optional)
    for future in futures:
        future.result()
