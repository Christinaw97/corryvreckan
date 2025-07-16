#!/usr/bin/env python3

import uproot
import numpy as np
import matplotlib.pyplot as plt
import sys
from pathlib import Path
import subprocess
import argparse


cmap = plt.cm.viridis
new_cmap = cmap.copy()
new_cmap.set_under('white')  # any value below vmin will be white

BASE_PATH = subprocess.Popen(['git', 'rev-parse', '--show-toplevel'], stdout=subprocess.PIPE).communicate()[0].rstrip().decode('utf-8')


parser = argparse.ArgumentParser(description="Read and display the contents of a file.")
parser.add_argument("--input_file", help="Path to the input file")
parser.add_argument("--group", help="Tracking4D or Clustering4D")
parser.add_argument("--postfix", help="Postfix to add to plot name")

args = parser.parse_args()
file_path = args.input_file
group = args.group
postfix = args.postfix

if group not in ["Tracking4D", "Clustering4D"]:
    print("Error: group must be either 'Tracking4D' or 'Clustering4D'")
    sys.exit(1)
if group == "Tracking4D": hist_name = "global_intersect"
else: hist_name = "clusterPositionGlobal"

histograms = []

try:
    with uproot.open(file_path) as file:
        if group not in file:
            print(f"Error: Group '{group}' not found in the ROOT file.")
            sys.exit(1)

        tree = file[group]
        for i in range(6):
            key = f"MIMOSA26_{i}"
            if key not in tree:
                print(f"Warning: '{key}' not found in '{group}'. Skipping.")
                histograms.append(None)
                continue
            subdir = tree[key]
            if hist_name not in subdir:
                print(f"Warning: {hist_name} not found under '{group}/{key}'. Skipping.")
                histograms.append(None)
                continue
            histograms.append(subdir[hist_name])
except Exception as e:
    print(f"Error opening file or reading histograms: {e}")
    sys.exit(1)

fig, axs = plt.subplots(2, 3, figsize=(18, 10), constrained_layout=True)
axs = axs.flatten()

for i, (hist, ax) in enumerate(zip(histograms, axs)):
    if hist is None:
        ax.set_visible(False)
        continue

    vals, edges_x, edges_y = hist.to_numpy()
    X, Y = np.meshgrid(edges_x, edges_y)
    pcm = ax.pcolormesh(X, Y, vals.T, cmap=new_cmap, vmin=1e-3)

    ax.set_title(f"{group} MIMOSA26_{i} Global Intersect")
    ax.set_xlabel("X [mm]")
    ax.set_ylabel("Y [mm]")
    fig.colorbar(pcm, ax=ax, label='Entries')

plt.suptitle(f"{group} Global Intersect TH2F Histograms", fontsize=20)
from pathlib import Path

output_file = f"{BASE_PATH}/plots/{group.lower()}_global_intersect_{postfix}.png"
plt.savefig(output_file, format="png")
print(f"Plot saved to {output_file}")

