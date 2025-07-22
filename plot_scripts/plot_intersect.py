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
parser.add_argument("--output_dir", help="Path to the output plots")
parser.add_argument("--postfix", help="Postfix to add to plot name")

args = parser.parse_args()
file_path = args.input_file
group = args.group
postfix = args.postfix
output_dir = args.output_dir

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

output_file = f"{BASE_PATH}/{output_dir}/{group.lower()}_global_intersect2D_{postfix}.png"
plt.savefig(output_file, format="png")
print(f"2D projection plot saved to {output_file}")


# ---------- ADDITION: Create both X and Y projections from TH2 and plot with stats ----------
def compute_stats(values, centers):
    if len(centers) == 0 or len(values) == 0: return 0, 0
    if np.sum(values) == 0: return 0, 0
    mean = np.average(centers, weights=values)
    std = np.sqrt(np.average((centers - mean) ** 2, weights=values))
    return mean, std

fig_proj, axs_proj = plt.subplots(2, 6, figsize=(24, 8), constrained_layout=True)

for i, hist in enumerate(histograms):
    ax_x = axs_proj[0, i]
    ax_y = axs_proj[1, i]

    if hist is None:
        ax_x.set_visible(False)
        ax_y.set_visible(False)
        continue

    try:
        vals, edges_x, edges_y = hist.to_numpy()

        # X projection
        proj_x = np.sum(vals, axis=1)
        centers_x = 0.5 * (edges_x[:-1] + edges_x[1:])
        mean_x, std_x = compute_stats(proj_x, centers_x)

        ax_x.bar(centers_x, proj_x, width=np.diff(edges_x), color='dodgerblue', alpha=0.8)
        ax_x.set_title(f"MIMOSA26_{i} X Projection")
        ax_x.set_xlabel("X [mm]")
        ax_x.set_ylabel("Entries")
        ax_x.grid(True)
        ax_x.text(0.95, 0.95,
                  f"$\mu$ = {mean_x:.2f} mm\n$\sigma$ = {std_x:.2f} mm",
                  transform=ax_x.transAxes,
                  fontsize=10,
                  ha='right', va='top',
                  bbox=dict(boxstyle="round", facecolor="white", alpha=0.7))

        # Y projection
        proj_y = np.sum(vals, axis=0)
        centers_y = 0.5 * (edges_y[:-1] + edges_y[1:])
        mean_y, std_y = compute_stats(proj_y, centers_y)

        ax_y.bar(centers_y, proj_y, width=np.diff(edges_y), color='darkorange', alpha=0.8)
        ax_y.set_title(f"MIMOSA26_{i} Y Projection")
        ax_y.set_xlabel("Y [mm]")
        ax_y.set_ylabel("Entries")
        ax_y.grid(True)
        ax_y.text(0.95, 0.95,
                  f"$\mu$ = {mean_y:.2f} mm\n$\sigma$ = {std_y:.2f} mm",
                  transform=ax_y.transAxes,
                  fontsize=10,
                  ha='right', va='top',
                  bbox=dict(boxstyle="round", facecolor="white", alpha=0.7))

    except Exception as e:
        print(f"Error projecting histogram {i}: {e}")
        ax_x.set_visible(False)
        ax_y.set_visible(False)

fig_proj.suptitle(f"{group} Global Intersect 1D Projections with Stats", fontsize=22)

output_file_proj = f"{BASE_PATH}/{output_dir}/{group.lower()}_global_intersect1D_xy_{postfix}.png"
fig_proj.savefig(output_file_proj, format="png")
print(f"1D X/Y projection plot with stats saved to {output_file_proj}")
