import argparse
import uproot
import matplotlib.pyplot as plt
import numpy as np
import os
import subprocess
parser = argparse.ArgumentParser(description="Read and display the contents of a file.")
parser.add_argument("--input_file", help="Path to the input file")
parser.add_argument("--postfix", help="Postfix to add to plot name")

args = parser.parse_args()
file_path = args.input_file
postfix = args.postfix

BASE_PATH = subprocess.Popen(['git', 'rev-parse', '--show-toplevel'], stdout=subprocess.PIPE).communicate()[0].rstrip().decode('utf-8')

if not os.path.isfile(file_path):
    print(f"Error: '{file_path}' is not a valid file.")
    exit()


# Create a 2-row, 6-column canvas
fig, axs = plt.subplots(2, 6, figsize=(20, 6), constrained_layout=True)
fig.suptitle("Global Residuals for MIMOSA Planes", fontsize=16)
# Open the ROOT file and the Tracking4D tree
with uproot.open(file_path) as file:
    tree = file["Tracking4D"]
    for i in range(6):
        hist_x = tree[f"MIMOSA26_{i}"]["global_residuals"]["GlobalResidualsX"]
        hist_y = tree[f"MIMOSA26_{i}"]["global_residuals"]["GlobalResidualsY"]

       # Convert to numpy
        values_x, edges_x = hist_x.to_numpy()
        values_y, edges_y = hist_y.to_numpy()
        centers_x = 0.5 * (edges_x[:-1] + edges_x[1:])
        centers_y = 0.5 * (edges_y[:-1] + edges_y[1:])

        # Calculate mean and std for X
        if len(centers_x) == 0: mean_x, std_x = 0, 0
        elif np.sum(values_x) == 0: mean_x, std_x = 0, 0
        else:
            mean_x = np.average(centers_x, weights=values_x)
            std_x = np.sqrt(np.average((centers_x - mean_x) ** 2, weights=values_x))

        # Plot X residuals (top row)
        ax_x = axs[0, i]
        ax_x.bar(centers_x, values_x, width=np.diff(edges_x), alpha=0.7, color='skyblue')
        ax_x.set_title(f"MIMOSA26_{i} - X")
        ax_x.set_xlabel("Residual X [mm]")
        ax_x.set_ylabel("Entries")
        ax_x.grid(True)
        ax_x.text(0.05, 0.95,
                  f"$\mu$ = {mean_x:.2e}\n$\sigma$ = {std_x:.2e}",
                  transform=ax_x.transAxes,
                  fontsize=10,
                  verticalalignment='top',
                  bbox=dict(boxstyle="round", facecolor="white", alpha=0.7))

        # Calculate mean and std for Y
        if len(values_y) or np.sum(values_y) == 0: mean_y, std_y = 0, 0
        else:
            mean_y = np.average(centers_y, weights=values_y)
            std_y = np.sqrt(np.average((centers_y - mean_y) ** 2, weights=values_y))

        # Plot Y residuals (bottom row)
        ax_y = axs[1, i]
        ax_y.bar(centers_y, values_y, width=np.diff(edges_y), alpha=0.7, color='salmon')
        ax_y.set_title(f"MIMOSA26_{i} - Y")
        ax_y.set_xlabel("Residual Y [mm]")
        ax_y.set_ylabel("Entries")
        ax_y.grid(True)
        ax_y.text(0.05, 0.95,
                  f"$\mu$ = {mean_x:.2e}\n$\sigma$ = {std_x:.2e}",
                  transform=ax_y.transAxes,
                  fontsize=10,
                  verticalalignment='top',
                  bbox=dict(boxstyle="round", facecolor="white", alpha=0.7))

plt.savefig(f"{BASE_PATH}/plots/global_residuals_{postfix}.png", format="png")

