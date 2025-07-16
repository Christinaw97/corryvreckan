import uproot
import matplotlib.pyplot as plt
import numpy as np
import argparse
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

# Open ROOT file and extract histograms
with uproot.open(file_path) as file:
    hist_chi2 = file["Tracking4D/trackChi2ndof"]
    hist_clusters = file["Tracking4D/clustersPerTrack"]
    hist_tracks = file["Tracking4D/tracksPerEvent"]

    # Convert to numpy
    vals_chi2, edges_chi2 = hist_chi2.to_numpy()
    vals_clusters, edges_clusters = hist_clusters.to_numpy()
    vals_tracks, edges_tracks = hist_tracks.to_numpy()

    centers_chi2 = 0.5 * (edges_chi2[:-1] + edges_chi2[1:])
    centers_clusters = 0.5 * (edges_clusters[:-1] + edges_clusters[1:])
    centers_tracks = 0.5 * (edges_tracks[:-1] + edges_tracks[1:])

    # Compute stats
    def compute_stats(values, centers):
        if len(centers) == 0 or len(values)==0: return 0,0
        if np.sum(values) == 0: return 0,0
        mean = np.average(centers, weights=values)
        std = np.sqrt(np.average((centers - mean)**2, weights=values))
        return mean, std

    mean_chi2, std_chi2 = compute_stats(vals_chi2, centers_chi2)
    mean_clusters, std_clusters = compute_stats(vals_clusters, centers_clusters)
    mean_tracks, std_tracks = compute_stats(vals_tracks, centers_tracks)

# Create the 3-subplot figure
fig, axs = plt.subplots(3, 1, figsize=(10, 10), constrained_layout=True)
fig.suptitle("Tracking4D Summary Histograms", fontsize=16)

# Plot 1: trackChi2ndof
axs[0].bar(centers_chi2, vals_chi2, width=np.diff(edges_chi2), color='mediumseagreen', alpha=0.7)
axs[0].set_title("Track χ² / ndof")
axs[0].set_xlabel("χ² / ndof")
axs[0].set_ylabel("Entries")
axs[0].grid(True)
axs[0].set_xticks(np.linspace(edges_chi2[0], edges_chi2[-1], 11))  # 11 ticks spaced evenly
axs[0].text(0.95, 0.95,
            f"$\mu$ = {mean_chi2:.2e}\n$\sigma$ = {std_chi2:.2e}",
            transform=axs[0].transAxes,
            fontsize=10,
            ha='right', va='top',
            bbox=dict(boxstyle="round", facecolor="white", alpha=0.7))

# Plot 2: clustersPerTrack
axs[1].bar(centers_clusters, vals_clusters, width=np.diff(edges_clusters), color='cornflowerblue', alpha=0.7)
axs[1].set_title("Clusters per Track")
axs[1].set_xlabel("Number of Clusters")
axs[1].set_ylabel("Entries")
axs[1].grid(True)
axs[1].set_xticks(np.linspace(edges_clusters[0], edges_clusters[-1], 6))  # 11 ticks spaced evenly
axs[1].text(0.95, 0.95,
            f"$\mu$ = {mean_clusters:.2e}\n$\sigma$ = {std_clusters:.2e}",
            transform=axs[1].transAxes,
            fontsize=10,
            ha='right', va='top',
            bbox=dict(boxstyle="round", facecolor="white", alpha=0.7))

# Plot 3: tracksPerEvent
axs[2].bar(centers_tracks, vals_tracks, width=np.diff(edges_tracks), color='indianred', alpha=0.7)
axs[2].set_title("Track Multiplicity per Event")
axs[2].set_xlabel("Tracks per Event")
axs[2].set_ylabel("Entries")
axs[2].grid(True)
axs[2].set_xticks(np.linspace(edges_tracks[0], edges_tracks[-1], 11))  # 11 ticks spaced evenly
axs[2].text(0.95, 0.95,
            f"$\mu$ = {mean_tracks:.2e}\n$\sigma$ = {std_tracks:.2e}",
            transform=axs[2].transAxes,
            fontsize=10,
            ha='right', va='top',
            bbox=dict(boxstyle="round", facecolor="white", alpha=0.7))

# Save and show
plt.savefig(f"{BASE_PATH}/plots/trackings_{postfix}.png", format="png")
