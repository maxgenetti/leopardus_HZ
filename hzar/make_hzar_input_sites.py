# Script to create HZAR input files based on genomic site posteriors.
# Trindade F.T., 2025.

import pandas as pd
from numpy import median
import os

# Step 1: Generate Range from BED
def generate_range(bed_file, posterior_file, window):
    """
    Create the range of sites of which create HZAR input files.
    The range is based on given position (bed_file), window size (window) and chromosome boundaries.
    """
    bed = pd.read_csv(bed_file, sep=" ", header=None, names=["chrom", "start", "end"])
    window = int(window)
    
    # Read the first posterior file to determine chromosome boundaries
    posterior = pd.read_csv(posterior_file, sep="\t")
    chrom_max_pos = posterior.groupby("chrom")["position"].max().to_dict()

    ranges = []
    for _, row in bed.iterrows():
        chrom, start, end = row["chrom"], row["start"], row["end"]
        start = max(start - window, 1)
        end = min(end + window, chrom_max_pos.get(chrom, end))
        ranges.append((chrom, start, end))
    return ranges

# Step 2: Filter Posterior Files
def filter_posterior_files(ranges, posterior_dir, samples, primary_column):
    """
    Filter posterior files to contain only information for target range of sites.
    Get the average posterior for the required ancestry (primary_column).
    """
    freq_data = []
    primary_column = str(primary_column)

    for sample in samples:
        posterior_file = os.path.join(posterior_dir, f"{sample}.posterior")
        if not os.path.exists(posterior_file):
            print(f"Warning: Posterior file for {sample} not found.")
            continue
        
        posterior = pd.read_csv(posterior_file, sep="\t")

        # Calculate average posterior of one ancestry
        for chrom, start, end in ranges:
            filtered = posterior[(posterior["chrom"] == chrom) & (posterior["position"].between(start, end))]
            if primary_column == "2,0":
                filtered["FREQ"] = round(filtered["2,0"] + (filtered["1,1"] / 2), 6)
            elif primary_column == "0,2":
                filtered["FREQ"] = round(filtered["0,2"] + (filtered["1,1"] / 2), 6)
            filtered = filtered[["chrom", "position", "FREQ"]]
            filtered["sample"] = sample
            freq_data.append(filtered)

    return pd.concat(freq_data) if freq_data else pd.DataFrame()

# Step 3: Average Site Position
def average_site_pos(freq_data, site_dist=1000):
    """
    Adjust genomic positions in freq_data by averaging positions within a given distance (site_dist).
    If positions are closer than site_dist, they will be averaged.
    """
    freq_data = freq_data.sort_values(by=["chrom", "position"]).reset_index(drop=True)
    i = 0
    while i < len(freq_data):
        current_chrom = freq_data.loc[i, "chrom"]
        group = [i]
        for j in range(i + 1, len(freq_data)):
 
            if freq_data.loc[j, "chrom"] != current_chrom:
                break
            # If positions are close enough, group them together
            if abs(freq_data.loc[j, "position"] - freq_data.loc[i, "position"]) <= site_dist:
                group.append(j)
            else:
                break

        # If more than one position in the group, average it and update the df
        if len(group) > 1:
            avg_position = round(freq_data.loc[group, "position"].mean())
            freq_data.loc[group, "position"] = avg_position

        i = group[-1] + 1

    return freq_data

# Step 4: Calculate median FREQ in windows
def calculate_mean_freq(transect, freq_data):
    transect = transect[transect["ID"].isin(freq_data["sample"])]
    freq_data = freq_data[freq_data["sample"].isin(transect["ID"])]
    transect = transect.sort_values("DIST").reset_index(drop=True)

    results = []
    i = 0
    while i < len(transect):
        start_idx = i
        window_ids = [transect.loc[start_idx, "ID"]]
        window_distances = [transect.loc[start_idx, "DIST"]]
        window_freqs = [freq_data[freq_data["sample"] == transect.loc[start_idx, "ID"]]["FREQ"].mean()]

        # Expand the location range (km) to improve sampling
        for max_dist in [1,5,10,20]:
            for j in range(i + 1, len(transect)):
                if transect.loc[j, "DIST"] - transect.loc[start_idx, "DIST"] < max_dist:
                    window_ids.append(transect.loc[j, "ID"])
                    window_distances.append(transect.loc[j, "DIST"])
                    window_freqs.append(freq_data[freq_data["sample"] == transect.loc[j, "ID"]]["FREQ"].mean())
                else:
                    break
            if len(window_ids) >= 2:
                break

        mean_freq = f"{round(median(window_freqs), 6):.6f}"
        mean_dist = round(median(window_distances), 2)
        sampling = len(window_ids) * 2

        results.append({"ID": "_".join(window_ids), "DIST": mean_dist, "FREQ": mean_freq, "SAMPLING": sampling})
        i += len(window_ids)
    return pd.DataFrame(results)

# Step 5: Generate final HZAR input files
def generate_files(transect_file, freq_data, output_dir, dist_window=None):
    """
    Generate hzar input files for each position and save them to the output directory.
    Only includes IDs with DIST within the specified dist_window (if provided).
    """
    transect = pd.read_csv(transect_file, sep=";")
    for position in freq_data["position"].unique():
        position_freqs = freq_data[freq_data["position"] == position]
        results = calculate_mean_freq(transect, position_freqs)

        output_file = os.path.join(output_dir, f"{position}.hzar")
        if results.empty:
            print(f"No results for position {position}")
        else:
            if dist_window:
                            before = results[results["DIST"] < dist_window[0]]
                            after = results[results["DIST"] > dist_window[1]]
                            within = results[(results["DIST"] >= dist_window[0]) & (results["DIST"] <= dist_window[1])]

                            high_freq = results[results["FREQ"].astype(float) > 0.95]
                            low_freq = results[results["FREQ"].astype(float) < 0.05]

                            if (
                                len(before) >= 3 and
                                len(after) >= 3 and
                                len(within) >= 3 and
                                len(high_freq) >= 2 and
                                len(low_freq) >= 2
                            ):
                                results.to_csv(output_file, sep=",", index=False)
                                print(f"Saved: {output_file}")
                            else:
                                print(f"Excluded {position}: does not meet filtering criteria")

# Main Function
def main(bed_file, window, sample_list_file, transect_file, posterior_dir, primary_column, output_dir, site_dist=1000, dist_window=None):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    samples = pd.read_csv(sample_list_file, header=None).squeeze("columns").tolist()
    posterior_file = os.path.join(posterior_dir, f"{samples[0]}.posterior")

    # Step 1: Generate Range from BED
    ranges = generate_range(bed_file, posterior_file, window)

    # Step 2: Filter posterior files
    freq_data = filter_posterior_files(ranges, posterior_dir, samples, primary_column)

    # Step 3: Average Site Position
    freq_data = average_site_pos(freq_data, site_dist)

    # Step 4 and 5: Calculate median freq and create final HZAR input files
    generate_files(transect_file, freq_data, output_dir, dist_window)


# Example usage
main(
    bed_file="region.bed", 
    window="500000",
    sample_list_file="samples.txt", 
    transect_file="transect.csv", 
    posterior_dir="dir/", 
    primary_column="2,0", 
    output_dir="dir/",
    site_dist=1000,
    dist_window=(1500, 2500),
)
