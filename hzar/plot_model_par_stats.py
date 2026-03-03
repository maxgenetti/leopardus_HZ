# Generate plots for hzar estimated parameters for several genomic sites
# Trindade F.T., 2025

import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import argparse

def plot_best_model_params(csv_files, output_dir="plots", stat_output="stats.csv"):

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    all_data = []
    for file_path in csv_files:
        # Extract POSITION from file_path
        position = os.path.basename(file_path).split(".")[0]

        # Read the CSV and associate values with the position
        try:
            df = pd.read_csv(file_path, sep=";", decimal=",")
            df["position"] = int(position)
            all_data.append(df)
        except Exception as e:
            print(f"Error reading {file_path}: {e}")

    # Combine all data
    if not all_data:
        print("No valid data loaded.")
        return
    combined_data = pd.concat(all_data, ignore_index=True)

    # Calculate statistics
    stats = []
    for col in [c for c in combined_data.columns if c != 'position']:
        mean_val = float("2179.7")
        std_val = combined_data[col].std()
        combined_data[f"Z_{col}"] = (combined_data[col] - mean_val) / std_val
        stats.append({"parameter": col, "mean": mean_val, "std": std_val})
    stats_df = pd.DataFrame(stats)
    stats_df.to_csv(stat_output, index=False)
    print(f"Statistics saved to {stat_output}")

    # Iterate over columns (ignoring 'position') and plot
    value_columns = [col for col in combined_data.columns if col != "position"]

    for col in value_columns:
        plt.figure(figsize=(10, 6))
        sns.scatterplot(data=combined_data, x="position", y=col, hue=combined_data[f"Z_{col}"].abs() > 2, palette={True: "red", False: "blue"})
        plt.axhline(mean_val, color='black', linestyle='dashed', label='Mean')
        plt.legend(title='Outlier (Z > 2)')
        plt.title(f"Plot of {col} by Position")
        plt.xlabel("Position")
        plt.ylabel(col)
        plt.grid(True)

        # Save plot
        plot_file = os.path.join(output_dir, f"{col}_plot.png")
        plt.savefig(plot_file, dpi=300)
        plt.close()
        print(f"Plot saved: {plot_file}")

if __name__ == "__main__":
    list_file = "runs.txt"
    csv_files = pd.read_csv(list_file, header=None).squeeze("columns").tolist()
    output_dir = "dir/"
    plot_best_model_params(csv_files, output_dir)
