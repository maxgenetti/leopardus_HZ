# Script used to create hzar input for mtDNA or average genetic proportion
# Trindade F.T., 2024.

import pandas as pd
from numpy import median

def generate_hzar_input(input_file, output_file, target_attr):
    df = pd.read_csv(input_file, sep=";")
    df = df.dropna(subset=[df.columns[2]])
    df.columns = ["ID", "DIST", "ATTR"]
    df = df.sort_values(by="DIST").reset_index(drop=True)
    

    results = []
    i = 0
    while i < len(df):
        start_idx = i
        window_ids = [df.loc[start_idx, "ID"]]
        window_distances = [df.loc[start_idx, "DIST"]]
        window_attr = [df.loc[start_idx, "ATTR"]]

        # Expand the location range (km) to improve sampling
        for max_dist in [1,5,10,20]:
            for j in range(i + 1, len(df)):
                if df.loc[j, "DIST"] - df.loc[start_idx, "DIST"] < max_dist:
                    window_ids.append(df.loc[j, "ID"])
                    window_distances.append(df.loc[j, "DIST"])
                    window_attr.append(df.loc[j, "ATTR"])
                else:
                    break
            if len(window_ids) >= 2:
                break

        # Calculate FREQ and SAMPLING
        if target_attr == "freq":         
            window_attr = list(map(float, window_attr))
            sampling_count = len(window_ids) * 2
            freq = f"{round(median(window_attr), 6):.6f}"

        else:
            sampling_count = len(window_ids)
            target_count = window_attr.count(target_attr)
            freq = target_count / sampling_count

        # Update ID and DIST
        combined_id = "_".join(window_ids)
        mean_dist = round(median(window_distances), 2)

        results.append({
            "ID": combined_id,
            "DIST": mean_dist,
            "FREQ": float(freq),
            "SAMPLING": sampling_count})

        i += len(window_ids)


    output_df = pd.DataFrame(results)
    output_df["DIST"] = output_df["DIST"].map(lambda x: f"{x:.3f}")
    output_df["FREQ"] = output_df["FREQ"].map(lambda x: f"{x:.3f}")
    output_df.to_csv(output_file, sep=",", index=False)
    print(f"Output file saved to: {output_file}")


input_file = "posterior.csv"
output_file = "file.csv" # this will be the hzar input file
target_attr = "freq"  # mtDNA attribute or "freq"

generate_hzar_input(input_file, output_file, target_attr)
