import pandas as pd
import subprocess
import argparse
import os

def get_args():
    parser = argparse.ArgumentParser(
        description='Merge per-sample AHMM-S panel files into a single cohort panel'
    )
    parser.add_argument('--file_path', type=str, required=True,
                        help='Directory containing per-sample .panel files')
    parser.add_argument('--out_path', type=str, required=True,
                        help='Output directory')
    parser.add_argument('--sample_file', type=str, required=True,
                        help='File listing sample IDs to include (one per line)')
    parser.add_argument('--label', type=str, default='Parental',
                        help='Panel label (default: Parental)')
    return parser.parse_args()


def main():
    args = get_args()
    file_path = args.file_path
    out_path = args.out_path
    label = args.label

    os.makedirs(out_path, exist_ok=True)

    sample_file = open(f"{out_path}/merged.{label}.sample", 'w')
    columns = {"Chromosome": str, "Position": int, "P1_A": int, "P1_a": int, "P2_A": int, "P2_a": int}
    column_order = ["Chromosome", "Position", "P1_A", "P1_a", "P2_A", "P2_a", "Morgans"]

    # Read shared sites across all samples
    command = (
        f"ls {file_path}*.panel | grep -f {args.sample_file} | "
        f"xargs cat | cut -f 1-6 | sort -k1,1 -k2,2n -u > {out_path}/admixed.{label}.sites"
    )
    subprocess.run(command, shell=True, check=True)

    panel = pd.read_csv(f"{out_path}/admixed.{label}.sites", sep="\t", header=None)
    panel = panel.rename(columns={0: "Chromosome", 1: "Position", 2: "P1_A",
                                   3: "P1_a", 4: "P2_A", 5: "P2_a"})

    file_list = []
    with open(args.sample_file) as f:
        for line in f:
            sample = line.rstrip()
            sample_file.write(f"{sample}\t2\n")
            file_list.append((sample, f"{file_path}/{sample}.panel"))
    sample_file.close()

    for sample, f in file_list:
        sample_panel = pd.read_csv(f, sep="\t", header=None)
        sample_panel.columns = ["Chromosome", "Position", "P1_A", "P1_a",
                                  "P2_A", "P2_a", "Distance",
                                  f"{sample}_A", f"{sample}_a"]
        columns[f"{sample}_A"] = int
        columns[f"{sample}_a"] = int
        column_order += [f"{sample}_A", f"{sample}_a"]
        sample_panel = sample_panel.drop(columns=["P1_A", "P1_a", "P2_A", "P2_a", "Distance"])
        panel = pd.merge(panel, sample_panel, on=['Chromosome', 'Position'], how='left')

    panel.fillna(0, inplace=True)

    # Add morgans column
    positions = [0] + panel["Position"].to_list()
    morgans = [(max(positions[i + 1] - positions[i], 1)) * 0.00000002
               for i in range(len(positions) - 1)]
    panel['Morgans'] = morgans

    panel = panel[column_order]
    panel = panel.astype(columns)

    # Write forward panel
    panel.to_csv(f"{out_path}/merged.{label}.panel", sep="\t", header=False, index=False)

    # Write reversed panel (P2 before P1) for reverse-direction AHMM-S runs
    column_order_rev = column_order[0:2] + column_order[4:6] + column_order[2:4] + column_order[6:]
    panel[column_order_rev].to_csv(f"{out_path}/merged.{label}.panelR", sep="\t", header=False, index=False)


if __name__ == "__main__":
    main()
