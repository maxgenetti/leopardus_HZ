#!/usr/bin/env python3
"""
Identify and visualize significant selection peaks from AHMM-S output.
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy
from scipy.signal import find_peaks, peak_widths
import sys
import argparse
import os


CHROMS = {
    'A1': 239694388, 'B1': 205836458, 'C1': 222052948,
    'D1': 115783437, 'E1': 61591894,  'A2': 169709481,
    'B2': 152606360, 'C2': 158996348, 'D2': 88472993,
    'E2': 62031975,  'A3': 140630183, 'B3': 148130213,
    'C3': 153706148, 'D3': 94884351,  'E3': 42047855,
    'B4': 142838203, 'D4': 94461142
}

OFFSET = {
    'A1': 0,          'B1': 239694388,  'C1': 445530846,
    'D1': 667583794,  'E1': 783367231,  'A2': 844959125,
    'B2': 1014668606, 'C2': 1167274966, 'D2': 1326271314,
    'E2': 1414744307, 'A3': 1476776282, 'B3': 1617406465,
    'C3': 1765536678, 'D3': 1919242826, 'E3': 2014127177,
    'B4': 2056175032, 'D4': 2199013235, 'g':  2293474377
}


def get_args():
    parser = argparse.ArgumentParser(
        description='Identify significant selection peaks from AHMM-S output'
    )
    parser.add_argument('--in_dir', type=str, required=True,
                        help='Directory containing AHMM-S output files')
    parser.add_argument('--out_dir', type=str, required=True,
                        help='Directory for output files')
    parser.add_argument('--populations', type=str, default='ABCD',
                        help='Populations to analyze (e.g. ABCD or A)')
    parser.add_argument('--directions', type=str, default='both',
                        choices=['both', 'forw', 'rev'],
                        help='Directions to analyze (default: both)')
    parser.add_argument('--write', action='store_true', default=False,
                        help='Write output files and plots (default: print to stdout)')
    parser.add_argument('--sel_cutoff', type=float, default=0.05,
                        help='Selection coefficient cutoff (default: 0.05)')
    return parser.parse_args()


def get_x_ticks():
    ticks = [[], []]
    for key, value in OFFSET.items():
        if key in CHROMS:
            ticks[0].append(value + CHROMS[key] / 2)
            ticks[1].append(key)
    return ticks


def count_sites(in_dir, population, direction):
    total = 0
    for chrom in CHROMS:
        path = f"{in_dir}/ahmms_merged.{population}.{chrom}-{direction}.out"
        total += sum(1 for _ in open(path))
    return total


def process_population(population, direction, in_dir, out_dir, write_file, sel_cutoff):
    x_ticks = get_x_ticks()
    plt.figure(figsize=(40, 10))
    plt.axis('off')

    num_sites = count_sites(in_dir, population, direction)
    lnl_cutoff = scipy.stats.chi2.ppf(1 - (0.05 / num_sites), 1) / 2
    plt.title(f'{population} - {direction}')

    peaks = []
    peaks_plot = []
    peaks_wide = []
    full_lnl = []
    full_sel = []
    full_pos = []

    for chrom in CHROMS:
        pos_values, gen_pos, sel_values, lnl_values = [], [], [], []

        with open(f"{in_dir}/ahmms_merged.{population}.{chrom}-{direction}.out") as f:
            for line in f:
                l = line.split()
                pos_values.append(int(l[0]))
                gen_pos.append(int(l[0]) + OFFSET[chrom])
                sel_values.append(float(l[1]))
                lnl_values.append(max(float(l[2]), 0))

        lnl_prom = 0.5 * np.array(lnl_values)
        p1, p2 = find_peaks(lnl_values, height=lnl_cutoff,
                             prominence=lnl_prom, distance=100, width=51)

        for p in range(len(p1)):
            diff_a = int((pos_values[int(p2["left_ips"][p]) + 1] -
                          pos_values[int(p2["left_ips"][p])]) * (p2["left_ips"][p] % 1))
            diff_b = int((pos_values[int(p2["right_ips"][p]) + 1] -
                          pos_values[int(p2["right_ips"][p])]) * (p2["right_ips"][p] % 1))
            a = int(p2["left_ips"][p])
            b = int(p2["right_ips"][p])

            peaks_wide.append([chrom, pos_values[p1[p]], sel_values[p1[p]],
                                lnl_values[p1[p]], pos_values[a] + diff_a,
                                pos_values[b] + diff_b])

            if lnl_values[p1[p]] > scipy.stats.chi2.ppf(1 - (0.05 / num_sites), 1):
                prom = p2['prominences'][p]
                lnl = lnl_values[p1[p]]
                rh = min(1 - (lnl_cutoff - (lnl - prom)) / prom, 5000)
                pw = peak_widths(lnl_values, [p1[p]], rel_height=rh)
                diff_c = int((pos_values[int(pw[2][0]) + 1] -
                               pos_values[int(pw[2][0])]) * (pw[2][0] % 1))
                diff_d = int((pos_values[int(pw[3][0]) + 1] -
                               pos_values[int(pw[3][0])]) * (pw[3][0] % 1))
                c = int(pw[2][0])
                d = int(pw[3][0])

                peaks.append([chrom, pos_values[p1[p]], sel_values[p1[p]],
                               lnl_values[p1[p]], pos_values[a] + diff_a,
                               pos_values[b] + diff_b, pos_values[c] + diff_c,
                               pos_values[d] + diff_d])
                peaks_plot.append([chrom, gen_pos[p1[p]], sel_values[p1[p]],
                                    lnl_values[p1[p]], gen_pos[a] + diff_a,
                                    gen_pos[b] + diff_b, gen_pos[c] + diff_c,
                                    gen_pos[d] + diff_d])

        full_sel += sel_values
        full_lnl += lnl_values
        full_pos += gen_pos

        for i in range(1, 3):
            plt.subplot(2, 1, i)
            plt.axvline(x=OFFSET[chrom], color='k', linestyle='-')

    if write_file:
        os.makedirs(out_dir, exist_ok=True)

        with open(f"{out_dir}/peaks.{population}.{direction}.out", "w") as f:
            for peak in peaks:
                f.write(f'{peak[0]}\t{peak[6]}\t{peak[1]}\t{peak[7]}\t{peak[2]:.3f}\t{peak[3]}\n')

        with open(f"{out_dir}/peaks.{population}.{direction}.wide.out", "w") as f:
            for peak in peaks_wide:
                f.write(f'{peak[0]}\t{peak[4]}\t{peak[1]}\t{peak[5]}\t{peak[2]:.3f}\t{peak[3]}\n')

        plt.subplot(2, 1, 1)
        plt.plot(full_pos, full_lnl, color='#0000FF')
        plt.xlim(0, OFFSET['g'])
        plt.ylim(0, max(full_lnl + [0.05]) * 1.2)
        plt.xticks(x_ticks[0], x_ticks[1], fontsize=24)
        plt.ylabel(f'LNL ratio - {population} {direction}')
        plt.axhline(y=lnl_cutoff, color='#FF0000', linestyle='-')
        for peak in peaks_plot:
            plt.scatter(peak[1], peak[3], color='#FF0000')

        plt.subplot(2, 1, 2)
        plt.xticks(x_ticks[0], x_ticks[1], fontsize=24)
        plt.plot(full_pos, full_sel, color='#0000FF')
        plt.xlim(0, OFFSET['g'])
        plt.ylim(0, max(full_sel + [0.05]) * 1.2)
        plt.ylabel(f'Selection coefficient - {population} {direction}')
        plt.axhline(y=sel_cutoff, color='#FF0000', linestyle='-')
        for peak in peaks_plot:
            plt.scatter(peak[1], peak[2], color='#FF0000')

        plt.subplots_adjust(left=.1, bottom=.1, right=.95, top=.95, wspace=.01, hspace=0.1)
        plt.savefig(f"{out_dir}/peaks.{population}.{direction}.jpg", dpi=300)

    else:
        print(f"{population}-{direction} : {len(peaks)} peaks")
        for peak in peaks:
            print(f'{peak[0]}\t{peak[6]}\t{peak[1]}\t{peak[7]}\t{peak[2]:.3f}\t{peak[3]:.3f}')
        print(f"Wide {population}-{direction} : {len(peaks_wide)} peaks")
        for peak in peaks_wide:
            print(f'{peak[0]}\t{peak[4]}\t{peak[1]}\t{peak[5]}\t{peak[2]:.3f}\t{peak[3]:.3f}')


def main():
    args = get_args()
    populations = list(args.populations.upper())
    directions = ['forw', 'rev'] if args.directions == 'both' else [args.directions]

    for population in populations:
        for direction in directions:
            process_population(population, direction, args.in_dir,
                               args.out_dir, args.write, args.sel_cutoff)


if __name__ == "__main__":
    main()
