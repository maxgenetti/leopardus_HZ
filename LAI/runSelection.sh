#!/usr/bin/env bash
# Run AHMM-S for ancestry-specific selection inference
# Runs both forward (L. geoffroyi ancestry) and reverse (L. guttulus ancestry) directions
# Usage: bash runAHMMSS.sh <out_dir>

set -euo pipefail

if [ $# -lt 1 ]; then
    echo "Usage: bash runAHMMSS.sh <out_dir>"
    exit 1
fi

outDir="$1"

declare -A prop
prop[A]=0.941
prop[B]=0.856
prop[C]=0.574
prop[D]=0.203

declare -A gen
gen[A]=140.78
gen[B]=75.9482
gen[C]=69.1105
gen[D]=154.115

chroms="A1 A2 A3 B1 B2 B3 B4 C1 C2 C3 D1 D2 D3 D4 E1 E2 E3"

for P in A B C D; do
    p=${prop[$P]}
    q=$(echo "1 - $p" | bc)
    g=${gen[$P]}

    for chrom in $chroms; do
        # Forward direction: selection on L. geoffroyi ancestry
        ahmm-s \
            -i "${outDir}/merged.${P}.panel" \
            -s "${outDir}/merged.${P}.sample" \
            --chr "${chrom}" \
            -p 1 10000 "${p}" \
            -p 0 "${g}" "${q}" \
            --window m 0.01 \
            --gss 2 980000 1 0.001 0.5 \
            --ne 10000 \
            --full_selection_space \
            1>"${outDir}/ahmms_merged.${P}.${chrom}-forw.out" \
            2>"${outDir}/ahmms_merged.${P}.${chrom}-forw.err"

        # Reverse direction: selection on L. guttulus ancestry
        ahmm-s \
            -i "${outDir}/merged.${P}.panelR" \
            -s "${outDir}/merged.${P}.sampleR" \
            --chr "${chrom}" \
            -p 1 10000 "${q}" \
            -p 0 "${g}" "${p}" \
            --window m 0.01 \
            --gss 2 980000 1 0.001 0.5 \
            --ne 10000 \
            --full_selection_space \
            1>"${outDir}/ahmms_merged.${P}.${chrom}-rev.out" \
            2>"${outDir}/ahmms_merged.${P}.${chrom}-rev.err"
    done
done
