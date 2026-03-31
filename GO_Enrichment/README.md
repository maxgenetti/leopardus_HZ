# GOEnrichment

Permutation-based GO term enrichment analysis for ancestry-specific selection peaks identified by AHMM-S.

---

**`runAnnotate.sh`**
Runs `annotateGenes.py` for all directions (forw, rev) and GO categories (bio, mol, loc) in parallel.

```
Usage:  bash runAnnotate.sh <pop> <annotations_dir>
```

**`annotateGenes.py`**
Matches each selection peak to nearby genes using an interval tree over GO annotation ranges. Produces three gene lists per run differing in how genes are assigned to peaks: all genes within 5 Mb of the peak center, the single closest gene, and all genes overlapping the measured peak window.

```
Input:  <pop>/peaks.B.<direction>.out
        go_ranges.<go>.txt
Output: <pop>/real-<direction>.<go>.0.txt  (all genes within 5 Mb)
        <pop>/real-<direction>.<go>.1.txt  (closest gene only)
        <pop>/real-<direction>.<go>.2.txt  (all genes within peak window)
Args:   --pop, --d forw/rev, --go bio/mol/loc, --annotations
```

---

**`runSimulate.sh`**
Prints `simulateNull.py` job commands to stdout for cluster submission. Run 10 blocks per condition.

```
Usage:  bash runSimulate.sh <pop> <annotations_dir> | bash
        bash runSimulate.sh <pop> <annotations_dir> | xargs -I {} sbatch --wrap="{}"
```

**`simulateNull.py`**
Generates a null distribution for GO enrichment by randomly placing peaks across the AIM site universe, preserving peak count and width. Results are stored as pickles and read by `compareEnrichment.py`.

```
Input:  <pop>/*.panel files (defines AIM site universe)
        <pop>/peaks.B.<direction>.out
        go_ranges.<go>.txt
Output: ./pickles/<pop>-<direction>.<go>.<window>.<block>.pkl
Args:   -p <pop>, -i forw/rev/both, -o bio/mol/loc, -s <n_sims>,
        -g <block_index>, -d <window_kb>, -n (nearest gene only), --annotations
```

---

**`runCompare.sh`**
Runs `compareEnrichment.py` for all directions and GO categories. Pass `--window-only` to run windowed analyses (100, 500, 2500 kb) instead of the default peak-width analysis.

```
Usage:  bash runCompare.sh <pop> <annotations_dir> [--window-only]
```

**`compareEnrichment.py`**
Compares observed GO term counts in peaks to the null distributions from `simulateNull.py` and computes enrichment p-values. Reads the gene lists from `annotateGenes.py` and the pickled simulations, then writes enrichment TSVs to `./bio/`, `./mol/`, and `./loc/`.

```
Input:  <pop>/real-<direction>.<go>.<suffix>.txt
        ./pickles/<pop>-<direction>.<go>.<window>.<block>.pkl
        go_ranges.<go>.txt
Output: ./<go>/<pop>-<direction>.<go>.<window>.tsv
Args:   --pop, --d forw/rev/both, --go bio/mol/loc, --window <kb>,
        --dup <n_blocks>, -n (nearest gene only), --annotations
```

---

**`runCorrection.sh`**
Runs `fdr.py` on all enrichment TSVs in bio/, mol/, and loc/.

```
Usage:  bash runCorrection.sh <enrichment_dir>
```

**`fdr.py`**
Applies Benjamini-Hochberg and Bonferroni correction to enrichment p-values. Writes two output files per method: `.full` (all significant terms) and `.subset` (terms observed in more than one peak).

```
Input:  <go>/<pop>-<direction>.<go>.<window>.tsv
Output: <go>/<pop>-<direction>.<go>.<window>.bh.full
        <go>/<pop>-<direction>.<go>.<window>.bh.subset
        <go>/<pop>-<direction>.<go>.<window>.bf.full
        <go>/<pop>-<direction>.<go>.<window>.bf.subset
Args:   <tsv_path_without_extension> (positional)
```

---

**`cleanResults.sh`**
Removes all intermediate TSV files from GO category subdirectories.

```
Usage:  bash cleanResults.sh <enrichment_dir>
```
