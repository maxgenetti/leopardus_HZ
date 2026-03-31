# LAI and Selection

Local ancestry inference and ancestry-specific selection for the *L. geoffroyi* × *L. guttulus* hybrid zone. The LAI pipeline runs first and produces the panels and ancestry estimates that feed into the Selection pipeline.

---

## LAI

**`vcf2aHMM.py`**
Converts a per-chromosome LD-pruned VCF to Ancestry_HMM panel format, filtering for ancestry informative markers (AIMs) based on allele frequency difference between parental panels, minimum inter-site distance, and admixed sample coverage. Run once per chromosome; output is concatenated before running Ancestry_HMM.

```
Input:  <chrom>.final.vcf.gz
        <sample>.pop  (parent1/parent2/admixed assignment file)
Output: <sample>.<chrom>.panel  (written to stdout)
Args:   --vcf, --pop, --rate <cM/Mb>, --dist <min_bp>, --minDif, --minGT, --minDP, -geno
```

**`run.py`**
Runs a 7-step iterative Ancestry_HMM pipeline for all samples. Generates panels via `vcf2aHMM.py`, runs Ancestry_HMM, updates admixture fractions from posteriors, identifies and masks putatively introgressed sites in parental panels, and repeats. Steps 3 and 5 perform parental site masking (posterior > 0.1 of opposite-ancestry origin). Step 6 runs the final ancestry estimation and Viterbi decoding. Step 7 compiles per-sample summary statistics. Use `--step` to resume from a specific step. Parallelism is handled via GNU parallel.

```
Input:  chroms.txt, admixed_samples.txt, geoffroyi_samples.txt, guttulus_samples.txt
        <vcf_dir>/<chrom>.final.vcf.gz
Output: <out_dir>/<run>/step6/<sample>.posterior  (final per-site ancestry posteriors)
        <out_dir>/<run>/stepV/<sample>.posterior  (Viterbi-decoded tracts)
        <out_dir>/<run>/summary.step6.csv         (ancestry proportion, generations since admixture)
Args:   --chroms, --directory <run_name>, --out-dir, --vcf-dir,
        --samples, --parent1, --parent2,
        --threads, --dist, --diff, --minGT, --step <1-7>
```

Admixed individuals are then stratified into cohorts A–D by *L. geoffroyi* ancestry proportion (A: 0.10–0.40, B: 0.40–0.70, C: 0.70–0.95, D: 0.90–0.98) using `summary.step6.csv`. AHMM-S is run on cohorts A and B.

---

## Selection

**`mergePanels.py`**
Merges per-sample Ancestry_HMM panel files into a single cohort panel for AHMM-S input. Produces both a forward panel (P1 = *L. geoffroyi*) and a reversed panel (P1 = *L. guttulus*) for bidirectional selection inference.

```
Input:  <panel_dir>/<sample>.panel files
        sample list file (one sample per line)
Output: <out_dir>/merged.<label>.panel
        <out_dir>/merged.<label>.panelR  (reversed)
        <out_dir>/merged.<label>.sample
Args:   --file_path, --out_path, --sample_file, --label
```

**`runSelection.sh`**
Runs AHMM-S on all cohorts (A–D) and chromosomes in both forward and reverse directions. Selection coefficients are inferred between 0.001 and 0.5 using golden section search. Time since admixture and admixture fraction are taken from `run.py` step 6 estimates per cohort.

```
Input:  <out_dir>/merged.<cohort>.panel and .panelR
Output: <out_dir>/ahmms_merged.<cohort>.<chrom>-forw.out
        <out_dir>/ahmms_merged.<cohort>.<chrom>-rev.out
Usage:  bash runSelection.sh <out_dir>
```

**`findPeaks.py`**
Identifies significant selection peaks in AHMM-S log-likelihood ratio scores using a Bonferroni-corrected genome-wide threshold (α = 0.05). Peak boundaries are estimated from the LNL ratio profile. Optionally writes peak files and genome-wide plots. Peak output files feed directly into the GOEnrichment pipeline.

```
Input:  <in_dir>/ahmms_merged.<cohort>.<chrom>-<direction>.out
Output: <out_dir>/peaks.<cohort>.<direction>.out       (significant peaks)
        <out_dir>/peaks.<cohort>.<direction>.wide.out  (all candidate peaks)
        <out_dir>/peaks.<cohort>.<direction>.jpg       (genome-wide plots, with --write)
Args:   --in_dir, --out_dir, --populations <e.g. ABCD>, --directions forw/rev/both,
        --sel_cutoff, --write
```
