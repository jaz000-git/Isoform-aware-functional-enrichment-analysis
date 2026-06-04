# Isoform-Aware Functional Enrichment Analysis Pipeline

This repository provides a pipeline for **isoform-level differential transcript usage (DTU) classification** followed by **GO term over-representation analysis (ORA)**. The method is designed to capture functional enrichment driven specifically by **isoform switching events**, rather than gene-level changes.

---

## Overview

The pipeline performs three main steps:

1. **DTU classification**
   - Classifies transcripts into:
     - `cond` (condition-specific upregulated isoforms)
     - `ref` (reference-specific upregulated isoforms)
     - `non_sig` (not significant)

2. **GO annotation merging**
   - Maps transcripts to GO terms using a provided annotation file.

3. **GO enrichment analysis**
   - Performs over-representation analysis separately for `cond` and `ref`
   - Uses:
     - Chi-square test (default)
     - Fisher’s exact test (if expected counts < 5)
   - Computes:
     - Odds ratio (Haldane correction)
     - FDR (Benjamini–Hochberg correction)

---

## Example Run

```bash
python pipeline/isoform_ora_pipeline.py \
  --dtu data/DTU_DCM_vs_NF.csv \
  --prefix dcm \
  --go_file data/IPS_go_map.txt
```

## ⚠️ Work in Progress (WIP)

This repository is part of an ongoing thesis project. It is incomplete and subject to future refinement, but the pipeline is stable for the analyses presented in the thesis.

---
