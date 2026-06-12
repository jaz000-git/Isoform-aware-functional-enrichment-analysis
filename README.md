# Isoform-Aware Functional Enrichment Analysis Pipeline

This repository provides a pipeline for **isoform-level differential transcript usage (DTU) classification** followed by **GO term over-representation analysis (ORA)**. The method is designed to capture functional enrichment driven specifically by **isoform-level variation**, rather than gene-level changes.

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
  --go_file data/IPS_go_map_clean.txt
```

---

## Input Files

### DTU file (`--dtu`)

A CSV file containing differential transcript usage results. The following columns are required:

| Column     | Description                   |
| ---------- | ----------------------------- |
| isoform_id | Ensembl transcript identifier |
| padj       | FDR adjusted p-value          |
| logFC      | Log 2 fold-change             |

Example:

```csv
isoform_id,padj,logFC
ENST00000335137,0.001,1.2
ENST00000423372,0.04,-0.8
```

### GO mapping file (`--go_file`)

A tab-separated file containing transcript-to-GO mappings.

Required columns:

| Column     | Description           |
| ---------- | --------------------- |
| transcript | Ensembl transcript identifier |
| gene       | Ensembl gene identifier       |
| GO         | Gene Ontology term    |

Example:

```text
GO	transcript	gene
GO:0090336	ENST00000496770	ENSG00000160097
```

---

## Output Files

Running the pipeline produces three output files:

### `{prefix}_dtu.tsv`

Classified transcript table containing:

| Column     | Description           |
| ---------- | --------------------- |
| transcript | Transcript identifier |
| condition  | cond, ref, or non_sig |

### `{prefix}_cond.tsv`

GO enrichment results for transcripts enriched in the condition of interest.

### `{prefix}_ref.tsv`

GO enrichment results for transcripts enriched in the reference condition.

Enrichment output columns:

| Column          | Description                            |
| --------------- | -------------------------------------- |
| GO              | GO term                                |
| fg_hits         | Foreground transcript hits             |
| bg_hits         | Background transcript hits             |
| TranscriptRatio | Foreground hit proportion              |
| BgRatio         | Background hit proportion              |
| OR              | Odds ratio                             |
| p               | Raw p-value                            |
| FDR             | Benjamini–Hochberg adjusted p-value    |
| test            | Statistical test used (chi2 or fisher) |
| genes           | Genes contributing to enrichment       |

---

## Classification Criteria

Transcripts are classified according to the following thresholds:

* **cond**: `padj < 0.05` and `logFC > 0.5`
* **ref**: `padj < 0.05` and `logFC < -0.5`
* **non_sig**: all other transcripts

---

## Statistical Methods

For each GO term, enrichment is assessed using a 2×2 contingency table comparing foreground and background transcript sets.

* Chi-square test is used when expected counts are sufficient.
* Fisher's exact test is used when any expected cell counts are less than 5.
* Odds ratios are calculated using a Haldane correction.
* Multiple testing correction is performed using the Benjamini–Hochberg false discovery rate (FDR) procedure.

---

## Notes

* Transcript version suffixes (e.g. `.1`, `.2`) should be removed prior to analysis.

---

