# AML WES Variant Processing Pipeline

Pipeline for processing and analyzing somatic variants from acute myeloid leukemia (AML) whole-exome sequencing data. Developed in the [Eisfeld Lab](https://cancer.osu.edu/find-a-researcher/search-researcher-directory/ann-kate-eisfeld) at The Ohio State University College of Medicine.

RNA-seq expression analysis and RAS pathway activation scoring live in the companion repository: [aml-rnaseq-analysis](https://github.com/ethanhamp/aml-rnaseq-analysis).

---

## What This Pipeline Does

1. **Ingests** per-sample somatic and germline variant files (Excel format from Varhouse/Mutect2)
2. **Filters** variants by read depth, VAF, and an artifact gene exclusion list
3. **Combines** samples across sequencing runs into a unified variant table
4. **Analyzes** mutation burden, VAF dynamics (diagnosis → relapse), oncoplots, and domain-level variant distributions

---

## Repository Structure

```
├── config/
│   ├── hypermutators.example.txt      # Template: patient IDs to flag as hypermutated
│   └── excluded_samples.example.txt   # Template: sample IDs to exclude from analysis
├── Exomes/
│   └── Scripts/
│       ├── build_master_variants.R    # Main pipeline: ingest, filter, combine all runs
│       ├── mutation_burden.R          # Dx vs Relapse mutation count and VAF analysis
│       ├── exomes_maftools.R          # MAF canonicalization, oncoplots, gene summaries
│       └── exome_filtering3.R         # VAF filtering, pairing, and delta computation
├── Variants/
│   └── Scripts/                       # ARHG domain localization and effect prediction
├── General/
│   └── artifact_genes.txt             # Curated list of recurrently artifact-prone genes
│                                      # excluded from somatic variant analysis (e.g.,
│                                      # genes with known mapping artifacts, common
│                                      # germline contamination, or sequencer noise)
└── cBioPortal/
    └── Scripts/                       # Public database mining across cancer studies
```

---

## Dependencies

All analysis is done in **R**. Install required packages before running:

```r
install.packages(c("tidyverse", "readxl", "writexl", "ggplot2", "ggrepel", "readr"))

# maftools via Bioconductor
if (!require("BiocManager")) install.packages("BiocManager")
BiocManager::install("maftools")
```

---

## Setup

### 1. Set your project root

The scripts resolve all file paths relative to a single `BASE_DIR`. Set this before running any script by adding the following to your `~/.Renviron` file:

```
ARHG_BASE_DIR=/path/to/your/project/root
```

Then restart R. Alternatively, set it at the top of any script with `Sys.setenv(ARHG_BASE_DIR = "/path/to/project")`.

### 2. Configure sample exclusions

Copy the template files and fill in your own IDs:

```bash
cp config/hypermutators.example.txt config/hypermutators.txt
cp config/excluded_samples.example.txt config/excluded_samples.txt
```

- `hypermutators.txt` — patient IDs with abnormally high variant burden to flag (but not remove) during processing
- `excluded_samples.txt` — sample IDs to drop entirely (e.g., duplicate timepoints, QC failures)

Both files support `#` comments. These files are gitignored and never leave your machine.

---

## Usage

### Build the master variant table

Reads per-sample somatic and germline Excel files across all configured sequencing runs, applies inclusion filters, and outputs a combined variant table per run and a single master file.

```r
Rscript Exomes/Scripts/build_master_variants.R
```

Edit the `RUN_CONFIG` block at the top of the script to point to your own sequencing run directories.

### Analyze mutation burden

Compares mutation counts and VAF distributions between diagnosis and relapse timepoints.

```r
Rscript Exomes/Scripts/mutation_burden.R
```

### Generate oncoplots

Converts the combined variant table to MAF format and produces oncoplots using `maftools`.

```r
Rscript Exomes/Scripts/exomes_maftools.R
```

---

## Variant Filtering Criteria

| Filter | Threshold |
|--------|-----------|
| Read depth (AD Total) | > 20 |
| Variant allele frequency | ≥ 2% |
| FILTER field | PASS only |
| Artifact genes | Excluded via `General/artifact_genes.txt` |
| Hypermutated patients | Flagged via `config/hypermutators.txt` |

---

## Input Data Format

Variant calls originate from **Nationwide Children's Hospital's Churchill WES pipeline** — a Mutect2-based somatic caller that produces per-sample Excel workbooks through the Varhouse reporting system. Churchill calls variants in tumor-normal pairs: the tumor sample is peripheral blood or bone marrow, and the matched germline normal is sorted B and T cells from the same patient. Standard FILTER flags are applied, and one Excel file per sample is exported containing:

- Somatic and germline calls in separate sheets
- Per-variant columns: Gene, HGVS notation (MANE transcript), Alt Percentage, AD Total, T-VAF, N-VAF, FILTER status, and functional effect annotation
- Sample metadata in the filename (patient ID, timepoint, run ID)

The `build_master_variants.R` script reads these Excel files directly — no VCF conversion required. It handles two column naming formats used across different Churchill pipeline versions automatically.

This pipeline is gene-agnostic — it processes all PASS variants across the exome. Downstream scripts can filter to any gene or gene set of interest.

---

## Contact

Ethan Hamp — [ethanhamp@gmail.com](mailto:ethanhamp@gmail.com)  
Eisfeld Lab, The Ohio State University College of Medicine
