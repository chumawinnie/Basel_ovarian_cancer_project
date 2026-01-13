# PhD Candidate Test: Single-Cell CNV Analysis

## Context

You are provided with **single-cell RNA-seq data from High-Grade Serous Ovarian Cancer (HGSOC)**. For each patient, tumor cells were collected from:

- **Solid tumor tissue**
- **Ascites (fluid) samples**

Non-malignant cells (e.g., immune cells) are included to serve as reference for copy number variation (CNV) inference.

The goal of this task is to evaluate your ability to perform **inferCNV analysis**, compare CNV profiles between sample types, and identify recurrent genomic alterations.

---

## Data Preparation

The dataset has been prepared as follows:

1. Original single-cell RNA-seq data were obtained from **Vasquez-Garcia et al.**
2. Only patients with cells available from **both Solid tumor tissue and Ascites** were retained.
3. For each patient and sample type:
   - Up to **1,500 tumor cells** were randomly selected (if available).
   - **1,000 immune cells** were randomly selected to serve as reference/baseline for inferCNV.
4. Cell metadata includes:
   - `cell_id`
   - `patient_id`
   - `tumor_site` (`Solid` / `Ascites`)
   - `cell_type` (`Malignant` / `Reference`)

This ensures a **balanced dataset per patient and sample type** for CNV comparison, while keeping the analysis feasible on a personal computer.

---

## Data Provided

1. A **SingleCellExperiment (SCE) object** containing:
   - Raw count matrix
   - Cell metadata
2. A gene annotation file: `gene_positions.tsv`, containing chromosome and genomic position information.

---

## Tasks

You are expected to perform the following:

### 1. InferCNV Analysis

- Use the `inferCNV` R package.
- The following parameters are suggested to keep the runtime manageable:
``` cutoff = 0.1,
    cluster_by_groups = TRUE,
    denoise = TRUE,
    HMM = FALSE,
    diagnostics = FALSE,
    plot_steps = FALSE,
    inspect_subclusters = FALSE,
    no_prelim_plot = TRUE
```

### 2. Compare CNV profiles between Solid and Ascites samples

- Compute CNV metrics per cell and/or per chromosome.
- Compare CNV profiles between Solid and Ascites samples.
- Visualize differences using appropriate plots (e.g., chromosome-level CNV line plots).

### 3. Identify frequently altered regions

- Identify genomic regions or chromosomes that are **most consistently amplified or deleted**
- Distinguish between **shared alterations** and **site-specific alterations**.
- You may use simple statistics (mean CNV, frequency of gain/loss per region) and/or visual inspection of plots

### 4. Short presentation

- Prepare a **short** presentation summarizing your analysis.

### 5. Code and Reproducibility

- All code used for the analysis **must be shared** via a Git repository (e.g., GitHub, GitLab).