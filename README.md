# inferCNV Analysis: Solid vs Ascites in HGSOC

Single-cell CNV analysis comparing solid tumor and ascites compartments in high-grade serous ovarian cancer (HGSOC), using inferCNV on scRNA-seq data from the SPECTRUM cohort.

## Dataset

- **Source**: Vasquez-Garcia et al. (*Nature*, 2022) — SPECTRUM cohort
- **Size**: 5,000 genes x 70,437 cells from 25 HGSOC patients
- **Composition**: 45,437 malignant cells (34,500 Solid + 10,937 Ascites) and 25,000 immune reference cells
- **Paired samples**: 19 patients with malignant cells from both tumor sites

## Analysis overview

1. **inferCNV** run with immune cells as diploid reference (cutoff=0.1, denoise=TRUE, HMM=FALSE)
2. Per-chromosome CNV scores computed for each cell, then aggregated by tumor site
3. Paired Wilcoxon signed-rank tests (n=19) with BH correction to compare Solid vs Ascites
4. Sensitivity analysis excluding 7 patients with <50 ascites malignant cells
5. Gene-level classification of shared vs site-specific alterations (gain >1.05, loss <0.95)

## Key findings

- **High overall correlation** (Pearson r = 0.94) between Solid and Ascites CNV profiles, consistent with shared clonal origin
- **Shared gains**: Chr20 (77% of genes), Chr8 (51%), Chr3 (41%)
- **Shared losses**: Chr22 (100% of genes), Chr6 (51% — HLA/MHC region)
- **Site-specific**: Chr20 is significantly more amplified in Ascites (p_adj = 0.045, n=19; confirmed in sensitivity analysis with p_adj = 0.038, n=12)
- Chr20 harbours metastasis-relevant genes including GNAS, LAMA5, WFDC2 (HE4), SOX18, MMP9, and UBE2C

## Output structure

```
Infer_copy_number_results/
├── infercnv_run/             # raw inferCNV output
├── figures/
│   ├── 01_genome_wide_cnv_profile.png
│   ├── 02_chr_barplot_comparison.png
│   ├── 03_chr_alteration_frequency.png
│   ├── 04a_heatmap_solid.png
│   ├── 04b_heatmap_ascites.png
│   ├── 05_heatmap_difference.png
│   ├── 06_correlation_scatter.png
│   └── 07_chr20_boxplot.png
└── tables/
    ├── chr_comparison.csv
    ├── chr_paired_tests.csv
    ├── chr_paired_tests_sensitivity.csv
    ├── chr_alteration_summary.csv
    ├── gene_level_cnv_results.csv
    └── patient_chr_data.csv
```

## How to run

```r
# from the project root directory
source("Basel_ovarian_cancer_project.R")
```

Requires R >= 4.2 with: `SingleCellExperiment`, `infercnv`, `dplyr`, `tidyr`, `ggplot2`, `pheatmap`, `RColorBrewer`.

please note : The inferCNV step takes roughly 4-5 hours on 4 threads. If the output already exists it will be skipped automatically.

## Limitations

- 7 of 19 paired patients had fewer than 50 ascites malignant cells, making their per-patient CNV estimates noisy. Sensitivity analysis confirmed the main finding (chr20) is robust to their exclusion.
- inferCNV infers copy number from expression data, it is a proxy, not direct DNA-level measurement.
- Gain/loss thresholds (1.05/0.95) applied to mean CNV across cells are a reasonable but somewhat arbitrary choice for aggregated data.
