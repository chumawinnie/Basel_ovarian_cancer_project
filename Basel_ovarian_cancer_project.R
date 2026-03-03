#!/usr/bin/env Rscript

# ============================================================================
# inferCNV Analysis: Solid vs Ascites in HGSOC (scRNA-seq)
# 
# Author: Chukwuma Winner Obiora
# Date: March 2026
# 
# Dataset: Vasquez-Garcia et al. (Nature, 2022) - SPECTRUM cohort
# 5,000 genes x 70,437 cells from 25 HGSOC patients
#
# Task: Compare CNV profiles between solid tumor and ascites compartments
#       using inferCNV on single-cell RNA-seq data.
# ============================================================================

library(SingleCellExperiment)
library(infercnv)
library(dplyr)
library(tidyr)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)

# -- paths
data_dir <- "Data"
output_dir <- "Infer_copy_number_results"
dir.create(output_dir, showWarnings = FALSE)
dir.create(file.path(output_dir, "figures"), showWarnings = FALSE)
dir.create(file.path(output_dir, "tables"), showWarnings = FALSE)


# ============================================================================
# 1. Data loading and exploration
# ============================================================================

sce <- readRDS(file.path(data_dir, "sce_subset.RDS"))
gene_pos <- read.delim(file.path(data_dir, "gene_positions.tsv"),
                        header = FALSE, stringsAsFactors = FALSE)
colnames(gene_pos) <- c("gene", "chr", "start", "stop")

cat("SCE dimensions:", dim(sce), "\n")
cat("Gene positions:", nrow(gene_pos), "\n")

meta <- as.data.frame(colData(sce))

cat("\nCell type breakdown:\n")
print(table(meta$cell_type))
cat("\nTumor site breakdown:\n")
print(table(meta$tumor_site))

# which patients have malignant cells in both compartments?
mal_meta <- meta[meta$cell_type == "Malignant", ]
patient_site <- mal_meta %>%
  group_by(patient_id, tumor_site) %>%
  summarise(n = n(), .groups = "drop") %>%
  pivot_wider(names_from = tumor_site, values_from = n, values_fill = 0)

patients_both <- patient_site$patient_id[patient_site$Solid > 0 & patient_site$Ascites > 0]
cat("\nPatients with both Solid + Ascites malignant cells:", length(patients_both), "\n")

# flag patients with very few ascites cells -- revisited in sensitivity analysis
low_ascites <- patient_site %>%
  filter(patient_id %in% patients_both, Ascites < 50) %>%
  select(patient_id, Solid, Ascites)
cat("\nPatients with <50 Ascites malignant cells:\n")
print(as.data.frame(low_ascites))


# ============================================================================
# 2. Prepare inferCNV inputs
# ============================================================================

autosome_chrs <- paste0("chr", 1:22)
gene_pos_auto <- gene_pos[gene_pos$chr %in% autosome_chrs, ]

overlap_genes <- intersect(rownames(sce), gene_pos_auto$gene)
cat("\nGene overlap (SCE and autosomal positions):", length(overlap_genes), "\n")

gene_order_df <- gene_pos_auto[gene_pos_auto$gene %in% overlap_genes, ]
rownames(gene_order_df) <- gene_order_df$gene
gene_order_df <- gene_order_df[, c("chr", "start", "stop")]

# group malignant cells by patient + site, pool all reference into one group
annotations <- data.frame(
  cell_id = meta$cell_id,
  group = ifelse(meta$cell_type == "Malignant",
                 paste0(meta$patient_id, "_", meta$tumor_site),
                 "Reference"),
  stringsAsFactors = FALSE
)
rownames(annotations) <- annotations$cell_id

ann_file <- file.path(output_dir, "cell_annotations.txt")
write.table(annotations[, "group", drop = FALSE], ann_file,
            sep = "\t", quote = FALSE, col.names = FALSE)

counts <- as.matrix(counts(sce))


# ============================================================================
# 3. Run inferCNV
# ============================================================================
# Takes ~4-5 hours on 4 threads. Skips if output already exists.

run_infercnv <- !file.exists(file.path(output_dir, "infercnv_run", "run.final.infercnv_obj"))

if (run_infercnv) {
  
  infercnv_obj <- CreateInfercnvObject(
    raw_counts_matrix = counts,
    annotations_file = ann_file,
    gene_order_file = gene_order_df,
    ref_group_names = c("Reference")
  )
  
  infercnv_obj <- infercnv::run(
    infercnv_obj,
    cutoff = 0.1,
    out_dir = file.path(output_dir, "infercnv_run"),
    cluster_by_groups = TRUE,
    denoise = TRUE,
    HMM = FALSE,
    diagnostics = FALSE,
    plot_steps = FALSE,
    inspect_subclusters = FALSE,
    no_prelim_plot = TRUE,
    num_threads = 4
  )
  
  cat("inferCNV run complete.\n")
  
} else {
  cat("inferCNV output found, skipping run.\n")
}


# ============================================================================
# 4. Extract CNV results
# ============================================================================

infercnv_final <- readRDS(file.path(output_dir, "infercnv_run", "run.final.infercnv_obj"))
cnv_matrix <- infercnv_final@expr.data

cat("\nCNV matrix:", dim(cnv_matrix)[1], "genes x", dim(cnv_matrix)[2], "cells\n")
cat("Value range:", round(range(cnv_matrix), 4), "\n")

gene_order <- infercnv_final@gene_order
gene_chr <- data.frame(
  gene = rownames(gene_order),
  chr = gene_order[, 1],
  start = gene_order[, 2],
  stop = gene_order[, 3],
  stringsAsFactors = FALSE
)
gene_chr <- gene_chr[gene_chr$chr %in% autosome_chrs, ]

# malignant cells only
malignant_cells <- meta$cell_id[meta$cell_type == "Malignant"]
malignant_in_cnv <- intersect(malignant_cells, colnames(cnv_matrix))
genes_to_use <- intersect(rownames(cnv_matrix), gene_chr$gene)

cnv_mal <- cnv_matrix[genes_to_use, malignant_in_cnv]
cat("Working matrix:", dim(cnv_mal)[1], "genes x", dim(cnv_mal)[2], "malignant cells\n")

meta_mal <- meta[meta$cell_id %in% malignant_in_cnv, ]
rownames(meta_mal) <- meta_mal$cell_id
solid_cells <- meta_mal$cell_id[meta_mal$tumor_site == "Solid"]
ascites_cells <- meta_mal$cell_id[meta_mal$tumor_site == "Ascites"]

cat("Solid:", length(solid_cells), "| Ascites:", length(ascites_cells), "\n")


# ============================================================================
# 5. Per-chromosome CNV scores
# ============================================================================

gene_chr_map <- gene_chr[gene_chr$gene %in% genes_to_use, ]
rownames(gene_chr_map) <- gene_chr_map$gene

chr_cnv_per_cell <- matrix(NA, nrow = 22, ncol = ncol(cnv_mal))
rownames(chr_cnv_per_cell) <- autosome_chrs
colnames(chr_cnv_per_cell) <- colnames(cnv_mal)

for (ch in autosome_chrs) {
  g <- gene_chr_map$gene[gene_chr_map$chr == ch]
  g <- intersect(g, rownames(cnv_mal))
  if (length(g) > 1) {
    chr_cnv_per_cell[ch, ] <- colMeans(cnv_mal[g, , drop = FALSE])
  } else if (length(g) == 1) {
    chr_cnv_per_cell[ch, ] <- cnv_mal[g, ]
  }
}

mean_solid <- rowMeans(chr_cnv_per_cell[, solid_cells, drop = FALSE])
mean_ascites <- rowMeans(chr_cnv_per_cell[, ascites_cells, drop = FALSE])

chr_comparison <- data.frame(
  chromosome = autosome_chrs,
  chr_num = 1:22,
  mean_solid = round(mean_solid, 6),
  mean_ascites = round(mean_ascites, 6),
  difference = round(mean_ascites - mean_solid, 6),
  stringsAsFactors = FALSE
)

overall_cor <- cor(mean_solid, mean_ascites, method = "pearson")
cat("\nPearson r (Solid vs Ascites):", round(overall_cor, 4), "\n")

write.csv(chr_comparison, file.path(output_dir, "tables", "chr_comparison.csv"), row.names = FALSE)


# ============================================================================
# 6. Per-patient paired analysis
# ============================================================================

patient_chr_list <- list()

for (p in patients_both) {
  for (site in c("Solid", "Ascites")) {
    cells <- meta_mal$cell_id[meta_mal$patient_id == p & meta_mal$tumor_site == site]
    cells <- intersect(cells, colnames(chr_cnv_per_cell))
    n_cells <- length(cells)
    
    if (n_cells == 0) next
    
    if (n_cells == 1) {
      means <- chr_cnv_per_cell[, cells]
    } else {
      means <- rowMeans(chr_cnv_per_cell[, cells, drop = FALSE])
    }
    
    patient_chr_list[[paste0(p, "_", site)]] <- data.frame(
      patient_id = p, tumor_site = site, n_cells = n_cells,
      chromosome = names(means),
      chr_num = as.numeric(gsub("chr", "", names(means))),
      mean_cnv = means,
      stringsAsFactors = FALSE, row.names = NULL
    )
  }
}

patient_chr_df <- bind_rows(patient_chr_list)

# reusable function for paired tests
run_paired_tests <- function(patients_to_use, patient_data) {
  results <- data.frame()
  for (ch in autosome_chrs) {
    solid_vals <- c(); ascites_vals <- c()
    for (p in patients_to_use) {
      s <- patient_data[patient_data$patient_id == p &
                          patient_data$tumor_site == "Solid" &
                          patient_data$chromosome == ch, "mean_cnv"]
      a <- patient_data[patient_data$patient_id == p &
                          patient_data$tumor_site == "Ascites" &
                          patient_data$chromosome == ch, "mean_cnv"]
      if (length(s) == 1 & length(a) == 1) {
        solid_vals <- c(solid_vals, s)
        ascites_vals <- c(ascites_vals, a)
      }
    }
    if (length(solid_vals) >= 3) {
      test <- wilcox.test(ascites_vals, solid_vals, paired = TRUE)
      results <- rbind(results, data.frame(
        chromosome = ch,
        chr_num = as.numeric(gsub("chr", "", ch)),
        n_patients = length(solid_vals),
        mean_solid = round(mean(solid_vals), 6),
        mean_ascites = round(mean(ascites_vals), 6),
        mean_diff = round(mean(ascites_vals - solid_vals), 6),
        p_value = round(test$p.value, 6),
        stringsAsFactors = FALSE
      ))
    }
  }
  results$p_adj <- round(p.adjust(results$p_value, method = "BH"), 6)
  results[order(results$chr_num), ]
}

# main test with all 19 patients
chr_tests <- run_paired_tests(patients_both, patient_chr_df)

cat("\n--- Paired Wilcoxon Results (all 19 patients) ---\n")
print(chr_tests, row.names = FALSE)

write.csv(chr_tests, file.path(output_dir, "tables", "chr_paired_tests.csv"), row.names = FALSE)
write.csv(patient_chr_df, file.path(output_dir, "tables", "patient_chr_data.csv"), row.names = FALSE)


# ============================================================================
# 7. Sensitivity analysis
# ============================================================================
# 7 patients had <50 ascites malignant cells (some as few as 1-3 cells).
# inferCNV estimates from so few cells are unreliable, so we check whether
# our findings hold after excluding them.

ascites_counts <- meta_mal %>%
  filter(tumor_site == "Ascites") %>%
  group_by(patient_id) %>%
  summarise(n_ascites = n(), .groups = "drop")

reliable_patients <- ascites_counts$patient_id[ascites_counts$n_ascites >= 50]
excluded <- setdiff(patients_both, reliable_patients)

cat("\n--- Sensitivity Analysis ---\n")
cat("Excluded (< 50 ascites cells):", paste(excluded, collapse = ", "), "\n")
cat("Remaining:", length(reliable_patients), "patients\n")

chr_tests_sens <- run_paired_tests(reliable_patients, patient_chr_df)

# chr20 is the main finding -- check if it holds
main_chr20 <- chr_tests[chr_tests$chromosome == "chr20", ]
sens_chr20 <- chr_tests_sens[chr_tests_sens$chromosome == "chr20", ]
cat(sprintf("\nchr20: Main p_adj=%.4f (n=%d) -> Sensitivity p_adj=%.4f (n=%d)\n",
            main_chr20$p_adj, main_chr20$n_patients,
            sens_chr20$p_adj, sens_chr20$n_patients))

write.csv(chr_tests_sens, file.path(output_dir, "tables", "chr_paired_tests_sensitivity.csv"),
          row.names = FALSE)


# ============================================================================
# 8. Identify frequently altered regions
# ============================================================================

gain_threshold <- 1.05
loss_threshold <- 0.95

gene_mean_cnv <- rowMeans(cnv_mal)
gene_mean_solid <- rowMeans(cnv_mal[, solid_cells, drop = FALSE])
gene_mean_ascites <- rowMeans(cnv_mal[, ascites_cells, drop = FALSE])

gene_results <- data.frame(
  gene = names(gene_mean_cnv),
  mean_cnv = gene_mean_cnv,
  chr = gene_chr_map[names(gene_mean_cnv), "chr"],
  start = gene_chr_map[names(gene_mean_cnv), "start"],
  stop = gene_chr_map[names(gene_mean_cnv), "stop"],
  mean_solid = gene_mean_solid[names(gene_mean_cnv)],
  mean_ascites = gene_mean_ascites[names(gene_mean_cnv)],
  stringsAsFactors = FALSE
)

gene_results$alteration <- ifelse(gene_results$mean_cnv > gain_threshold, "Gain",
                                  ifelse(gene_results$mean_cnv < loss_threshold, "Loss", "Neutral"))

gene_results$solid_status <- ifelse(gene_results$mean_solid > gain_threshold, "Gain",
                                    ifelse(gene_results$mean_solid < loss_threshold, "Loss", "Neutral"))
gene_results$ascites_status <- ifelse(gene_results$mean_ascites > gain_threshold, "Gain",
                                      ifelse(gene_results$mean_ascites < loss_threshold, "Loss", "Neutral"))

gene_results$shared <- (gene_results$solid_status == gene_results$ascites_status) &
  (gene_results$solid_status != "Neutral")
gene_results$site_specific <- (gene_results$solid_status != gene_results$ascites_status) &
  (gene_results$solid_status != "Neutral" | gene_results$ascites_status != "Neutral")

cat("\n--- Gene-level alteration counts ---\n")
cat("Total gains:", sum(gene_results$alteration == "Gain"), "\n")
cat("Total losses:", sum(gene_results$alteration == "Loss"), "\n")
cat("Shared gains:", sum(gene_results$shared & gene_results$solid_status == "Gain"), "\n")
cat("Shared losses:", sum(gene_results$shared & gene_results$solid_status == "Loss"), "\n")
cat("Site-specific:", sum(gene_results$site_specific), "\n")

# per-chromosome summary
chr_alt_summary <- gene_results %>%
  mutate(chr_num = as.numeric(gsub("chr", "", chr))) %>%
  group_by(chr, chr_num) %>%
  summarise(
    n_genes = n(),
    n_gains = sum(alteration == "Gain"),
    n_losses = sum(alteration == "Loss"),
    pct_gains = round(100 * n_gains / n(), 1),
    pct_losses = round(100 * n_losses / n(), 1),
    mean_cnv = round(mean(mean_cnv), 4),
    n_shared_gains = sum(shared & solid_status == "Gain"),
    n_shared_losses = sum(shared & solid_status == "Loss"),
    n_site_specific = sum(site_specific),
    .groups = "drop"
  ) %>%
  arrange(chr_num)

print(as.data.frame(chr_alt_summary))

# top altered genes
cat("\nTop 15 amplified:\n")
print(head(gene_results[order(-gene_results$mean_cnv),
                         c("gene","chr","mean_cnv","mean_solid","mean_ascites")], 15),
      row.names = FALSE)

cat("\nTop 15 deleted:\n")
print(head(gene_results[order(gene_results$mean_cnv),
                         c("gene","chr","mean_cnv","mean_solid","mean_ascites")], 15),
      row.names = FALSE)

write.csv(gene_results, file.path(output_dir, "tables", "gene_level_cnv_results.csv"), row.names = FALSE)
write.csv(chr_alt_summary, file.path(output_dir, "tables", "chr_alteration_summary.csv"), row.names = FALSE)


# ============================================================================
# 9. Figures
# ============================================================================

# Fig 1: Genome-wide CNV profile

gene_results_ordered <- gene_results %>%
  mutate(chr_num = as.numeric(gsub("chr", "", chr))) %>%
  arrange(chr_num, start)
gene_results_ordered$genome_pos <- seq_len(nrow(gene_results_ordered))

chr_boundaries <- gene_results_ordered %>%
  group_by(chr_num) %>%
  summarise(start_pos = min(genome_pos), end_pos = max(genome_pos),
            mid_pos = mean(genome_pos), .groups = "drop")

p1 <- ggplot(gene_results_ordered) +
  geom_rect(data = chr_boundaries,
            aes(xmin = start_pos, xmax = end_pos, ymin = -Inf, ymax = Inf,
                fill = factor(chr_num %% 2)), alpha = 0.1) +
  scale_fill_manual(values = c("0" = "white", "1" = "grey90"), guide = "none") +
  geom_point(aes(x = genome_pos, y = mean_solid), color = "#2166AC", alpha = 0.4, size = 0.5) +
  geom_point(aes(x = genome_pos, y = mean_ascites), color = "#B2182B", alpha = 0.4, size = 0.5) +
  geom_hline(yintercept = 1, color = "black", linewidth = 0.3) +
  geom_hline(yintercept = c(gain_threshold, loss_threshold),
             linetype = "dashed", color = "grey50", linewidth = 0.3) +
  scale_x_continuous(breaks = chr_boundaries$mid_pos, labels = chr_boundaries$chr_num) +
  labs(x = "Chromosome", y = "inferCNV Score",
       title = "Genome-wide CNV Profile: Solid vs Ascites",
       subtitle = "Blue = Solid | Red = Ascites | Dashed = gain/loss thresholds") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 8),
        plot.title = element_text(face = "bold", size = 14),
        panel.grid.minor = element_blank())

ggsave(file.path(output_dir, "figures", "01_genome_wide_cnv_profile.png"),
       p1, width = 14, height = 5, dpi = 300)


# Fig 2: Per-chromosome bar plot

chr_comp_long <- chr_comparison %>%
  pivot_longer(cols = c(mean_solid, mean_ascites),
               names_to = "site", values_to = "mean_cnv") %>%
  mutate(site = ifelse(site == "mean_solid", "Solid", "Ascites"))

p2 <- ggplot(chr_comp_long, aes(x = factor(chr_num), y = mean_cnv - 1, fill = site)) +
  geom_bar(stat = "identity", position = position_dodge(0.7), width = 0.6) +
  geom_hline(yintercept = 0, linewidth = 0.5) +
  scale_fill_manual(values = c("Solid" = "#2166AC", "Ascites" = "#B2182B")) +
  labs(x = "Chromosome", y = "Mean CNV Deviation from Baseline",
       title = "Per-Chromosome CNV: Solid vs Ascites",
       fill = "Tumor Site") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"), legend.position = "top")

ggsave(file.path(output_dir, "figures", "02_chr_barplot_comparison.png"),
       p2, width = 12, height = 5, dpi = 300)


# Fig 3: Alteration frequency

chr_alt_long <- chr_alt_summary %>%
  select(chr_num, pct_gains, pct_losses) %>%
  pivot_longer(cols = c(pct_gains, pct_losses),
               names_to = "type", values_to = "percent") %>%
  mutate(type = ifelse(type == "pct_gains", "Gains", "Losses"),
         percent = ifelse(type == "Losses", -percent, percent))

p3 <- ggplot(chr_alt_long, aes(x = factor(chr_num), y = percent, fill = type)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_hline(yintercept = 0, linewidth = 0.5) +
  scale_fill_manual(values = c("Gains" = "#D73027", "Losses" = "#4575B4")) +
  labs(x = "Chromosome", y = "% Genes Altered",
       title = "Copy Number Alteration Frequency by Chromosome", fill = "") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"), legend.position = "top")

ggsave(file.path(output_dir, "figures", "03_chr_alteration_frequency.png"),
       p3, width = 12, height = 5, dpi = 300)


# Fig 4a/4b: Per-patient heatmaps

build_patient_chr_matrix <- function(df, site) {
  sub <- df[df$tumor_site == site, ]
  mat <- sub %>%
    select(patient_id, chromosome, mean_cnv) %>%
    pivot_wider(names_from = chromosome, values_from = mean_cnv) %>%
    as.data.frame()
  rownames(mat) <- mat$patient_id
  mat$patient_id <- NULL
  chr_present <- intersect(paste0("chr", 1:22), colnames(mat))
  as.matrix(mat[, chr_present])
}

mat_solid <- build_patient_chr_matrix(patient_chr_df, "Solid")
mat_ascites <- build_patient_chr_matrix(patient_chr_df, "Ascites")

hm_breaks <- seq(0.85, 1.15, length.out = 100)
hm_colors <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(99)

png(file.path(output_dir, "figures", "04a_heatmap_solid.png"),
    width = 10, height = 8, units = "in", res = 300)
pheatmap(mat_solid, cluster_cols = FALSE, color = hm_colors, breaks = hm_breaks,
         main = "Per-Patient CNV: Solid Tumors",
         labels_col = gsub("chr", "", colnames(mat_solid)), fontsize = 10)
dev.off()

png(file.path(output_dir, "figures", "04b_heatmap_ascites.png"),
    width = 10, height = 8, units = "in", res = 300)
pheatmap(mat_ascites, cluster_cols = FALSE, color = hm_colors, breaks = hm_breaks,
         main = "Per-Patient CNV: Ascites",
         labels_col = gsub("chr", "", colnames(mat_ascites)), fontsize = 10)
dev.off()


# Fig 5: Difference heatmap

common_patients <- intersect(rownames(mat_solid), rownames(mat_ascites))
common_chrs <- intersect(colnames(mat_solid), colnames(mat_ascites))
mat_diff <- mat_ascites[common_patients, common_chrs] - mat_solid[common_patients, common_chrs]

png(file.path(output_dir, "figures", "05_heatmap_difference.png"),
    width = 10, height = 8, units = "in", res = 300)
pheatmap(mat_diff, cluster_cols = FALSE,
         color = colorRampPalette(rev(brewer.pal(11, "PiYG")))(99),
         breaks = seq(-0.1, 0.1, length.out = 100),
         main = "CNV Difference (Ascites - Solid)",
         labels_col = gsub("chr", "", colnames(mat_diff)), fontsize = 10)
dev.off()


# Fig 6: Correlation scatter

p6 <- ggplot(chr_comparison, aes(x = mean_solid, y = mean_ascites)) +
  geom_point(size = 3, color = "#333333") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  geom_text(aes(label = chr_num), nudge_x = 0.003, nudge_y = 0.003, size = 3) +
  labs(x = "Mean CNV (Solid)", y = "Mean CNV (Ascites)",
       title = "Solid vs Ascites CNV Correlation",
       subtitle = paste0("Pearson r = ", round(overall_cor, 3))) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"))

ggsave(file.path(output_dir, "figures", "06_correlation_scatter.png"),
       p6, width = 7, height = 6, dpi = 300)


# Fig 7: Chr20 boxplot -- the main site-specific finding that survives
# both the biological threshold check and the sensitivity analysis

chr20_data <- patient_chr_df[patient_chr_df$chromosome == "chr20" &
                               patient_chr_df$patient_id %in% patients_both, ]

p7 <- ggplot(chr20_data, aes(x = "Chr20", y = mean_cnv, fill = tumor_site)) +
  geom_boxplot(position = position_dodge(0.8), width = 0.6, outlier.size = 1) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
  annotate("text", x = 1, y = max(chr20_data$mean_cnv) + 0.02,
           label = "p_adj = 0.045 (n=19)", size = 3.5, fontface = "italic") +
  scale_fill_manual(values = c("Solid" = "#2166AC", "Ascites" = "#B2182B")) +
  labs(x = NULL, y = "Mean CNV Score (per patient)",
       title = "Chr20 Amplification: Solid vs Ascites",
       subtitle = "Paired Wilcoxon, BH-adjusted",
       fill = "Tumor Site") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 13), legend.position = "top")

ggsave(file.path(output_dir, "figures", "07_chr20_boxplot.png"),
       p7, width = 6, height = 6, dpi = 300)


cat("\nAll figures saved to", file.path(output_dir, "figures"), "\n")
cat("All tables saved to", file.path(output_dir, "tables"), "\n")
cat("Done.\n")
