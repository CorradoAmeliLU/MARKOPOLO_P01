# ==============================================================================
# P01 Fig2
# ==============================================================================
# Manuscript project: Acute co-exposure to particulate matter and aircraft noise in the lung-brain-heart axis
# Target journal: Environmental Pollution
# Author / analyst: Corrado Ameli
# Date: 2026-05-22
#
# Purpose:
#   Generate Figure 2 lung heatmap of representative pulmonary remodeling modules.
#
# Workflow:
#   1. Read the Figure 2 lung matrix exported by the main analysis script.
#   2. Aggregate protein groups to gene-level representative abundance.
#   3. Median-normalize, z-score by gene, order samples/markers, and export PNG/PDF heatmaps.
#
# ==============================================================================

# ------------------------------------------------------------------------------
# PACKAGES
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(proDA)
  library(tibble)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

# ------------------------------------------------------------------------------
# SETTINGS
# ------------------------------------------------------------------------------
obj = readRDS("./Current/P01/Results/Fig2/P_spread_subset_lung_Fig2.rds")
out_dir = "./Current/P01/Results/Fig2/"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
mm_to_in = function(x) x / 25.4
# Match Fig3 font settings
font_family = "sans"
title_size = 10
axis_text_size = 8
legend_text_size = 7
condition_labels = c(
  CTRL = "CTRL",
  NOISE = "NOISE",
  NIST = "PM",
  NN = "PM + NOISE"
)
condition_order = c("CTRL", "NOISE", "PM", "PM + NOISE")
condition_cols = c(
  "CTRL" = "#BDBDBD",
  "PM" = "#A67C52",
  "NOISE" = "#6C8E7B",
  "PM + NOISE" = "#8C7AA9"
)
col_fun = circlize::colorRamp2(
  c(-2.0, 0, 2.0),
  c("#4f6383", "lightgrey", "#8a424a")
)
# ------------------------------------------------------------------------------
# PREPARE DATA
# ------------------------------------------------------------------------------
values_lu = obj$LU$values
annotation_lu = obj$LU$annotation
marker_annotation_lu = obj$LU$marker_annotation %>%
  filter(Module != "Biosynthetic / nuclear regulation")
sample_anno = values_lu %>%
  select(Sample, Condition) %>%
  mutate(
    Condition = recode(Condition, !!!condition_labels),
    Condition = factor(Condition, levels = condition_order)
  ) %>%
  arrange(Condition, Sample)
gene_order = marker_annotation_lu$GeneSymbol
module_order = unique(marker_annotation_lu$Module)
mat_raw = values_lu %>%
  select(Sample, all_of(unique(annotation_lu$protein_col))) %>%
  pivot_longer(
    cols = -Sample,
    names_to = "protein_col",
    values_to = "abundance"
  ) %>%
  left_join(
    annotation_lu %>% select(protein_col, GeneSymbol),
    by = "protein_col"
  ) %>%
  filter(
    !is.na(GeneSymbol),
    GeneSymbol %in% gene_order
  ) %>%
  group_by(Sample, GeneSymbol) %>%
  summarise(
    abundance = ifelse(
      all(is.na(abundance)),
      NA_real_,
      median(abundance, na.rm = TRUE)
    ),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = Sample,
    values_from = abundance
  ) %>%
  mutate(
    GeneSymbol = factor(GeneSymbol, levels = gene_order)
  ) %>%
  arrange(GeneSymbol)
mat = mat_raw %>%
  column_to_rownames("GeneSymbol") %>%
  as.matrix()
mat = mat[, sample_anno$Sample, drop = FALSE]
mat[is.na(mat)] = 1
mat_log2 = log2(mat)
mat_norm = proDA::median_normalization(mat_log2)
mat_z = t(scale(t(mat_norm)))
keep = rowSums(!is.na(mat_z)) >= 2 &
  apply(mat_norm, 1, sd, na.rm = TRUE) > 0
mat_z = mat_z[keep, , drop = FALSE]
row_info = marker_annotation_lu %>%
  filter(GeneSymbol %in% rownames(mat_z)) %>%
  mutate(
    GeneSymbol = factor(GeneSymbol, levels = gene_order),
    Module = factor(Module, levels = module_order)
  ) %>%
  arrange(Module, GeneSymbol)
mat_z = mat_z[as.character(row_info$GeneSymbol), , drop = FALSE]
gene_order_within_module = unlist(
  lapply(
    split(seq_len(nrow(mat_z)), row_info$Module),
    function(idx) {
      submat = mat_z[idx, , drop = FALSE]
      submat_for_order = submat
      submat_for_order[is.na(submat_for_order)] = 0
      idx[hclust(dist(submat_for_order), method = "complete")$order]
    }
  )
)
mat_z = mat_z[gene_order_within_module, , drop = FALSE]
row_info = row_info[gene_order_within_module, , drop = FALSE]
mat_z_t = t(mat_z)
sample_anno = sample_anno %>%
  mutate(
    Sample = factor(Sample, levels = rownames(mat_z_t))
  ) %>%
  arrange(Sample)
col_info = row_info %>%
  mutate(
    GeneSymbol = factor(
      as.character(GeneSymbol),
      levels = colnames(mat_z_t)
    )
  ) %>%
  arrange(GeneSymbol)
# ------------------------------------------------------------------------------
# HEATMAP
# ------------------------------------------------------------------------------
left_ha = rowAnnotation(
  Condition = sample_anno$Condition,
  col = list(Condition = condition_cols),
  show_annotation_name = FALSE,
  show_legend = FALSE,
  width = unit(3, "mm")
)
ht = Heatmap(
  mat_z_t,
  name = "Z-score",
  col = col_fun,
  na_col = "white",
  left_annotation = left_ha,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_dend = FALSE,
  show_column_dend = FALSE,
  show_row_names = FALSE,
  show_column_names = TRUE,
  column_names_side = "bottom",
  column_names_rot = 45,
  column_names_gp = gpar(
    fontsize = axis_text_size,
    fontfamily = font_family
  ),
  row_split = factor(sample_anno$Condition, levels = condition_order),
  cluster_row_slices = FALSE,
  row_title_side = "left",
  row_title_rot = 0,
  row_title_gp = gpar(
    fontsize = axis_text_size,
    fontface = "plain",
    fontfamily = font_family
  ),
  column_split = factor(col_info$Module, levels = module_order),
  cluster_column_slices = FALSE,
  column_gap = unit(1.5, "mm"),
  column_title = NULL,
  show_heatmap_legend = TRUE,
  heatmap_legend_param = list(
    title = "Z-score",
    at = c(-2, 0, 2),
    title_gp = gpar(
      fontsize = legend_text_size,
      fontfamily = font_family
    ),
    labels_gp = gpar(
      fontsize = legend_text_size,
      fontfamily = font_family
    )
  ),
  border = FALSE
)
# =========================
# EXPORTS
# =========================
png(
  filename = file.path(out_dir, "Fig2.png"),
  width = mm_to_in(190),
  height = mm_to_in(75),
  units = "in",
  res = 600,
  bg = "white"
)
draw(
  ht,
  column_title = "Representative pulmonary remodeling modules",
  column_title_gp = gpar(
    fontsize = title_size,
    fontface = "bold",
    fontfamily = font_family
  ),
  heatmap_legend_side = "right",
  annotation_legend_side = "right"
)
dev.off()
pdf(
  file = file.path(out_dir, "Fig2.pdf"),
  width = mm_to_in(190),
  height = mm_to_in(75),
  bg = "white"
)
draw(
  ht,
  column_title = "Representative pulmonary remodeling modules",
  column_title_gp = gpar(
    fontsize = title_size,
    fontface = "bold",
    fontfamily = font_family
  ),
  heatmap_legend_side = "right",
  annotation_legend_side = "right"
)
dev.off()
