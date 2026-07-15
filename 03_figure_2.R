# ==============================================================================
# Figure 2: representative pulmonary remodeling modules
# ==============================================================================

# Packages ---------------------------------------------------------------------

library(circlize)
library(ComplexHeatmap)
library(dplyr)
library(grid)
library(openxlsx)
library(proDA)
library(tibble)
library(tidyr)

# Paths ------------------------------------------------------------------------

project_dir <- "."

input_file <- file.path(
  project_dir,
  "results",
  "figure2",
  "figure2_lung_heatmap_input.rds"
)
figure_dir <- file.path(project_dir, "results", "figure2")
table_dir <- file.path(project_dir, "results", "tables")

output_png <- file.path(figure_dir, "Figure2.png")
output_pdf <- file.path(figure_dir, "Figure2.pdf")
output_table <- file.path(table_dir, "figure2_zscores.xlsx")

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

# Parameters -------------------------------------------------------------------

condition_labels <- c(
  CTRL = "CTRL",
  NOISE = "NOISE",
  NIST = "PM",
  NN = "PM + NOISE"
)
condition_order <- c("CTRL", "NOISE", "PM", "PM + NOISE")

condition_colors <- c(
  "CTRL" = "#BDBDBD",
  "PM" = "#A67C52",
  "NOISE" = "#6C8E7B",
  "PM + NOISE" = "#8C7AA9"
)

heatmap_color_function <- circlize::colorRamp2(
  c(-2, 0, 2),
  c("#4f6383", "lightgrey", "#8a424a")
)

font_family <- "sans"
title_size <- 10
axis_text_size <- 8
legend_text_size <- 7

millimeters_to_inches <- function(x) {
  x / 25.4
}

# Input ------------------------------------------------------------------------

if (!file.exists(input_file)) {
  stop("Missing required input file: ", input_file)
}

heatmap_input <- readRDS(input_file)

values_lung <- heatmap_input$LU$values
annotation_lung <- heatmap_input$LU$annotation
marker_annotation <- heatmap_input$LU$marker_annotation

# Prepare matrix ---------------------------------------------------------------

sample_annotation <- values_lung %>%
  select(Sample, Condition) %>%
  mutate(
    Condition = recode(Condition, !!!condition_labels),
    Condition = factor(Condition, levels = condition_order)
  ) %>%
  arrange(Condition, Sample)

gene_order <- marker_annotation$GeneSymbol
module_order <- unique(marker_annotation$Module)

abundance_long <- values_lung %>%
  select(Sample, all_of(unique(annotation_lung$protein_col))) %>%
  pivot_longer(
    cols = -Sample,
    names_to = "protein_col",
    values_to = "abundance"
  ) %>%
  left_join(
    annotation_lung %>% select(protein_col, GeneSymbol),
    by = "protein_col"
  ) %>%
  filter(
    !is.na(GeneSymbol),
    GeneSymbol %in% gene_order
  ) %>%
  group_by(Sample, GeneSymbol) %>%
  summarise(
    abundance = if (all(is.na(abundance))) {
      NA_real_
    } else {
      median(abundance, na.rm = TRUE)
    },
    .groups = "drop"
  )

abundance_matrix <- abundance_long %>%
  pivot_wider(names_from = Sample, values_from = abundance) %>%
  mutate(GeneSymbol = factor(GeneSymbol, levels = gene_order)) %>%
  arrange(GeneSymbol) %>%
  column_to_rownames("GeneSymbol") %>%
  as.matrix()

abundance_matrix <- abundance_matrix[
  ,
  sample_annotation$Sample,
  drop = FALSE
]

normalized_matrix <- proDA::median_normalization(log2(abundance_matrix))
zscore_matrix <- t(scale(t(normalized_matrix)))

keep_rows <- rowSums(!is.na(zscore_matrix)) >= 2 &
  apply(normalized_matrix, 1, sd, na.rm = TRUE) > 0

zscore_matrix <- zscore_matrix[keep_rows, , drop = FALSE]

column_information <- marker_annotation %>%
  filter(GeneSymbol %in% rownames(zscore_matrix)) %>%
  mutate(
    GeneSymbol = factor(GeneSymbol, levels = gene_order),
    Module = factor(Module, levels = module_order)
  ) %>%
  arrange(Module, GeneSymbol)

zscore_matrix <- zscore_matrix[
  as.character(column_information$GeneSymbol),
  ,
  drop = FALSE
]

# Order genes within each module without displaying dendrograms.
gene_order_within_module <- unlist(
  lapply(
    split(seq_len(nrow(zscore_matrix)), column_information$Module),
    function(index) {
      submatrix <- zscore_matrix[index, , drop = FALSE]
      submatrix[is.na(submatrix)] <- 0
      index[hclust(dist(submatrix), method = "complete")$order]
    }
  )
)

zscore_matrix <- zscore_matrix[
  gene_order_within_module,
  ,
  drop = FALSE
]
column_information <- column_information[
  gene_order_within_module,
  ,
  drop = FALSE
]

plot_matrix <- t(zscore_matrix)

sample_annotation <- sample_annotation %>%
  mutate(Sample = factor(Sample, levels = rownames(plot_matrix))) %>%
  arrange(Sample)

column_information <- column_information %>%
  mutate(
    GeneSymbol = factor(
      as.character(GeneSymbol),
      levels = colnames(plot_matrix)
    )
  ) %>%
  arrange(GeneSymbol)

# Heatmap ----------------------------------------------------------------------

condition_annotation <- rowAnnotation(
  Condition = sample_annotation$Condition,
  col = list(Condition = condition_colors),
  show_annotation_name = FALSE,
  show_legend = FALSE,
  width = unit(3, "mm")
)

heatmap <- Heatmap(
  plot_matrix,
  name = "Z-score",
  col = heatmap_color_function,
  na_col = "white",
  left_annotation = condition_annotation,
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
  row_split = factor(
    sample_annotation$Condition,
    levels = condition_order
  ),
  cluster_row_slices = FALSE,
  row_title_side = "left",
  row_title_rot = 0,
  row_title_gp = gpar(
    fontsize = axis_text_size,
    fontface = "plain",
    fontfamily = font_family
  ),
  column_split = factor(
    column_information$Module,
    levels = module_order
  ),
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

# Export -----------------------------------------------------------------------

openxlsx::write.xlsx(
  as.data.frame(plot_matrix) %>%
    rownames_to_column("Sample"),
  output_table,
  overwrite = TRUE
)

png(
  filename = output_png,
  width = millimeters_to_inches(190),
  height = millimeters_to_inches(75),
  units = "in",
  res = 600,
  bg = "white"
)

draw(
  heatmap,
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
  file = output_pdf,
  width = millimeters_to_inches(190),
  height = millimeters_to_inches(75),
  bg = "white"
)

draw(
  heatmap,
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

message("Figure 2 completed.")
