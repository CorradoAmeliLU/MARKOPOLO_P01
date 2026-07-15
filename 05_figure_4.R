# ==============================================================================
# Figure 4: transcriptome-proteome concordance in lung
# ==============================================================================

# Packages ---------------------------------------------------------------------

library(dplyr)
library(ggplot2)
library(ggrepel)
library(openxlsx)
library(ragg)
library(stringr)
library(tidyr)

# Paths ------------------------------------------------------------------------

project_dir <- "."

table_dir <- file.path(project_dir, "results", "tables")
figure_dir <- file.path(project_dir, "results", "figure4")
transcriptomics_dir <- file.path(project_dir, "data", "processed")

proteomics_file <- file.path(
  table_dir,
  "proteomics_differential_abundance.xlsx"
)
transcriptomics_file <- file.path(
  transcriptomics_dir,
  "transcriptomics_differential_expression.xlsx"
)

output_png <- file.path(figure_dir, "Figure4.png")
output_pdf <- file.path(figure_dir, "Figure4.pdf")
integration_file <- file.path(
  table_dir,
  "transcriptome_proteome_integration.xlsx"
)

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

# Parameters -------------------------------------------------------------------

conditions <- c("NIST", "NOISE", "NN")
fold_change_threshold <- 1
p_value_threshold <- 0.05

labels_per_condition <- 8
noise_nominal_labels <- 6

condition_labels <- c(
  NIST = "PM",
  NOISE = "NOISE",
  NN = "PM + NOISE"
)
condition_order <- c("PM", "NOISE", "PM + NOISE")

class_colors <- c(
  "Other" = "grey82",
  "Nominally significant" = "#B88746",
  "Concordant significant" = "#8a424a"
)
legend_order <- c(
  "Nominally significant",
  "Concordant significant",
  "Other"
)

font_family <- "sans"
title_size <- 10
axis_title_size <- 8
axis_text_size <- 8
legend_text_size <- 7
strip_text_size <- 8
dpi_export <- 600

# Helper functions -------------------------------------------------------------

check_file_exists <- function(path) {
  if (!file.exists(path)) {
    stop("Missing required input file: ", path)
  }
}

check_columns <- function(data, required_columns, object_name) {
  missing_columns <- setdiff(required_columns, names(data))

  if (length(missing_columns) > 0) {
    stop(
      object_name,
      " is missing required column(s): ",
      paste(missing_columns, collapse = ", ")
    )
  }
}

standardize_region <- function(region) {
  case_when(
    region %in% c("LU", "Lu", "Lung", "Lungs") ~ "LU",
    region %in% c("BR", "Br", "Brain") ~ "BR",
    region %in% c("HR", "Hr", "Heart") ~ "HR",
    TRUE ~ toupper(as.character(region))
  )
}

standardize_comparison <- function(comparison) {
  case_when(
    comparison %in% c("Ni", "NIST", "PM") ~ "NIST",
    comparison %in% c("No", "NOISE") ~ "NOISE",
    comparison %in% c(
      "NiNo",
      "NN",
      "NISTNOISE",
      "PMNOISE",
      "PM_NOISE"
    ) ~ "NN",
    TRUE ~ as.character(comparison)
  )
}

standardize_reference <- function(reference) {
  case_when(
    reference %in% c("Ctr", "CTRL", "Control") ~ "CTRL",
    TRUE ~ as.character(reference)
  )
}

millimeters_to_inches <- function(x) {
  x / 25.4
}

# Input ------------------------------------------------------------------------

check_file_exists(proteomics_file)
check_file_exists(transcriptomics_file)

proteomics <- openxlsx::read.xlsx(proteomics_file)
transcriptomics <- openxlsx::read.xlsx(transcriptomics_file)

check_columns(
  proteomics,
  c("Region", "ident.1", "ident.2", "GeneSymbol", "diff", "pval", "adj_pval"),
  "proteomics"
)
check_columns(
  transcriptomics,
  c("Region", "ident1", "ident2", "gene_symbol", "logFC", "P.Value", "adj.P.Val"),
  "transcriptomics"
)

# Harmonize omics layers -------------------------------------------------------

proteomics_clean <- proteomics %>%
  filter(
    ident.2 == "CTRL",
    ident.1 %in% conditions
  ) %>%
  transmute(
    Region = standardize_region(Region),
    Comparison = ident.1,
    GeneSymbol = str_trim(GeneSymbol),
    protein_logFC = diff,
    protein_pval = pval,
    protein_adj_pval = adj_pval
  ) %>%
  separate_rows(GeneSymbol, sep = ";") %>%
  mutate(GeneSymbol = str_trim(GeneSymbol)) %>%
  filter(!is.na(GeneSymbol), GeneSymbol != "") %>%
  group_by(Region, Comparison, GeneSymbol) %>%
  slice_min(order_by = protein_pval, n = 1, with_ties = FALSE) %>%
  ungroup()

transcriptomics_clean <- transcriptomics %>%
  mutate(
    Region = standardize_region(Region),
    Comparison = standardize_comparison(ident1),
    Reference = standardize_reference(ident2)
  ) %>%
  filter(
    Reference == "CTRL",
    Comparison %in% conditions
  ) %>%
  transmute(
    Region,
    Comparison,
    GeneSymbol = str_trim(gene_symbol),
    transcript_logFC = logFC,
    transcript_pval = P.Value,
    transcript_adj_pval = adj.P.Val
  ) %>%
  filter(!is.na(GeneSymbol), GeneSymbol != "") %>%
  group_by(Region, Comparison, GeneSymbol) %>%
  slice_min(order_by = transcript_pval, n = 1, with_ties = FALSE) %>%
  ungroup()

integration <- inner_join(
  proteomics_clean,
  transcriptomics_clean,
  by = c("Region", "Comparison", "GeneSymbol")
) %>%
  filter(
    Comparison %in% conditions,
    is.finite(protein_logFC),
    is.finite(transcript_logFC)
  ) %>%
  mutate(
    Comparison_label = factor(
      condition_labels[Comparison],
      levels = condition_order
    ),
    Regulation_class = case_when(
      protein_pval < p_value_threshold &
        transcript_pval < p_value_threshold &
        abs(protein_logFC) >= fold_change_threshold &
        abs(transcript_logFC) >= fold_change_threshold &
        sign(protein_logFC) == sign(transcript_logFC) ~
        "Concordant significant",
      protein_pval < p_value_threshold &
        transcript_pval < p_value_threshold ~
        "Nominally significant",
      TRUE ~ "Other"
    ),
    Regulation_class = factor(
      Regulation_class,
      levels = c(
        "Other",
        "Nominally significant",
        "Concordant significant"
      )
    ),
    label_score = abs(protein_logFC) + abs(transcript_logFC)
  )

openxlsx::write.xlsx(
  integration,
  integration_file,
  overwrite = TRUE
)

# Lung figure data -------------------------------------------------------------

lung_data <- integration %>%
  filter(Region == "LU")

main_labels <- lung_data %>%
  filter(Regulation_class == "Concordant significant") %>%
  group_by(Comparison_label) %>%
  slice_max(
    order_by = label_score,
    n = labels_per_condition,
    with_ties = FALSE
  ) %>%
  ungroup()

noise_labels <- lung_data %>%
  filter(
    Comparison_label == "NOISE",
    Regulation_class == "Nominally significant"
  ) %>%
  slice_max(
    order_by = label_score,
    n = noise_nominal_labels,
    with_ties = FALSE
  )

label_data <- bind_rows(main_labels, noise_labels) %>%
  distinct(Comparison_label, GeneSymbol, .keep_all = TRUE)

correlation_labels <- lung_data %>%
  group_by(Comparison_label) %>%
  summarise(
    correlation = suppressWarnings(
      cor(
        transcript_logFC,
        protein_logFC,
        method = "pearson"
      )
    ),
    .groups = "drop"
  ) %>%
  mutate(
    label = paste0("r = ", round(correlation, 2)),
    x = -Inf,
    y = Inf
  )

# Plot -------------------------------------------------------------------------

figure_4 <- ggplot(
  lung_data,
  aes(x = transcript_logFC, y = protein_logFC)
) +
  geom_hline(
    yintercept = 0,
    color = "grey70",
    linewidth = 0.35
  ) +
  geom_vline(
    xintercept = 0,
    color = "grey70",
    linewidth = 0.35
  ) +
  geom_abline(
    slope = 1,
    intercept = 0,
    color = "grey55",
    linewidth = 0.45,
    linetype = "dashed"
  ) +
  geom_point(
    data = lung_data %>% filter(Regulation_class == "Other"),
    aes(color = Regulation_class),
    size = 1.4,
    alpha = 0.55
  ) +
  geom_point(
    data = lung_data %>%
      filter(Regulation_class == "Nominally significant"),
    aes(color = Regulation_class),
    size = 1.6,
    alpha = 0.80
  ) +
  geom_point(
    data = lung_data %>%
      filter(Regulation_class == "Concordant significant"),
    aes(color = Regulation_class),
    size = 1.9,
    alpha = 0.95
  ) +
  geom_text_repel(
    data = label_data,
    aes(label = GeneSymbol),
    size = 2.5,
    color = "black",
    max.overlaps = Inf,
    min.segment.length = 0,
    segment.color = "grey60",
    segment.linewidth = 0.25,
    box.padding = 0.30,
    point.padding = 0.15,
    show.legend = FALSE,
    seed = 123
  ) +
  geom_text(
    data = correlation_labels,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust = -0.15,
    vjust = 1.15,
    size = 2.7,
    color = "black"
  ) +
  facet_wrap(~ Comparison_label, nrow = 1) +
  scale_color_manual(
    values = class_colors,
    breaks = legend_order,
    drop = FALSE
  ) +
  labs(
    title = "Transcriptome-proteome concordance in lung",
    x = expression(Transcript~log[2]~fold~change),
    y = expression(Protein~log[2]~fold~change)
  ) +
  theme_classic(base_size = 9, base_family = font_family) +
  theme(
    text = element_text(family = font_family, color = "black"),
    plot.title = element_text(
      face = "bold",
      size = title_size,
      hjust = 0.5
    ),
    axis.title = element_text(size = axis_title_size),
    axis.text = element_text(size = axis_text_size),
    strip.background = element_blank(),
    strip.text = element_text(
      face = "bold",
      size = strip_text_size
    ),
    legend.title = element_blank(),
    legend.text = element_text(size = legend_text_size),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box = "horizontal",
    plot.margin = margin(5, 5, 5, 5)
  ) +
  guides(
    color = guide_legend(
      nrow = 1,
      byrow = TRUE,
      override.aes = list(size = 2.4, alpha = 1)
    )
  )

# Export -----------------------------------------------------------------------

ggsave(
  output_png,
  plot = figure_4,
  width = millimeters_to_inches(180),
  height = millimeters_to_inches(85),
  units = "in",
  dpi = dpi_export,
  device = ragg::agg_png,
  bg = "white"
)

ggsave(
  output_pdf,
  plot = figure_4,
  width = millimeters_to_inches(180),
  height = millimeters_to_inches(85),
  units = "in",
  device = cairo_pdf,
  bg = "white"
)

message("Figure 4 completed.")
