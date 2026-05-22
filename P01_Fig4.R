# ==============================================================================
# P01 Fig4
# ==============================================================================
# Manuscript project: Acute co-exposure to particulate matter and aircraft noise in the lung-brain-heart axis
# Target journal: Environmental Pollution
# Author / analyst: Corrado Ameli
# Date: 2026-05-22
#
# Purpose:
#   Generate Figure 4 transcriptome-proteome fold-change scatter plots for lung contrasts.
#
# Workflow:
#   1. Read proteomics and transcriptomics differential-expression outputs.
#   2. Harmonize region and condition labels across omics layers.
#   3. Classify concordant and nominally significant genes, add selected labels, and export scatter plots.
#
# ==============================================================================

# ------------------------------------------------------------------------------
# PACKAGES
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(openxlsx)
  library(stringr)
  library(ggrepel)
})

# ---- Input data ----

P = read.xlsx("./Current/P01/Results/DiffExpr_ProDA.xlsx")
T = read.xlsx("./PreviousStudies/Results/4d/DEGS_unfiltered.xlsx")
# =========================
# P01 FIGURE 4
# Transcriptome-proteome fold-change scatter
# Lung only
# Requires existing objects: P and T
# =========================
# ------------------------------------------------------------------------------
# SETTINGS
# ------------------------------------------------------------------------------
out_dir = "./Current/P01/Results/Fig4/"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
mm_to_in = function(x) x / 25.4
target_region = "LU"
conditions_to_plot = c("NIST", "NOISE", "NN")
fc_threshold = 1
p_threshold = 0.05
label_top_n_per_facet = 8
label_noise_nominal_top_n = 6
# Match Fig3 font sizes
font_family = "sans"
title_size = 10
axis_title_size = 8
axis_text_size = 8
legend_text_size = 7
strip_text_size = 8
# ------------------------------------------------------------------------------
# HARMONIZE PROTEOMICS
# ------------------------------------------------------------------------------
P_clean = P %>%
  filter(
    ident.2 == "CTRL",
    ident.1 %in% conditions_to_plot
  ) %>%
  transmute(
    Region = case_when(
      Region %in% c("LU", "Lu", "Lung", "Lungs") ~ "LU",
      Region %in% c("BR", "Br", "Brain") ~ "BR",
      Region %in% c("HR", "Hr", "Heart") ~ "HR",
      TRUE ~ toupper(as.character(Region))
    ),
    Comparison = ident.1,
    GeneSymbol = str_trim(GeneSymbol),
    protein_logFC = diff,
    protein_pval = pval,
    protein_adj_pval = adj_pval
  ) %>%
  filter(
    !is.na(GeneSymbol),
    GeneSymbol != ""
  )
# ------------------------------------------------------------------------------
# HARMONIZE TRANSCRIPTOMICS
# ------------------------------------------------------------------------------
T_clean = T %>%
  mutate(
    Region = case_when(
      Region %in% c("LU", "Lu", "Lung", "Lungs") ~ "LU",
      Region %in% c("BR", "Br", "Brain") ~ "BR",
      Region %in% c("HR", "Hr", "Heart") ~ "HR",
      TRUE ~ toupper(as.character(Region))
    ),
    Comparison = case_when(
      ident1 %in% c("Ni", "NIST", "PM") ~ "NIST",
      ident1 %in% c("No", "NOISE") ~ "NOISE",
      ident1 %in% c("NiNo", "NN", "NISTNOISE", "PMNOISE", "PM_NOISE") ~ "NN",
      TRUE ~ as.character(ident1)
    ),
    Reference = case_when(
      ident2 %in% c("Ctr", "CTRL", "Control") ~ "CTRL",
      TRUE ~ as.character(ident2)
    )
  ) %>%
  filter(
    Reference == "CTRL",
    Comparison %in% conditions_to_plot
  ) %>%
  transmute(
    Region,
    Comparison,
    GeneSymbol = str_trim(gene_symbol),
    transcript_logFC = logFC,
    transcript_pval = P.Value,
    transcript_adj_pval = adj.P.Val
  ) %>%
  filter(
    !is.na(GeneSymbol),
    GeneSymbol != ""
  )
# ------------------------------------------------------------------------------
# JOIN OMICS LAYERS
# ------------------------------------------------------------------------------
omics_df = inner_join(
  P_clean,
  T_clean,
  by = c("Region", "Comparison", "GeneSymbol")
) %>%
  filter(
    Region == target_region,
    Comparison %in% conditions_to_plot,
    is.finite(protein_logFC),
    is.finite(transcript_logFC)
  ) %>%
  mutate(
    Comparison = factor(
      Comparison,
      levels = c("NIST", "NOISE", "NN"),
      labels = c("PM", "NOISE", "PM + NOISE")
    ),
    regulation_class = case_when(
      protein_pval < p_threshold &
        transcript_pval < p_threshold &
        abs(protein_logFC) >= fc_threshold &
        abs(transcript_logFC) >= fc_threshold &
        sign(protein_logFC) == sign(transcript_logFC) ~
        "Concordant significant",
      protein_pval < p_threshold &
        transcript_pval < p_threshold ~
        "Nominally significant",
      TRUE ~ "Other"
    ),
    regulation_class = factor(
      regulation_class,
      levels = c(
        "Other",
        "Nominally significant",
        "Concordant significant"
      )
    ),
    label_score = abs(protein_logFC) + abs(transcript_logFC)
  )
# ------------------------------------------------------------------------------
# LABELS
# ------------------------------------------------------------------------------
label_df_main = omics_df %>%
  filter(regulation_class == "Concordant significant") %>%
  group_by(Comparison) %>%
  slice_max(
    order_by = label_score,
    n = label_top_n_per_facet,
    with_ties = FALSE
  ) %>%
  ungroup()
label_df_noise_nominal = omics_df %>%
  filter(
    Comparison == "NOISE",
    regulation_class == "Nominally significant"
  ) %>%
  slice_max(
    order_by = label_score,
    n = label_noise_nominal_top_n,
    with_ties = FALSE
  )
label_df = bind_rows(
  label_df_main,
  label_df_noise_nominal
) %>%
  distinct(Comparison, GeneSymbol, .keep_all = TRUE)
# ------------------------------------------------------------------------------
# CORRELATION LABELS
# ------------------------------------------------------------------------------
cor_df = omics_df %>%
  group_by(Comparison) %>%
  summarise(
    r = suppressWarnings(cor(transcript_logFC, protein_logFC, method = "pearson")),
    .groups = "drop"
  ) %>%
  mutate(
    label = paste0("r = ", round(r, 2)),
    x = -Inf,
    y = Inf
  )
# ------------------------------------------------------------------------------
# COLORS
# ------------------------------------------------------------------------------
class_cols = c(
  "Other" = "grey82",
  "Nominally significant" = "#B88746",
  "Concordant significant" = "#8a424a"
)
legend_order = c(
  "Nominally significant",
  "Concordant significant",
  "Other"
)
# =========================
# PLOT
# =========================
p_fig4_lung = ggplot(
  omics_df,
  aes(
    x = transcript_logFC,
    y = protein_logFC
  )
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
    data = omics_df %>% filter(regulation_class == "Other"),
    aes(color = regulation_class),
    size = 1.4,
    alpha = 0.55
  ) +
  geom_point(
    data = omics_df %>% filter(regulation_class == "Nominally significant"),
    aes(color = regulation_class),
    size = 1.6,
    alpha = 0.80
  ) +
  geom_point(
    data = omics_df %>% filter(regulation_class == "Concordant significant"),
    aes(color = regulation_class),
    size = 1.9,
    alpha = 0.95
  ) +
  geom_text_repel(
    data = label_df,
    aes(label = GeneSymbol),
    size = 2.5,
    color = "black",
    max.overlaps = Inf,
    min.segment.length = 0,
    segment.color = "grey60",
    segment.linewidth = 0.25,
    box.padding = 0.30,
    point.padding = 0.15,
    show.legend = FALSE
  ) +
  geom_text(
    data = cor_df,
    aes(
      x = x,
      y = y,
      label = label
    ),
    inherit.aes = FALSE,
    hjust = -0.15,
    vjust = 1.15,
    size = 2.7,
    color = "black"
  ) +
  facet_wrap(
    ~Comparison,
    nrow = 1
  ) +
  scale_color_manual(
    values = class_cols,
    breaks = legend_order,
    drop = FALSE
  ) +
  theme_classic(base_size = 9) +
  theme(
    text = element_text(family = font_family, color = "black"),
    plot.title = element_text(
      family = font_family,
      face = "bold",
      size = title_size,
      hjust = 0.5,
      color = "black"
    ),
    axis.title = element_text(
      family = font_family,
      color = "black",
      size = axis_title_size
    ),
    axis.text = element_text(
      family = font_family,
      color = "black",
      size = axis_text_size
    ),
    strip.background = element_blank(),
    strip.text = element_text(
      family = font_family,
      color = "black",
      face = "bold",
      size = strip_text_size
    ),
    legend.title = element_blank(),
    legend.text = element_text(
      family = font_family,
      color = "black",
      size = legend_text_size
    ),
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
  ) +
  labs(
    title = "Transcriptome-proteome concordance in lung",
    x = "Transcript log2 fold change",
    y = "Protein log2 fold change"
  )
p_fig4_lung
# ------------------------------------------------------------------------------
# EXPORT
# ------------------------------------------------------------------------------
png(
  filename = file.path(out_dir, "Fig4.png"),
  width = mm_to_in(180),
  height = mm_to_in(85),
  units = "in",
  res = 600,
  bg = "white"
)
print(p_fig4_lung)
dev.off()
pdf(
  file = file.path(out_dir, "Fig4.pdf"),
  width = mm_to_in(180),
  height = mm_to_in(85),
  bg = "white"
)
print(p_fig4_lung)
dev.off()
# ------------------------------------------------------------------------------
# CHECK TABLES
# ------------------------------------------------------------------------------
omics_df %>%
  count(Comparison, regulation_class)
label_df %>%
  arrange(Comparison, desc(label_score)) %>%
  select(
    Comparison,
    GeneSymbol,
    transcript_logFC,
    protein_logFC,
    transcript_pval,
    protein_pval,
    regulation_class,
    label_score
  )
