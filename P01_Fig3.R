# ==============================================================================
# P01 Fig3
# ==============================================================================
# Manuscript project: Acute co-exposure to particulate matter and aircraft noise in the lung-brain-heart axis
# Target journal: Environmental Pollution
# Author / analyst: Corrado Ameli
# Date: 2026-05-22
#
# Purpose:
#   Generate Figure 3 additivity / co-exposure response classification panels.
#
# Workflow:
#   1. Read marginal and factorial proDA outputs.
#   2. Calculate expected additivity, observed co-exposure response, and interaction/deviation classes.
#   3. Summarize response classes and functional categories for Figure 3 export.
# ==============================================================================

# ------------------------------------------------------------------------------
# PACKAGES
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readxl)
  library(openxlsx)
  library(cellranger)
  library(stringr)
  library(patchwork)
  library(grid)
})

# ------------------------------------------------------------------------------
# SETTINGS
# ------------------------------------------------------------------------------
use_all_regions = FALSE
region_to_plot = "LU"
out_dir = "./Current/P01/Results/Fig3/"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
factorial_file = "./Current/P01/Results/DiffExpr_ProDA_Factorial.xlsx"
marginal_file  = "./Current/P01/Results/DiffExpr_ProDA.xlsx"
excel_file     = "./Current/P01/Results/Fig3/DAPS_Additivity Filters-ARW-18.04.2026.xlsx"
mm_to_in = function(x) x / 25.4
# Classification thresholds from Methods / Excel logic
nn_response_threshold = 0.35
deviation_threshold = 0.25
interaction_fdr = 0.05
nn_adj_pval_threshold = 0.05
# Desired visual order, top to bottom
response_order_top_to_bottom = c(
  "Proportional ↑",
  "Proportional ↓",
  "Constrained ↑"
)
# ggplot horizontal bars are drawn bottom to top
response_levels = rev(response_order_top_to_bottom)
response_cols = c(
  "Proportional ↑" = "#8a424a",
  "Proportional ↓" = "#4f6383",
  "Constrained ↑"  = "lightgrey"
)
category_levels = c(
  "Ribosome / translation",
  "RNA processing",
  "Chromatin / nuclear",
  "Immune / inflammatory",
  "Metabolic",
  "Vascular / extracellular signalling"
)
category_cols = c(
  "Ribosome / translation" = "#5E6F7F",
  "RNA processing" = "#A76F5E",
  "Chromatin / nuclear" = "#7A6A8F",
  "Immune / inflammatory" = "#6F8A4E",
  "Metabolic" = "#B88746",
  "Vascular / extracellular signalling" = "#5F8A84"
)
font_family = "sans"
title_size = 10
axis_title_size = 8
axis_text_size = 8
legend_text_size = 7
common_theme = theme(
  text = element_text(family = font_family, color = "black"),
  plot.title = element_text(
    family = font_family,
    color = "black",
    face = "bold",
    size = title_size,
    hjust = 0.5
  ),
  axis.title.x = element_text(
    family = font_family,
    color = "black",
    size = axis_title_size
  ),
  axis.text.x = element_text(
    family = font_family,
    color = "black",
    size = axis_text_size
  ),
  axis.text.y = element_text(
    family = font_family,
    color = "black",
    size = axis_text_size
  ),
  axis.title.y = element_blank(),
  axis.line.y = element_blank(),
  axis.ticks.y = element_blank(),
  legend.title = element_blank(),
  legend.text = element_text(
    family = font_family,
    color = "black",
    size = legend_text_size
  )
)
# ------------------------------------------------------------------------------
# CREATE PLOT_DF
# ------------------------------------------------------------------------------
DAPS_marginal = openxlsx::read.xlsx(marginal_file)
DAPS_factorial = openxlsx::read.xlsx(factorial_file)
biol_category = read_excel(
  excel_file,
  sheet = "biol_category",
  range = cell_cols("P:Q"),
  col_types = c("text", "text")
) %>%
  rename(
    GeneSymbol = 1,
    biol_category = 2
  ) %>%
  filter(
    !is.na(GeneSymbol),
    !is.na(biol_category)
  ) %>%
  mutate(
    GeneSymbol = str_trim(GeneSymbol),
    biol_category = str_trim(biol_category)
  ) %>%
  distinct(GeneSymbol, .keep_all = TRUE)
protein_additivity = DAPS_marginal %>%
  filter(ident.1 %in% c("NIST", "NOISE", "NN")) %>%
  transmute(
    Region,
    name,
    GeneSymbol,
    Condition = ident.1,
    diff,
    pval,
    adj_pval
  ) %>%
  pivot_wider(
    names_from = Condition,
    values_from = c(diff, pval, adj_pval),
    names_glue = "{.value}_{Condition}"
  ) %>%
  left_join(
    DAPS_factorial %>%
      filter(ident.1 == "PM:NOISE") %>%
      transmute(
        Region,
        name,
        GeneSymbol,
        interaction_diff = diff,
        interaction_pval = pval,
        interaction_adj_pval = adj_pval
      ),
    by = c("Region", "name", "GeneSymbol")
  )
plot_df = protein_additivity %>%
  mutate(
    GeneSymbol = str_trim(GeneSymbol),
    expected_additive = diff_NIST + diff_NOISE,
    observed_NN = diff_NN,
    deviation_D = observed_NN - expected_additive,
    interaction_minus_D = interaction_diff - deviation_D,
    abs_deviation_D = abs(deviation_D),
    NN_response_strength = abs(diff_NN)
  ) %>%
  left_join(biol_category, by = "GeneSymbol") %>%
  mutate(
    biol_category = ifelse(
      is.na(biol_category),
      "Other / uncertain",
      biol_category
    ),
    Additivity = case_when(
      !is.finite(diff_NN) |
        !is.finite(deviation_D) |
        !is.finite(interaction_adj_pval) |
        !is.finite(adj_pval_NN) ~ NA_character_,
      NN_response_strength < nn_response_threshold ~ NA_character_,
      abs(deviation_D) < deviation_threshold &
        interaction_adj_pval >= interaction_fdr &
        adj_pval_NN < nn_adj_pval_threshold &
        diff_NN > 0 ~ "Proportional ↑",
      abs(deviation_D) < deviation_threshold &
        interaction_adj_pval >= interaction_fdr &
        adj_pval_NN < nn_adj_pval_threshold &
        diff_NN < 0 ~ "Proportional ↓",
      deviation_D <= -deviation_threshold &
        interaction_adj_pval < interaction_fdr &
        diff_NN > 0 ~ "Constrained ↑",
      deviation_D <= -deviation_threshold &
        interaction_adj_pval < interaction_fdr &
        diff_NN < 0 ~ "Constrained ↓",
      deviation_D >= deviation_threshold &
        interaction_adj_pval < interaction_fdr &
        diff_NN > 0 ~ "Amplified ↑",
      deviation_D >= deviation_threshold &
        interaction_adj_pval < interaction_fdr &
        diff_NN < 0 ~ "Amplified ↓",
      TRUE ~ NA_character_
    ),
    Additivity = factor(
      Additivity,
      levels = c(
        "Proportional ↑",
        "Proportional ↓",
        "Constrained ↑",
        "Constrained ↓",
        "Amplified ↑",
        "Amplified ↓"
      )
    )
  )
# ------------------------------------------------------------------------------
# BASE FIGURE DATA
# ------------------------------------------------------------------------------
fig3_df = plot_df
if (!use_all_regions) {
  fig3_df = fig3_df %>%
    filter(Region == region_to_plot)
}
fig3_df = fig3_df %>%
  filter(Additivity %in% response_order_top_to_bottom) %>%
  mutate(
    Additivity = factor(
      as.character(Additivity),
      levels = response_levels
    )
  )
# ------------------------------------------------------------------------------
# PANEL A: DISTRIBUTION OF COEXPOSURE RESPONSES
# ------------------------------------------------------------------------------
panel_a_counts = fig3_df %>%
  count(Additivity, name = "n") %>%
  complete(
    Additivity = factor(response_levels, levels = response_levels),
    fill = list(n = 0)
  ) %>%
  filter(n > 0) %>%
  mutate(
    Additivity = factor(Additivity, levels = response_levels)
  )
p_fig3a = ggplot(
  panel_a_counts,
  aes(
    x = n,
    y = Additivity,
    fill = Additivity
  )
) +
  geom_col(
    width = 0.68,
    color = "white",
    linewidth = 0.4
  ) +
  scale_fill_manual(
    values = response_cols,
    guide = "none"
  ) +
  scale_x_continuous(
    expand = expansion(mult = c(0, 0.12))
  ) +
  coord_cartesian(clip = "off") +
  theme_classic(base_size = 9) +
  common_theme +
  theme(
    plot.margin = margin(5, 10, 5, 5)
  ) +
  labs(
    title = "A. Distribution of co-exposure responses",
    x = "Number of proteins"
  )
# ------------------------------------------------------------------------------
# PANEL B: FUNCTIONAL COMPOSITION
# ------------------------------------------------------------------------------
panel_b_df = fig3_df %>%
  mutate(
    biol_category = case_when(
      biol_category %in% category_levels ~ biol_category,
      str_detect(biol_category, regex("ribosome|translation", ignore_case = TRUE)) ~
        "Ribosome / translation",
      str_detect(biol_category, regex("rna|splic|processing", ignore_case = TRUE)) ~
        "RNA processing",
      str_detect(biol_category, regex("chromatin|nuclear|histone", ignore_case = TRUE)) ~
        "Chromatin / nuclear",
      str_detect(biol_category, regex("immune|inflamm|macrophage|acute", ignore_case = TRUE)) ~
        "Immune / inflammatory",
      str_detect(biol_category, regex("metabolic|metabolism|lipid|mitochond", ignore_case = TRUE)) ~
        "Metabolic",
      str_detect(biol_category, regex("vascular|extracellular|endothelial|matrix|ecm|signalling|signaling", ignore_case = TRUE)) ~
        "Vascular / extracellular signalling",
      TRUE ~ NA_character_
    ),
    biol_category = factor(biol_category, levels = category_levels)
  ) %>%
  filter(!is.na(biol_category))
panel_b_counts = panel_b_df %>%
  count(Additivity, biol_category, name = "n") %>%
  complete(
    Additivity = factor(response_levels, levels = response_levels),
    biol_category = factor(category_levels, levels = category_levels),
    fill = list(n = 0)
  ) %>%
  group_by(Additivity) %>%
  mutate(
    total = sum(n),
    fraction = ifelse(total > 0, n / total, 0)
  ) %>%
  ungroup() %>%
  filter(total > 0)
p_fig3b = ggplot(
  panel_b_counts,
  aes(
    x = fraction,
    y = Additivity,
    fill = biol_category
  )
) +
  geom_col(
    width = 0.68,
    color = "white",
    linewidth = 0.25
  ) +
  scale_y_discrete(
    limits = response_levels
  ) +
  scale_fill_manual(
    values = category_cols,
    breaks = category_levels,
    drop = FALSE
  ) +
  scale_x_continuous(
    labels = function(x) paste0(round(x * 100), "%"),
    breaks = seq(0, 1, 0.2),
    limits = c(0, 1),
    expand = c(0, 0)
  ) +
  theme_classic(base_size = 9) +
  common_theme +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box = "horizontal",
    plot.margin = margin(5, 5, 5, 5)
  ) +
  guides(
    fill = guide_legend(
      nrow = 2,
      byrow = TRUE,
      keywidth = unit(4, "mm"),
      keyheight = unit(4, "mm")
    )
  ) +
  labs(
    title = "B. Functional composition",
    x = "Percentage of proteins"
  )
# =========================
# COMBINED FIGURE 3A-B
# Response-class labels shown once, on Panel A
# Functional-category legend at bottom
# =========================
p_fig3a_combined = p_fig3a +
  theme(
    plot.margin = margin(5, 8, 5, 5)
  )
p_fig3b_combined = p_fig3b +
  theme(
    axis.text.y = element_blank(),
    plot.margin = margin(5, 5, 5, 2)
  )
fig3_ab = p_fig3a_combined + p_fig3b_combined +
  plot_layout(
    widths = c(1.05, 1.45),
    guides = "collect"
  ) &
  theme(
    legend.position = "bottom",
    legend.justification = "center",
    legend.box = "horizontal",
    legend.margin = margin(t = 2, r = 0, b = 0, l = 0),
    legend.text = element_text(
      family = font_family,
      color = "black",
      size = legend_text_size
    )
  )
# ------------------------------------------------------------------------------
# EXPORTS
# ------------------------------------------------------------------------------
png(
  filename = file.path(out_dir, "Fig3A.png"),
  width = mm_to_in(100),
  height = mm_to_in(80),
  units = "in",
  res = 600,
  bg = "white"
)
print(p_fig3a)
dev.off()
pdf(
  file = file.path(out_dir, "Fig3A.pdf"),
  width = mm_to_in(100),
  height = mm_to_in(80),
  bg = "white"
)
print(p_fig3a)
dev.off()
png(
  filename = file.path(out_dir, "Fig3B.png"),
  width = mm_to_in(120),
  height = mm_to_in(80),
  units = "in",
  res = 600,
  bg = "white"
)
print(p_fig3b)
dev.off()
pdf(
  file = file.path(out_dir, "Fig3B.pdf"),
  width = mm_to_in(120),
  height = mm_to_in(80),
  bg = "white"
)
print(p_fig3b)
dev.off()
png(
  filename = file.path(out_dir, "Fig3AB_combined.png"),
  width = mm_to_in(190),
  height = mm_to_in(95),
  units = "in",
  res = 600,
  bg = "white"
)
print(fig3_ab)
dev.off()
pdf(
  file = file.path(out_dir, "Fig3AB_combined.pdf"),
  width = mm_to_in(190),
  height = mm_to_in(95),
  bg = "white"
)
print(fig3_ab)
dev.off()
# ------------------------------------------------------------------------------
# CHECK TABLES
# ------------------------------------------------------------------------------
panel_a_counts
panel_b_counts %>%
  mutate(percent = round(fraction * 100, 1)) %>%
  arrange(Additivity, biol_category)
