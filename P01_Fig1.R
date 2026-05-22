# ==============================================================================
# P01 Fig1
# ==============================================================================
# Manuscript project: Acute co-exposure to particulate matter and aircraft noise in the lung-brain-heart axis
# Target journal: Environmental Pollution
# Author / analyst: Corrado Ameli
# Date: 2026-05-22
#
# Purpose:
#   Generate Figure 1 summary panels describing response magnitude, condition effects, directionality, and lung functional categories.
#
# Workflow:
#   1. Read manually annotated differential-abundance workbook.
#   2. Apply region-specific selection thresholds.
#   3. Build Panels A-D and export publication-sized figure files.
#
# ==============================================================================

# ------------------------------------------------------------------------------
# PACKAGES
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(purrr)
  library(ragg)
})

# ------------------------------------------------------------------------------
# INPUT / OUTPUT
# ------------------------------------------------------------------------------
xlsx = "./../../../Documents/MARKOPOLO/P01/DE vs CTRL Biol Categories/DiffExpr_ProDA_Ranked-ARW - annotations.xlsx"
out_dir = "./Current/P01/Results/"
# ------------------------------------------------------------------------------
# THRESHOLDS
# ------------------------------------------------------------------------------
# Region-specific selection rule used in the manuscript figure:
# LU: adj_pval < 0.05; BR/HR: nominal pval < 0.01; all regions: abs(log2FC) >= 0.5
p_cut = 0.01
adj_p_cut = 0.05
lfc_cut = 0.5
regions = c("LU", "BR", "HR")
conditions = c("NIST", "NOISE", "NN")
directions = c("UP", "DOWN")
condition_labels = c(
  NIST = "PM",
  NOISE = "NOISE",
  NN = "PM + NOISE"
)
condition_order = unname(condition_labels[conditions])
condition_cols = c(
  "PM" = "#A67C52",
  "NOISE" = "#6C8E7B",
  "PM + NOISE" = "#8C7AA9"
)
direction_cols = c(
  "UP" = "#8a424a",
  "DOWN" = "#4f6383"
)
# ------------------------------------------------------------------------------
# HELPER FUNCTIONS
# ------------------------------------------------------------------------------
pick_col = function(dat, pattern) {
  hits = names(dat)[grepl(pattern, names(dat))]
  if (length(hits) == 0) stop("Column not found: ", pattern)
  hits[1]
}
read_de = function(sheet) {
  dat = read_excel(xlsx, sheet = sheet, .name_repair = "unique")
  diff_col = pick_col(dat, "^diff(\\.\\.\\.[0-9]+)?$")
  dat %>%
    mutate(
      pval = as.numeric(.data[["pval"]]),
      adj_pval = as.numeric(.data[["adj_pval"]]),
      diff = as.numeric(.data[[diff_col]]),
      Region = as.character(.data[["Region"]]),
      Condition = as.character(.data[["ident.1"]]),
      Significant_calc = case_when(
        Region == "LU" ~ adj_pval < adj_p_cut & abs(diff) >= lfc_cut,
        Region %in% c("BR", "HR") ~ pval < p_cut & abs(diff) >= lfc_cut,
        TRUE ~ FALSE
      ),
      Direction = if_else(diff > 0, "UP", "DOWN")
    )
}
theme_fig = theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.title = element_blank(),
    legend.position = "right"
  )
# ------------------------------------------------------------------------------
# PANELS AC
# ------------------------------------------------------------------------------
de_all = read_de("DE NIST, NOISE, NN vs CTR") %>%
  filter(
    Region %in% regions,
    Condition %in% conditions
  ) %>%
  mutate(
    Region = factor(Region, levels = regions),
    Condition = factor(Condition, levels = conditions),
    Condition_label = factor(
      condition_labels[as.character(Condition)],
      levels = condition_order
    ),
    Direction = factor(Direction, levels = directions)
  )
counts_ab = de_all %>%
  filter(Significant_calc) %>%
  count(Region, Condition_label, name = "n") %>%
  complete(
    Region = regions,
    Condition_label = condition_order,
    fill = list(n = 0)
  ) %>%
  mutate(
    Region = factor(Region, levels = regions),
    Condition_label = factor(Condition_label, levels = condition_order)
  )
# -------------------------
# Panel A
# -------------------------
fig_A = ggplot(counts_ab, aes(x = Region, y = n, fill = Condition_label)) +
  geom_col(width = 0.65) +
  scale_fill_manual(values = condition_cols) +
  labs(
    title = "A. Response magnitude",
    x = "Tissue",
    y = "Number of selected proteins"
  ) +
  theme_fig
# -------------------------
# Panel B
# -------------------------
fig_B = ggplot(counts_ab, aes(x = Region, y = n, fill = Condition_label)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.65) +
  scale_fill_manual(values = condition_cols) +
  labs(
    title = "B. Condition comparison",
    x = "Tissue",
    y = "Number of selected proteins"
  ) +
  theme_fig
# -------------------------
# Panel C
# -------------------------
counts_c = de_all %>%
  filter(Significant_calc, Condition == "NN") %>%
  count(Region, Direction, name = "n") %>%
  complete(
    Region = regions,
    Direction = directions,
    fill = list(n = 0)
  ) %>%
  mutate(
    Region = factor(Region, levels = regions),
    Direction = factor(Direction, levels = directions),
    n_signed = if_else(Direction == "DOWN", -n, n)
  )
fig_C = ggplot(counts_c, aes(x = Region, y = n_signed, fill = Direction)) +
  geom_col(width = 0.65) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  scale_fill_manual(values = direction_cols) +
  labs(
    title = "C. Directionality",
    x = "Tissue",
    y = "Protein count (PM + NOISE)"
  ) +
  theme_fig
# ------------------------------------------------------------------------------
# PANEL D
# ------------------------------------------------------------------------------
# Since Panel D uses LU sheets only, Significant_calc automatically applies:
# adj_pval < 0.05 and abs(diff) >= 0.5
lu_noise = read_de("DE - LU-NOISE")
lu_nist  = read_de("DE - LU-NIST")
lu_nn    = read_de("DE - LU-NN")
cat_col = "biol_category NEW"
if (!cat_col %in% names(lu_noise)) stop(cat_col, " not found in lu_noise")
if (!cat_col %in% names(lu_nist))  stop(cat_col, " not found in lu_nist")
if (!cat_col %in% names(lu_nn))    stop(cat_col, " not found in lu_nn")
lu_noise_cat = lu_noise %>%
  mutate(
    Condition = "NOISE",
    Condition_label = "NOISE",
    Category = trimws(as.character(.data[[cat_col]]))
  )
lu_nist_cat = lu_nist %>%
  mutate(
    Condition = "NIST",
    Condition_label = "PM",
    Category = trimws(as.character(.data[[cat_col]]))
  )
lu_nn_cat = lu_nn %>%
  mutate(
    Condition = "NN",
    Condition_label = "PM + NOISE",
    Category = trimws(as.character(.data[[cat_col]]))
  )
de_lung_cat = bind_rows(
  lu_noise_cat,
  lu_nist_cat,
  lu_nn_cat
)
category_order = c(
  "Immune / inflammatory",
  "Ribosome / translation",
  "RNA processing",
  "Chromatin / nuclear",
  "Metabolic",
  "Oxidative stress",
  "Epithelial / lung function",
  "ECM / structural",
  "Vascular / extracellular signalling"
)
counts_d = de_lung_cat %>%
  filter(
    Significant_calc,
    Category %in% category_order
  ) %>%
  count(Condition_label, Category, Direction, name = "n") %>%
  complete(
    Condition_label = c("NOISE", "PM", "PM + NOISE"),
    Category = category_order,
    Direction = directions,
    fill = list(n = 0)
  ) %>%
  mutate(
    Condition_label = factor(
      Condition_label,
      levels = c("NOISE", "PM", "PM + NOISE")
    ),
    Category = factor(Category, levels = rev(category_order)),
    Direction = factor(Direction, levels = directions)
  )
fig_D = ggplot(counts_d, aes(x = n, y = Category, fill = Direction)) +
  geom_col(width = 0.65) +
  facet_grid(. ~ Condition_label) +
  scale_fill_manual(values = direction_cols) +
  scale_x_continuous(limits = c(0, 90), breaks = seq(0, 90, 20)) +
  labs(
    title = "D. Functional organisation of lung proteomic response",
    x = "Number of selected proteins",
    y = NULL
  ) +
  theme_classic(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    strip.background = element_blank(),
    strip.text = element_text(face = "plain"),
    legend.title = element_blank(),
    legend.position = "bottom",
    panel.spacing.x = unit(1.2, "lines")
  )
# ------------------------------------------------------------------------------
# EXPORT SETTINGS
# ------------------------------------------------------------------------------
out_dir = "./Current/P01/Results/Fig1"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
font_family = "Arial"
base_size = 8.5
dpi_export = 1200
mm_to_in = function(x) x / 25.4
theme_export = theme_classic(base_size = base_size, base_family = font_family) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = base_size + 1),
    axis.title = element_text(size = base_size),
    axis.text = element_text(size = base_size - 1),
    strip.text = element_text(size = base_size, face = "plain"),
    legend.title = element_blank(),
    legend.text = element_text(size = base_size - 1)
  )
save_panel = function(plot, filename, width_mm, height_mm) {
  ggsave(
    filename = file.path(out_dir, filename),
    plot = plot,
    width = mm_to_in(width_mm),
    height = mm_to_in(height_mm),
    units = "in",
    dpi = dpi_export,
    device = ragg::agg_png,
    bg = "white"
  )
}
# ------------------------------------------------------------------------------
# EXPORT WITHOUT LEGENDS
# ------------------------------------------------------------------------------
fig_A_export = fig_A +
  theme_export +
  theme(legend.position = "none")
fig_B_export = fig_B +
  theme_export +
  theme(legend.position = "none")
fig_C_export = fig_C +
  theme_export +
  theme(legend.position = "none")
fig_D_export = fig_D +
  theme_export +
  theme(
    legend.position = "none",
    strip.background = element_rect(fill = NA, colour = NA)
  )
save_panel(fig_A_export, "Fig1A_response_magnitude_nolegend.png", 55, 55)
save_panel(fig_B_export, "Fig1B_condition_comparison_nolegend.png", 55, 55)
save_panel(fig_C_export, "Fig1C_directionality_nolegend.png", 55, 55)
save_panel(fig_D_export, "Fig1D_lung_functional_organisation_nolegend.png", 165, 50)
# ------------------------------------------------------------------------------
# EXPORT WITH LEGENDS
# ------------------------------------------------------------------------------
fig_A_export = fig_A + theme_export
fig_B_export = fig_B + theme_export
fig_C_export = fig_C + theme_export
fig_D_export = fig_D + theme_export
save_panel(fig_A_export, "Fig1A_response_magnitude.png", 55, 55)
save_panel(fig_B_export, "Fig1B_condition_comparison.png", 55, 55)
save_panel(fig_C_export, "Fig1C_directionality.png", 55, 55)
save_panel(fig_D_export, "Fig1D_lung_functional_organisation.png", 165, 50)
