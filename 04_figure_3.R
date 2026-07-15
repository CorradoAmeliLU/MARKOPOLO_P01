# ==============================================================================
# Figure 3: co-exposure response classes and pulmonary immune markers
# ==============================================================================

# Packages ---------------------------------------------------------------------

library(dplyr)
library(ggbreak)
library(ggplot2)
library(openxlsx)
library(patchwork)
library(purrr)
library(ragg)
library(stringr)
library(tidyr)

# Paths ------------------------------------------------------------------------

project_dir <- "."

table_dir <- file.path(project_dir, "results", "tables")
figure_dir <- file.path(project_dir, "results", "figure3")
annotation_dir <- file.path(project_dir, "data", "annotations")

marginal_file <- file.path(
  table_dir,
  "proteomics_differential_abundance.xlsx"
)
factorial_file <- file.path(
  table_dir,
  "proteomics_factorial_interaction.xlsx"
)
category_file <- file.path(
  annotation_dir,
  "biological_categories.xlsx"
)

output_png <- file.path(figure_dir, "Figure3.png")
output_pdf <- file.path(figure_dir, "Figure3.pdf")
classification_file <- file.path(
  table_dir,
  "interaction_additivity_classification.xlsx"
)
marker_file <- file.path(table_dir, "figure3_immune_markers.xlsx")

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

# Parameters -------------------------------------------------------------------

region_to_plot <- "LU"

coexposure_response_threshold <- 0.35
deviation_threshold <- 0.25
interaction_fdr_threshold <- 0.05
coexposure_fdr_threshold <- 0.05
maximum_markers <- 25

response_order <- c(
  "Proportional ↑",
  "Proportional ↓",
  "Constrained ↑"
)
response_levels <- rev(response_order)

category_levels <- c(
  "Ribosome / translation",
  "RNA processing",
  "Chromatin / nuclear",
  "Immune / inflammatory",
  "Metabolic",
  "Vascular / extracellular signalling"
)

category_colors <- c(
  "Ribosome / translation" = "#5E6F7F",
  "RNA processing" = "#A76F5E",
  "Chromatin / nuclear" = "#7A6A8F",
  "Immune / inflammatory" = "#6F8A4E",
  "Metabolic" = "#B88746",
  "Vascular / extracellular signalling" = "#5F8A84"
)

condition_colors <- c(
  "PM" = "#A67C52",
  "PM + NOISE" = "#8C7AA9"
)

candidate_marker_genes <- c(
  "Saa3", "Il4i1", "Marco", "Ctsk", "Chi3l1", "Steap4", "Arg2",
  "Lcn2", "Nox3", "Cd14", "Pigr", "Itgad", "Itgax", "Cd68",
  "Sftpd", "Hmox1", "Chil3", "Chil4", "Vnn1", "Cybb", "Acp5",
  "Cyba", "Tcirg1", "Fcer1g", "Mcemp1", "Itgam", "Fgr",
  "Oas1a", "Oas1g", "Sting1", "Ncf2", "Apbb1ip", "Ifit1",
  "Ly75", "Alox5ap", "Ifi204", "Itgb2", "Ifi205a", "Ncf1",
  "Isg20", "Nfkb2", "Myo1f", "Mpeg1", "Ifit2", "Ncf4",
  "Was", "Vav1", "Plcg2", "Lgals3bp", "Skap2", "Lcp2",
  "Lyz1", "Myo1g", "Tbxas1", "Ptges", "Pld3", "Cd63",
  "Plek", "Tmem176a", "Cxcl15", "Ctsc", "Ripk3", "Cd44",
  "S100a4", "Syk", "Ptpn2", "Sp100", "Hcls1", "Npc2",
  "Lpl", "Rab32", "Lyz2", "Lpcat2", "Phf11", "Ptpn6",
  "Cd74", "Hexb", "Ifitm3", "Nckap1l", "Pycard", "Muc1",
  "Grk2", "Grk6", "Arrb2", "Ctsz", "Samd9l", "Parp9",
  "Dock2", "Pla2g4a", "Sp110"
)

font_family <- "sans"
base_size <- 8.5
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

read_category_annotation <- function(path) {
  annotation <- openxlsx::read.xlsx(path)

  category_column <- intersect(
    c("biol_category", "biol_category NEW", "Category"),
    names(annotation)
  )[1]

  if (is.na(category_column)) {
    stop(
      "The category annotation must contain one of these columns: ",
      "biol_category, biol_category NEW, Category."
    )
  }

  annotation$Category <- str_trim(
    as.character(annotation[[category_column]])
  )

  if ("GeneSymbol" %in% names(annotation)) {
    annotation$GeneSymbol <- str_trim(
      as.character(annotation$GeneSymbol)
    )
  }

  if ("name" %in% names(annotation)) {
    return(
      annotation %>%
        select(name, Category) %>%
        filter(!is.na(name), name != "", !is.na(Category), Category != "") %>%
        distinct(name, .keep_all = TRUE)
    )
  }

  if ("GeneSymbol" %in% names(annotation)) {
    return(
      annotation %>%
        select(GeneSymbol, Category) %>%
        filter(
          !is.na(GeneSymbol),
          GeneSymbol != "",
          !is.na(Category),
          Category != ""
        ) %>%
        distinct(GeneSymbol, .keep_all = TRUE)
    )
  }

  stop("The category annotation requires either a name or GeneSymbol column.")
}

add_category_annotation <- function(results, annotation) {
  if ("name" %in% names(annotation)) {
    return(results %>% left_join(annotation, by = "name"))
  }

  if ("GeneSymbol" %in% names(annotation)) {
    return(results %>% left_join(annotation, by = "GeneSymbol"))
  }

  stop("The category annotation requires either a name or GeneSymbol column.")
}

millimeters_to_inches <- function(x) {
  x / 25.4
}

integer_breaks <- function(x) {
  range_x <- range(x, na.rm = TRUE)
  seq(floor(range_x[1]), ceiling(range_x[2]), by = 1)
}

map_broad_category <- function(category) {
  case_when(
    category %in% category_levels ~ category,
    str_detect(
      category,
      regex("ribosome|translation", ignore_case = TRUE)
    ) ~ "Ribosome / translation",
    str_detect(
      category,
      regex("rna|splic|processing", ignore_case = TRUE)
    ) ~ "RNA processing",
    str_detect(
      category,
      regex("chromatin|nuclear|histone", ignore_case = TRUE)
    ) ~ "Chromatin / nuclear",
    str_detect(
      category,
      regex("immune|inflamm|macrophage|acute", ignore_case = TRUE)
    ) ~ "Immune / inflammatory",
    str_detect(
      category,
      regex("metabolic|metabolism|lipid|mitochond", ignore_case = TRUE)
    ) ~ "Metabolic",
    str_detect(
      category,
      regex(
        paste0(
          "vascular|extracellular|endothelial|matrix|ecm|",
          "signalling|signaling"
        ),
        ignore_case = TRUE
      )
    ) ~ "Vascular / extracellular signalling",
    TRUE ~ NA_character_
  )
}

theme_figure <- theme_classic(
  base_size = base_size,
  base_family = font_family
) +
  theme(
    text = element_text(family = font_family, color = "black"),
    plot.title = element_text(
      family = font_family,
      color = "black",
      face = "bold",
      size = base_size + 1,
      hjust = 0.5
    ),
    axis.title = element_text(
      family = font_family,
      color = "black",
      size = base_size
    ),
    axis.text = element_text(
      family = font_family,
      color = "black",
      size = base_size - 1
    ),
    legend.title = element_blank(),
    legend.text = element_text(
      family = font_family,
      color = "black",
      size = base_size - 1
    )
  )

# Input ------------------------------------------------------------------------

walk(c(marginal_file, factorial_file, category_file), check_file_exists)

marginal_results <- openxlsx::read.xlsx(marginal_file)
factorial_results <- openxlsx::read.xlsx(factorial_file)
category_annotation <- read_category_annotation(category_file)

check_columns(
  marginal_results,
  c("Region", "name", "GeneSymbol", "ident.1", "diff", "pval", "adj_pval"),
  "marginal_results"
)
check_columns(
  factorial_results,
  c("Region", "name", "GeneSymbol", "ident.1", "diff", "pval", "adj_pval"),
  "factorial_results"
)

# Additivity classification ----------------------------------------------------

protein_additivity <- marginal_results %>%
  filter(ident.1 %in% c("NIST", "NOISE", "NN"), ident.2 == "CTRL") %>%
  transmute(
    Region,
    name,
    GeneSymbol = str_trim(GeneSymbol),
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
    factorial_results %>%
      filter(ident.1 == "PM:NOISE") %>%
      transmute(
        Region,
        name,
        GeneSymbol = str_trim(GeneSymbol),
        interaction_diff = diff,
        interaction_pval = pval,
        interaction_adj_pval = adj_pval
      ),
    by = c("Region", "name", "GeneSymbol")
  ) %>%
  add_category_annotation(category_annotation)

classification <- protein_additivity %>%
  mutate(
    expected_additive = diff_NIST + diff_NOISE,
    observed_coexposure = diff_NN,
    deviation = observed_coexposure - expected_additive,
    response_strength = abs(observed_coexposure),
    Category = replace_na(Category, "Other / uncertain"),
    Response_class = case_when(
      !is.finite(observed_coexposure) |
        !is.finite(deviation) |
        !is.finite(interaction_adj_pval) ~ NA_character_,
      response_strength < coexposure_response_threshold ~ NA_character_,
      abs(deviation) < deviation_threshold &
        interaction_adj_pval >= interaction_fdr_threshold &
        adj_pval_NN < coexposure_fdr_threshold &
        observed_coexposure > 0 ~ "Proportional ↑",
      abs(deviation) < deviation_threshold &
        interaction_adj_pval >= interaction_fdr_threshold &
        adj_pval_NN < coexposure_fdr_threshold &
        observed_coexposure < 0 ~ "Proportional ↓",
      deviation <= -deviation_threshold &
        interaction_adj_pval < interaction_fdr_threshold &
        observed_coexposure > 0 ~ "Constrained ↑",
      deviation <= -deviation_threshold &
        interaction_adj_pval < interaction_fdr_threshold &
        observed_coexposure < 0 ~ "Constrained ↓",
      deviation >= deviation_threshold &
        interaction_adj_pval < interaction_fdr_threshold &
        observed_coexposure > 0 ~ "Amplified ↑",
      deviation >= deviation_threshold &
        interaction_adj_pval < interaction_fdr_threshold &
        observed_coexposure < 0 ~ "Amplified ↓",
      TRUE ~ NA_character_
    )
  )

openxlsx::write.xlsx(
  classification,
  classification_file,
  overwrite = TRUE
)

figure_data <- classification %>%
  filter(
    Region == region_to_plot,
    Response_class %in% response_order
  ) %>%
  mutate(
    Response_class = factor(Response_class, levels = response_levels)
  )

# Panel A ----------------------------------------------------------------------

panel_a_counts <- figure_data %>%
  mutate(
    Broad_category = map_broad_category(Category),
    Broad_category = factor(Broad_category, levels = category_levels)
  ) %>%
  filter(!is.na(Broad_category)) %>%
  count(Response_class, Broad_category, name = "n") %>%
  complete(
    Response_class = response_levels,
    Broad_category = category_levels,
    fill = list(n = 0)
  ) %>%
  mutate(
    Response_class = factor(Response_class, levels = response_levels),
    Broad_category = factor(Broad_category, levels = category_levels)
  ) %>%
  group_by(Response_class) %>%
  mutate(total = sum(n)) %>%
  ungroup() %>%
  filter(total > 0)

panel_a_totals <- panel_a_counts %>%
  distinct(Response_class, total)

panel_a_limit <- max(panel_a_totals$total, na.rm = TRUE) * 1.04
label_size <- (base_size - 3) / ggplot2::.pt

panel_a <- ggplot(
  panel_a_counts,
  aes(x = n, y = Response_class, fill = Broad_category)
) +
  geom_col(width = 0.68, color = "white", linewidth = 0.25) +
  geom_text(
    data = panel_a_totals,
    aes(
      x = total + panel_a_limit * 0.008,
      y = Response_class,
      label = total
    ),
    inherit.aes = FALSE,
    hjust = 0,
    size = label_size,
    family = font_family
  ) +
  scale_fill_manual(
    values = category_colors,
    breaks = category_levels,
    drop = FALSE
  ) +
  scale_x_continuous(
    limits = c(0, panel_a_limit),
    expand = expansion(mult = c(0, 0))
  ) +
  coord_cartesian(clip = "off") +
  labs(
    title = "A. Co-exposure response classes",
    x = "Number of proteins",
    y = NULL
  ) +
  theme_figure +
  theme(
    legend.position = "right",
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.margin = margin(6, 4, 22, 4)
  )

# Panel B ----------------------------------------------------------------------

panel_b_markers <- classification %>%
  filter(Region == region_to_plot) %>%
  separate_rows(GeneSymbol, sep = ";") %>%
  mutate(GeneSymbol = str_trim(GeneSymbol)) %>%
  filter(
    GeneSymbol %in% candidate_marker_genes,
    is.finite(diff_NIST),
    is.finite(diff_NOISE),
    is.finite(diff_NN)
  ) %>%
  mutate(
    maximum_single_exposure = pmax(diff_NIST, diff_NOISE),
    gain_over_maximum_single = diff_NN - maximum_single_exposure,
    gain_over_pm = diff_NN - diff_NIST,
    gain_over_noise = diff_NN - diff_NOISE
  ) %>%
  filter(
    diff_NN > 0,
    gain_over_maximum_single > 0
  ) %>%
  group_by(GeneSymbol) %>%
  arrange(
    desc(gain_over_maximum_single),
    desc(diff_NN),
    .by_group = TRUE
  ) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  slice_max(
    order_by = gain_over_maximum_single,
    n = maximum_markers,
    with_ties = FALSE
  ) %>%
  arrange(desc(diff_NN), desc(diff_NIST)) %>%
  mutate(
    GeneSymbol = factor(GeneSymbol, levels = GeneSymbol),
    gain_label = paste0("+", sprintf("%.1f", gain_over_maximum_single))
  )

openxlsx::write.xlsx(
  panel_b_markers,
  marker_file,
  overwrite = TRUE
)

panel_b_long <- panel_b_markers %>%
  select(GeneSymbol, diff_NIST, diff_NN) %>%
  pivot_longer(
    cols = c(diff_NIST, diff_NN),
    names_to = "contrast",
    values_to = "effect"
  ) %>%
  mutate(
    Exposure = recode(
      contrast,
      diff_NIST = "PM",
      diff_NN = "PM + NOISE"
    ),
    Exposure = factor(Exposure, levels = c("PM", "PM + NOISE"))
  )

arrow_gap <- 0.10
panel_b_arrows <- panel_b_markers %>%
  mutate(
    effective_gap = pmin(arrow_gap, pmax(gain_over_pm, 0) / 3),
    arrow_start = diff_NIST + effective_gap,
    arrow_end = diff_NN - effective_gap
  )

panel_b <- ggplot(
  panel_b_long,
  aes(x = GeneSymbol, y = effect, color = Exposure)
) +
  geom_segment(
    data = panel_b_arrows,
    aes(
      x = GeneSymbol,
      xend = GeneSymbol,
      y = arrow_start,
      yend = arrow_end
    ),
    inherit.aes = FALSE,
    color = "grey65",
    linewidth = 0.35,
    arrow = grid::arrow(
      length = grid::unit(1.4, "mm"),
      type = "closed"
    )
  ) +
  geom_point(
    size = 2,
    alpha = 0.95,
    position = position_dodge(width = 0.15)
  ) +
  geom_text(
    data = panel_b_markers,
    aes(
      x = GeneSymbol,
      y = diff_NN,
      label = gain_label
    ),
    inherit.aes = FALSE,
    nudge_y = 0.30,
    vjust = 0,
    size = label_size,
    family = font_family,
    color = "black"
  ) +
  scale_color_manual(
    values = condition_colors,
    breaks = c("PM", "PM + NOISE"),
    drop = FALSE
  ) +
  scale_y_continuous(
    breaks = integer_breaks,
    labels = function(x) as.integer(x),
    limits = c(NA, 8.5),
    expand = expansion(mult = c(0.03, 0.18))
  ) +
  ggbreak::scale_y_break(
    c(3, 5.5),
    scales = 0.6,
    ticklabels = c(0, 1, 2, 3, 6, 7, 8)
  ) +
  coord_cartesian(clip = "off") +
  labs(
    title = "B. Immune markers peaking under PM + NOISE",
    x = NULL,
    y = expression("Fold change vs CTRL (" * log[2] * ")")
  ) +
  theme_figure +
  theme(
    legend.position = "right",
    axis.text.x = element_text(
      angle = 40,
      hjust = 1,
      vjust = 1
    ),
    axis.line.y.right = element_blank(),
    axis.ticks.y.right = element_blank(),
    axis.text.y.right = element_blank(),
    plot.margin = margin(6, 4, 22, 4)
  )

# Assemble and export ----------------------------------------------------------

figure_3 <- panel_a + panel_b +
  plot_layout(
    widths = c(1, 1.6),
    guides = "collect"
  ) &
  theme(
    legend.position = "right",
    legend.direction = "vertical",
    legend.box = "vertical"
  )

ggsave(
  output_png,
  plot = figure_3,
  width = millimeters_to_inches(190),
  height = millimeters_to_inches(95),
  units = "in",
  dpi = dpi_export,
  device = ragg::agg_png,
  bg = "white"
)

ggsave(
  output_pdf,
  plot = figure_3,
  width = millimeters_to_inches(190),
  height = millimeters_to_inches(95),
  units = "in",
  device = cairo_pdf,
  bg = "white"
)

message("Figure 3 completed.")
