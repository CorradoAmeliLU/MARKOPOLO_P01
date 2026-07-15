# ==============================================================================
# Figure 1: organ-level response magnitude and lung functional categories
# ==============================================================================

# Packages ---------------------------------------------------------------------

library(dplyr)
library(ggplot2)
library(openxlsx)
library(patchwork)
library(ragg)
library(stringr)
library(tidyr)

# Paths ------------------------------------------------------------------------

project_dir <- "."

table_dir <- file.path(project_dir, "results", "tables")
figure_dir <- file.path(project_dir, "results", "figure1")
annotation_dir <- file.path(project_dir, "data", "annotations")

marginal_file <- file.path(
  table_dir,
  "proteomics_differential_abundance.xlsx"
)
category_file <- file.path(
  annotation_dir,
  "biological_categories.xlsx"
)

output_png <- file.path(figure_dir, "Figure1.png")
output_pdf <- file.path(figure_dir, "Figure1.pdf")
output_tables <- file.path(table_dir, "figure1_counts.xlsx")

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

# Parameters -------------------------------------------------------------------

regions <- c("LU", "BR", "HR")
conditions <- c("NIST", "NOISE", "NN")
directions <- c("UP", "DOWN")

lung_fdr_threshold <- 0.05
distal_p_threshold <- 0.01
fold_change_threshold <- 0.50

condition_labels <- c(
  NIST = "PM",
  NOISE = "NOISE",
  NN = "PM + NOISE"
)
condition_order <- c("PM", "NOISE", "PM + NOISE")

condition_colors <- c(
  "PM" = "#A67C52",
  "NOISE" = "#6C8E7B",
  "PM + NOISE" = "#8C7AA9"
)

direction_colors <- c(
  "UP" = "#8a424a",
  "DOWN" = "#4f6383"
)

category_order <- c(
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

theme_figure <- theme_classic(
  base_size = base_size,
  base_family = font_family
) +
  theme(
    text = element_text(family = font_family, color = "black"),
    plot.title = element_text(
      face = "bold",
      hjust = 0.5,
      size = base_size + 1
    ),
    axis.title = element_text(size = base_size),
    axis.text = element_text(size = base_size - 1),
    strip.background = element_blank(),
    strip.text = element_text(size = base_size, face = "plain"),
    legend.title = element_blank(),
    legend.text = element_text(size = base_size - 1)
  )

# Input and classification -----------------------------------------------------

check_file_exists(marginal_file)
check_file_exists(category_file)

marginal_results <- openxlsx::read.xlsx(marginal_file)
category_annotation <- read_category_annotation(category_file)

check_columns(
  marginal_results,
  c("Region", "ident.1", "pval", "adj_pval", "diff", "GeneSymbol"),
  "marginal_results"
)

classified_results <- marginal_results %>%
  filter(
    Region %in% regions,
    ident.1 %in% conditions,
    ident.2 == "CTRL"
  ) %>%
  mutate(
    Significant = case_when(
      Region == "LU" ~
        adj_pval < lung_fdr_threshold &
        abs(diff) >= fold_change_threshold,
      Region %in% c("BR", "HR") ~
        pval < distal_p_threshold &
        abs(diff) >= fold_change_threshold,
      TRUE ~ FALSE
    ),
    Direction = if_else(diff > 0, "UP", "DOWN"),
    Condition = factor(ident.1, levels = conditions),
    Condition_label = factor(
      condition_labels[as.character(Condition)],
      levels = condition_order
    ),
    Region = factor(Region, levels = regions),
    Direction = factor(Direction, levels = directions)
  )

# Panels A and B ---------------------------------------------------------------

response_counts <- classified_results %>%
  filter(Significant) %>%
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

panel_a <- ggplot(
  response_counts,
  aes(x = Region, y = n, fill = Condition_label)
) +
  geom_col(width = 0.65) +
  scale_fill_manual(values = condition_colors, drop = FALSE) +
  labs(
    title = "A. Response magnitude",
    x = "Tissue",
    y = "Number of selected proteins"
  ) +
  theme_figure +
  theme(
    legend.position = "right",
    plot.margin = margin(5, 5, 5, 5)
  )

panel_b <- ggplot(
  response_counts,
  aes(x = Region, y = n, fill = Condition_label)
) +
  geom_col(
    position = position_dodge(width = 0.75),
    width = 0.65
  ) +
  scale_fill_manual(values = condition_colors, drop = FALSE) +
  labs(
    title = "B. Condition comparison",
    x = "Tissue",
    y = "Number of selected proteins"
  ) +
  theme_figure +
  theme(
    legend.position = "none",
    plot.margin = margin(5, 5, 5, 5)
  )

# Panel C ----------------------------------------------------------------------

direction_counts <- classified_results %>%
  filter(Significant, Condition == "NN") %>%
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

panel_c <- ggplot(
  direction_counts,
  aes(x = Region, y = n_signed, fill = Direction)
) +
  geom_col(width = 0.65) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  scale_fill_manual(values = direction_colors, drop = FALSE) +
  labs(
    title = "C. Directionality",
    x = "Tissue",
    y = "Protein count (PM + NOISE)"
  ) +
  theme_figure +
  theme(
    legend.position = "right",
    plot.margin = margin(5, 5, 5, 5)
  )

# Panel D ----------------------------------------------------------------------

lung_results <- classified_results %>%
  filter(Region == "LU") %>%
  add_category_annotation(category_annotation)

category_counts <- lung_results %>%
  filter(
    Significant,
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

panel_d <- ggplot(
  category_counts,
  aes(x = n, y = Category, fill = Direction)
) +
  geom_col(width = 0.65) +
  facet_grid(. ~ Condition_label) +
  scale_fill_manual(values = direction_colors, drop = FALSE) +
  scale_x_continuous(
    limits = c(0, 90),
    breaks = seq(0, 80, 20)
  ) +
  labs(
    title = "D. Functional organisation of lung proteomic response",
    x = "Number of selected proteins",
    y = NULL
  ) +
  theme_figure +
  theme(
    legend.position = "none",
    panel.spacing.x = grid::unit(1.2, "lines"),
    plot.margin = margin(5, 5, 5, 5)
  )

# Assemble and export ----------------------------------------------------------

figure_1 <- (panel_a | panel_b | panel_c) /
  panel_d +
  plot_layout(heights = c(1, 0.95))

openxlsx::write.xlsx(
  list(
    response_counts = response_counts,
    direction_counts = direction_counts,
    lung_category_counts = category_counts
  ),
  output_tables,
  overwrite = TRUE
)

ggsave(
  output_png,
  plot = figure_1,
  width = millimeters_to_inches(190),
  height = millimeters_to_inches(135),
  units = "in",
  dpi = dpi_export,
  device = ragg::agg_png,
  bg = "white"
)

ggsave(
  output_pdf,
  plot = figure_1,
  width = millimeters_to_inches(190),
  height = millimeters_to_inches(135),
  units = "in",
  device = cairo_pdf,
  bg = "white"
)

message("Figure 1 completed.")
