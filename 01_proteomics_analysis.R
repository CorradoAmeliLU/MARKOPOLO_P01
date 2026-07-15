# ==============================================================================
# Proteomics preprocessing and differential-abundance analysis
# ==============================================================================

set.seed(123)

# Packages ---------------------------------------------------------------------

library(AnnotationDbi)
library(clusterProfiler)
library(dplyr)
library(fgsea)
library(openxlsx)
library(org.Mm.eg.db)
library(purrr)
library(proDA)
library(stringr)
library(tibble)
library(tidyr)

# Paths ------------------------------------------------------------------------

project_dir <- "."

raw_data_dir <- file.path(project_dir, "data", "raw", "proteomics")
metadata_dir <- file.path(project_dir, "data", "metadata")
processed_dir <- file.path(project_dir, "results", "processed")
table_dir <- file.path(project_dir, "results", "tables")
qc_dir <- file.path(project_dir, "results", "qc")
figure2_dir <- file.path(project_dir, "results", "figure2")

brain_file <- file.path(
  raw_data_dir,
  "P01_directLFQ_unfiltered_BR.tsv"
)
heart_file <- file.path(
  raw_data_dir,
  "P01_directLFQ_unfiltered_HT.tsv"
)
lung_file <- file.path(
  raw_data_dir,
  "P01_directLFQ_unfiltered_LU.tsv"
)
metadata_file <- file.path(metadata_dir, "P01_metadata_sdrf.tsv")
gene_mapping_file <- file.path(
  metadata_dir,
  "idmapping_2025_08_25_directLFQ_v20.tsv"
)

walk(
  c(processed_dir, table_dir, qc_dir, figure2_dir),
  dir.create,
  recursive = TRUE,
  showWarnings = FALSE
)

# Parameters -------------------------------------------------------------------

regions <- c("LU", "BR", "HR")
condition_levels <- c("CTRL", "NOISE", "NIST", "NN")

replicate_cv_threshold <- 0.40
missingness_threshold <- 0.30
kegg_min_size <- 15
kegg_max_size <- 250

# Helper functions -------------------------------------------------------------

check_files_exist <- function(paths) {
  missing_paths <- paths[!file.exists(paths)]

  if (length(missing_paths) > 0) {
    stop(
      "Missing required input file(s):\n",
      paste0("- ", missing_paths, collapse = "\n")
    )
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

read_protein_table <- function(path, region) {
  data <- read.delim(
    path,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  names(data)[1] <- "protein"

  data %>%
    pivot_longer(
      cols = -protein,
      names_to = "sample_replicate",
      values_to = "intensity"
    ) %>%
    mutate(Region = region)
}

average_nonmissing <- function(x) {
  if (all(is.na(x))) {
    return(NA_real_)
  }

  mean(x, na.rm = TRUE)
}

map_protein_group_to_gene <- function(protein_group, mapping_vector) {
  accessions <- str_split(protein_group, ";", simplify = FALSE)[[1]]
  gene_symbols <- unname(mapping_vector[accessions])
  gene_symbols[is.na(gene_symbols)] <- accessions[is.na(gene_symbols)]
  paste(gene_symbols, collapse = ";")
}

filter_region_by_missingness <- function(region_data, threshold) {
  protein_columns <- setdiff(
    names(region_data),
    c("Sample", "Region", "Condition")
  )

  missingness <- region_data %>%
    group_by(Condition) %>%
    summarise(
      across(all_of(protein_columns), ~ mean(is.na(.x))),
      .groups = "drop"
    )

  remove_columns <- protein_columns[
    vapply(
      missingness[protein_columns],
      function(x) all(x > threshold),
      logical(1)
    )
  ]

  list(
    data = region_data %>% select(-all_of(remove_columns)),
    removed = tibble(
      Region = rep(unique(region_data$Region), length(remove_columns)),
      protein = remove_columns
    )
  )
}

prepare_expression_matrix <- function(abundance_data, region) {
  region_data <- abundance_data %>%
    filter(Region == region)

  sample_ids <- region_data$Sample
  expression_data <- region_data %>%
    select(-Sample, -Region, -Condition) %>%
    as.matrix()

  expression_matrix <- t(expression_data)
  colnames(expression_matrix) <- sample_ids

  expression_matrix <- expression_matrix[
    rowSums(!is.na(expression_matrix)) > 0,
    ,
    drop = FALSE
  ]

  proDA::median_normalization(log2(expression_matrix))
}

run_marginal_model <- function(abundance_data, region) {
  expression_matrix <- prepare_expression_matrix(abundance_data, region)

  sample_info <- abundance_data %>%
    filter(Region == region) %>%
    transmute(
      Condition = factor(Condition, levels = condition_levels)
    )

  fit <- proDA(
    data = expression_matrix,
    design = ~ Condition,
    col_data = sample_info,
    reference_level = "CTRL",
    verbose = TRUE
  )

  bind_rows(
    test_diff(fit, "ConditionNOISE") %>%
      mutate(ident.1 = "NOISE", ident.2 = "CTRL"),
    test_diff(fit, "ConditionNIST") %>%
      mutate(ident.1 = "NIST", ident.2 = "CTRL"),
    test_diff(fit, "ConditionNN") %>%
      mutate(ident.1 = "NN", ident.2 = "CTRL")
  ) %>%
    mutate(Region = region)
}

run_factorial_model <- function(abundance_data, region) {
  expression_matrix <- prepare_expression_matrix(abundance_data, region)

  sample_info <- abundance_data %>%
    filter(Region == region) %>%
    transmute(
      Condition = factor(Condition, levels = condition_levels),
      PM = factor(if_else(Condition %in% c("NIST", "NN"), 1, 0)),
      NOISE = factor(if_else(Condition %in% c("NOISE", "NN"), 1, 0))
    )

  fit <- proDA(
    data = expression_matrix,
    design = ~ PM * NOISE,
    col_data = sample_info,
    verbose = TRUE
  )

  bind_rows(
    test_diff(fit, "PM1") %>%
      mutate(ident.1 = "PM", ident.2 = "CTRL"),
    test_diff(fit, "NOISE1") %>%
      mutate(ident.1 = "NOISE", ident.2 = "CTRL"),
    test_diff(fit, "`PM1:NOISE1`") %>%
      mutate(ident.1 = "PM:NOISE", ident.2 = "additive")
  ) %>%
    mutate(Region = region)
}

prepare_kegg_gene_sets <- function(cache_file) {
  if (file.exists(cache_file)) {
    return(readRDS(cache_file))
  }

  kegg_data <- clusterProfiler::download_KEGG("mmu")
  kegg_table <- merge(
    kegg_data$KEGGPATHID2EXTID,
    kegg_data$KEGGPATHID2NAME,
    by = "from"
  )
  names(kegg_table) <- c("pathway_id", "ENTREZID", "pathway")

  gene_sets <- split(kegg_table$ENTREZID, kegg_table$pathway)
  saveRDS(gene_sets, cache_file)
  gene_sets
}

run_kegg_enrichment <- function(
  marginal_results,
  gene_map,
  gene_sets,
  region,
  condition
) {
  ranking_table <- marginal_results %>%
    filter(Region == region, ident.1 == condition) %>%
    separate_rows(GeneSymbol, sep = ";") %>%
    mutate(GeneSymbol = str_trim(GeneSymbol)) %>%
    filter(!is.na(GeneSymbol), GeneSymbol != "") %>%
    group_by(GeneSymbol) %>%
    slice_min(order_by = pval, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    left_join(gene_map, by = "GeneSymbol") %>%
    filter(!is.na(ENTREZID)) %>%
    mutate(rank = -log10(pval) * diff) %>%
    distinct(ENTREZID, .keep_all = TRUE)

  ranking_vector <- ranking_table$rank
  names(ranking_vector) <- ranking_table$ENTREZID
  ranking_vector <- sort(ranking_vector, decreasing = TRUE)

  fgsea::fgseaMultilevel(
    pathways = gene_sets,
    stats = ranking_vector,
    minSize = kegg_min_size,
    maxSize = kegg_max_size
  ) %>%
    as_tibble() %>%
    arrange(padj) %>%
    mutate(
      Region = region,
      Condition = condition,
      leadingEdge = map_chr(
        leadingEdge,
        ~ paste(
          unique(
            na.omit(
              gene_map$GeneSymbol[
                match(as.character(.x), gene_map$ENTREZID)
              ]
            )
          ),
          collapse = ";"
        )
      )
    )
}

# Input ------------------------------------------------------------------------

check_files_exist(
  c(
    brain_file,
    heart_file,
    lung_file,
    metadata_file,
    gene_mapping_file
  )
)

metadata <- read.delim(
  metadata_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

gene_mapping <- read.delim(
  gene_mapping_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

check_columns(
  metadata,
  c(
    "source.name",
    "comment.technical.replicate.",
    "characteristics.disease."
  ),
  "metadata"
)
check_columns(gene_mapping, c("From", "To"), "gene_mapping")

protein_to_gene <- setNames(gene_mapping$To, gene_mapping$From)

raw_abundance <- bind_rows(
  read_protein_table(brain_file, "BR"),
  read_protein_table(heart_file, "HR"),
  read_protein_table(lung_file, "LU")
) %>%
  mutate(intensity = na_if(intensity, 0))

metadata <- metadata %>%
  mutate(
    sample_replicate = paste(
      source.name,
      comment.technical.replicate.,
      sep = "_"
    )
  )

raw_abundance <- raw_abundance %>%
  left_join(
    metadata %>%
      transmute(
        sample_replicate,
        Condition = characteristics.disease.
      ),
    by = "sample_replicate"
  ) %>%
  mutate(
    Replicate = str_extract(sample_replicate, "[^_]+$"),
    Sample = str_remove(sample_replicate, "_[^_]+$")
  ) %>%
  select(protein, Sample, Replicate, Region, Condition, intensity)

if (any(is.na(raw_abundance$Condition))) {
  stop("Some proteomics columns could not be matched to the sample metadata.")
}

# Technical-replicate filtering ------------------------------------------------

replicate_cv <- raw_abundance %>%
  group_by(protein, Sample, Condition, Region) %>%
  summarise(
    replicate_mean = mean(intensity, na.rm = TRUE),
    replicate_sd = if_else(
      sum(!is.na(intensity)) > 1,
      sd(intensity, na.rm = TRUE),
      NA_real_
    ),
    replicate_cv = if_else(
      sum(!is.na(intensity)) > 1 & replicate_mean > 0,
      replicate_sd / replicate_mean,
      NA_real_
    ),
    .groups = "drop"
  ) %>%
  group_by(protein) %>%
  summarise(
    mean_cv = mean(replicate_cv, na.rm = TRUE),
    n_observed_pairs = sum(!is.na(replicate_cv)),
    .groups = "drop"
  )

high_cv_proteins <- replicate_cv %>%
  filter(is.finite(mean_cv), mean_cv > replicate_cv_threshold) %>%
  pull(protein)

filtered_replicates <- raw_abundance %>%
  filter(!protein %in% high_cv_proteins)

averaged_abundance <- filtered_replicates %>%
  group_by(protein, Sample, Region, Condition) %>%
  summarise(
    intensity = average_nonmissing(intensity),
    .groups = "drop"
  )

# Missing-value filtering ------------------------------------------------------

wide_abundance <- averaged_abundance %>%
  pivot_wider(names_from = protein, values_from = intensity)

region_filter_results <- wide_abundance %>%
  group_split(Region) %>%
  map(filter_region_by_missingness, threshold = missingness_threshold)

filtered_abundance <- region_filter_results %>%
  map("data") %>%
  bind_rows() %>%
  arrange(factor(Region, levels = regions), Sample)

missingness_removed <- region_filter_results %>%
  map("removed") %>%
  bind_rows()

saveRDS(
  filtered_abundance,
  file.path(processed_dir, "proteomics_filtered_abundance.rds")
)

# Quality-control outputs ------------------------------------------------------

identification_summary <- raw_abundance %>%
  group_by(Region) %>%
  summarise(n_identified = n_distinct(protein), .groups = "drop") %>%
  bind_rows(
    raw_abundance %>%
      summarise(Region = "All", n_identified = n_distinct(protein))
  )

filtering_summary <- tibble(
  filter_step = c(
    "High technical-replicate CV",
    "Missingness threshold"
  ),
  n_removed = c(
    length(high_cv_proteins),
    nrow(missingness_removed)
  )
)

pca_scores <- map_dfr(regions, function(region) {
  expression_matrix <- prepare_expression_matrix(filtered_abundance, region)

  imputed_matrix <- expression_matrix
  row_medians <- apply(imputed_matrix, 1, median, na.rm = TRUE)

  for (row_index in seq_len(nrow(imputed_matrix))) {
    missing_index <- is.na(imputed_matrix[row_index, ])
    imputed_matrix[row_index, missing_index] <- row_medians[row_index]
  }

  row_standard_deviation <- apply(imputed_matrix, 1, sd)
  imputed_matrix <- imputed_matrix[
    is.finite(row_standard_deviation) & row_standard_deviation > 0,
    ,
    drop = FALSE
  ]

  pca_fit <- prcomp(
    t(imputed_matrix),
    center = TRUE,
    scale. = TRUE
  )

  filtered_abundance %>%
    filter(Region == region) %>%
    select(Sample, Region, Condition) %>%
    mutate(
      PC1 = pca_fit$x[, 1],
      PC2 = pca_fit$x[, 2]
    )
})

sample_clustering <- setNames(
  map(regions, function(region) {
    expression_matrix <- prepare_expression_matrix(filtered_abundance, region)

    imputed_matrix <- expression_matrix
    row_medians <- apply(imputed_matrix, 1, median, na.rm = TRUE)

    for (row_index in seq_len(nrow(imputed_matrix))) {
      missing_index <- is.na(imputed_matrix[row_index, ])
      imputed_matrix[row_index, missing_index] <- row_medians[row_index]
    }

    row_standard_deviation <- apply(imputed_matrix, 1, sd)
    imputed_matrix <- imputed_matrix[
      is.finite(row_standard_deviation) & row_standard_deviation > 0,
      ,
      drop = FALSE
    ]

    hclust(dist(t(imputed_matrix)), method = "average")
  }),
  regions
)

write.xlsx(
  list(
    identification_summary = identification_summary,
    filtering_summary = filtering_summary,
    replicate_cv = replicate_cv,
    missingness_removed = missingness_removed,
    pca_scores = pca_scores
  ),
  file.path(qc_dir, "proteomics_qc.xlsx"),
  overwrite = TRUE
)

saveRDS(
  sample_clustering,
  file.path(qc_dir, "sample_hierarchical_clustering.rds")
)

# Figure 2 input ---------------------------------------------------------------

figure2_marker_annotation <- tribble(
  ~Region, ~Module, ~GeneSymbol,
  "LU", "Immune / Macrophage", "Marco",
  "LU", "Immune / Macrophage", "Saa3",
  "LU", "Immune / Macrophage", "Il4i1",
  "LU", "Immune / Macrophage", "Retnla",
  "LU", "Immune / Macrophage", "Chi3l1",
  "LU", "Immune / Macrophage", "Cd68",
  "LU", "Immune / Macrophage", "Itgax",
  "LU", "Iron Handling / Redox", "Fth1",
  "LU", "Iron Handling / Redox", "Hmox1",
  "LU", "Iron Handling / Redox", "Steap4",
  "LU", "Epithelial barrier", "Ces1e",
  "LU", "Epithelial barrier", "Sftpb",
  "LU", "Epithelial barrier", "Scgb1a1",
  "LU", "Epithelial barrier", "Cyp2f2",
  "LU", "Epithelial barrier", "Pigr",
  "LU", "Epithelial barrier", "Sftpd",
  "LU", "Endothelial / ECM", "Cav1",
  "LU", "Endothelial / ECM", "Thbd",
  "LU", "Endothelial / ECM", "Fbn1",
  "LU", "Endothelial / ECM", "Col4a1"
)

protein_columns <- setdiff(
  names(filtered_abundance),
  c("Sample", "Region", "Condition")
)

protein_key <- tibble(protein_col = protein_columns) %>%
  mutate(accession = str_split(protein_col, ";")) %>%
  unnest(accession) %>%
  left_join(gene_mapping, by = c("accession" = "From")) %>%
  filter(!is.na(To)) %>%
  transmute(
    protein_col,
    accession,
    GeneSymbol = To
  ) %>%
  distinct()

figure2_genes <- figure2_marker_annotation$GeneSymbol
figure2_key <- protein_key %>%
  filter(GeneSymbol %in% figure2_genes) %>%
  mutate(GeneSymbol = factor(GeneSymbol, levels = figure2_genes)) %>%
  arrange(GeneSymbol, protein_col)

figure2_values <- filtered_abundance %>%
  filter(Region == "LU") %>%
  select(
    Sample,
    Region,
    Condition,
    all_of(unique(figure2_key$protein_col))
  )

figure2_input <- list(
  LU = list(
    values = figure2_values,
    annotation = figure2_key,
    marker_annotation = figure2_marker_annotation,
    missing_genes = setdiff(
      figure2_genes,
      as.character(figure2_key$GeneSymbol)
    )
  )
)

saveRDS(
  figure2_input,
  file.path(figure2_dir, "figure2_lung_heatmap_input.rds")
)

# Differential-abundance models -----------------------------------------------

marginal_results <- map_dfr(
  regions,
  ~ run_marginal_model(filtered_abundance, .x)
) %>%
  mutate(
    GeneSymbol = map_chr(
      name,
      map_protein_group_to_gene,
      mapping_vector = protein_to_gene
    )
  ) %>%
  arrange(Region, ident.1, adj_pval)

factorial_results <- map_dfr(
  regions,
  ~ run_factorial_model(filtered_abundance, .x)
) %>%
  mutate(
    GeneSymbol = map_chr(
      name,
      map_protein_group_to_gene,
      mapping_vector = protein_to_gene
    )
  ) %>%
  arrange(Region, ident.1, adj_pval)

write.xlsx(
  marginal_results,
  file.path(table_dir, "proteomics_differential_abundance.xlsx"),
  overwrite = TRUE
)

write.xlsx(
  factorial_results,
  file.path(table_dir, "proteomics_factorial_interaction.xlsx"),
  overwrite = TRUE
)

# KEGG enrichment --------------------------------------------------------------

all_gene_symbols <- marginal_results %>%
  separate_rows(GeneSymbol, sep = ";") %>%
  transmute(GeneSymbol = str_trim(GeneSymbol)) %>%
  filter(!is.na(GeneSymbol), GeneSymbol != "") %>%
  distinct() %>%
  pull(GeneSymbol)

gene_map <- AnnotationDbi::select(
  org.Mm.eg.db,
  keys = all_gene_symbols,
  columns = "ENTREZID",
  keytype = "SYMBOL"
) %>%
  filter(!is.na(ENTREZID)) %>%
  distinct(SYMBOL, .keep_all = TRUE) %>%
  transmute(
    GeneSymbol = SYMBOL,
    ENTREZID = as.character(ENTREZID)
  )

kegg_gene_sets <- prepare_kegg_gene_sets(
  file.path(processed_dir, "mouse_kegg_gene_sets.rds")
)

kegg_results <- crossing(
  Region = regions,
  Condition = c("NIST", "NOISE", "NN")
) %>%
  pmap_dfr(
    ~ run_kegg_enrichment(
      marginal_results = marginal_results,
      gene_map = gene_map,
      gene_sets = kegg_gene_sets,
      region = ..1,
      condition = ..2
    )
  )

write.xlsx(
  kegg_results,
  file.path(table_dir, "kegg_enrichment.xlsx"),
  overwrite = TRUE
)

message("Proteomics analysis completed.")
