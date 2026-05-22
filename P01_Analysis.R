# ==============================================================================
# P01 Analysis
# ==============================================================================
# Manuscript project: Acute co-exposure to particulate matter and aircraft noise in the lung-brain-heart axis
# Target journal: Environmental Pollution
# Author / analyst: Corrado Ameli
# Date: 2026-05-22
#
# Purpose:
#   Import, quality-control, filter, normalize, and model directLFQ proteomics data; export proDA differential-abundance and enrichment results.
#
# Workflow:
#   1. Read directLFQ output for brain, heart, and lung.
#   2. Attach SDRF/sample metadata and map protein accessions to gene symbols.
#   3. Inspect technical replicates, remove high-CoV proteins, and average replicate measurements.
#   4. Apply missingness filtering, PCA/clustering QC, and export the Figure 2 matrix.
#   5. Fit proDA marginal and factorial models, then run MSigDB/KEGG enrichment.
#   6. Create exploratory manuscript and supplementary plots used during figure development.
#
# ==============================================================================

# ------------------------------------------------------------------------------
# PACKAGES
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(openxlsx)
  library(proDA)
  library(fgsea)
  library(ggsignif)
  library(ggrepel)
  library(msigdbr)
  library(corrplot)
  library(AnnotationDbi)
  library(org.Mm.eg.db)
  library(clusterProfiler)
  library(KEGGREST)
})

set.seed(123)

# ================================ PARAMETERS ==================================
input_file_BR =
  "./Data/Spectronaut_V20/P01_directLFQ_unfiltered_BR.tsv"
input_file_HR =
  "./Data/Spectronaut_V20/P01_directLFQ_unfiltered_HT.tsv"
input_file_LU =
  "./Data/Spectronaut_V20/P01_directLFQ_unfiltered_LU.tsv"
res_folder =
  "./Results/"
input_file_GeneSymbolMapping =
  "./Data/idmapping_2025_08_25_directLFQ_v20.tsv"
input_file_SampleInfo =
  "./Data/P01_metadata_sdrf.tsv"

# Set TRUE only when MSigDB collections should be downloaded/refreshed.
recompute_dbs = FALSE

dir.create(res_folder, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(res_folder, "fGSEA"), recursive = TRUE, showWarnings = FALSE)
# ------------------------------------------------------------------------------
# ANALYSIS
# ------------------------------------------------------------------------------
# ---- Importing datasets ----
P_BR = read.csv(input_file_BR, sep = "\t") %>% pivot_longer(cols = 2:41,
                                                names_to = "Sample",
                                                values_to = "n")
P_BR$Region = "BR"
P_HR = read.csv(input_file_HR, sep = "\t") %>% pivot_longer(cols = 2:41,
                                                names_to = "Sample",
                                                values_to = "n")
P_HR$Region = "HR"
P_LU = read.csv(input_file_LU, sep = "\t") %>% pivot_longer(cols = 2:41,
                                                names_to = "Sample",
                                                values_to = "n")
P_LU$Region = "LU"
P_raw = rbind(P_BR, P_HR, P_LU)
# ---- Assigning metadata ----
SampleInfo = read.csv(input_file_SampleInfo, sep = "\t")
SampleInfo$SampleReplicate = paste(SampleInfo$source.name,
                                   SampleInfo$comment.technical.replicate.,
                                   sep = "_")
P_raw$Condition = SampleInfo$characteristics.disease.[match(P_raw$Sample, SampleInfo$SampleReplicate)]
P_raw = P_raw %>% separate(Sample,
                           into = c("Sample", "Replicate"),
                           sep = "_")
P_raw$n[P_raw$n==0] = NA
colnames(P_raw)[1] = "protein"
# ---- Load gene names mapping ----
# Export all protein ids
#cat(str_replace_all(paste(unique(P_raw$PG.ProteinGroups), collapse = ";"), ";", "\n"),
#    file = "./Data/all_proteins_v20.txt")
ProteinToGene = read.csv(input_file_GeneSymbolMapping, sep = "\t")
protein_to_gene = setNames(ProteinToGene$To, ProteinToGene$From)
map_protein_to_gene = function(protein_ids_string) {
  protein_ids = unlist(str_split(protein_ids_string, ";"))
  gene_symbols = protein_to_gene[protein_ids]
  gene_symbols[is.na(gene_symbols)] = protein_ids[is.na(gene_symbols)]
  paste(gene_symbols, collapse = ";")
}
# ---- Inspect overall expression at a technical replicate level ----
P_raw %>%
  mutate(Group = paste(Sample, Replicate, sep = "_")) %>%
  mutate(ConditionRegion = paste(Condition, Region, sep = "_")) %>%
  ggplot(aes(x = Group, y = log2(n), fill = ConditionRegion)) +
  geom_boxplot() +
  labs(x = "Sample_Replicate", y = "n", title = "Protein intensities by technical replicate") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
# ---- Inspect correlation between technical replicates ----
# Joined
rep1 = P_raw %>%
  filter(Replicate == "1") %>%
  select(protein, Sample, n) %>%
  dplyr::rename(n_rep1 = n)
rep2 = P_raw %>%
  filter(Replicate == "2") %>%
  select(protein, Sample, n) %>%
  dplyr::rename(n_rep2 = n)
joined = inner_join(rep1, rep2, by = c("protein", "Sample")) %>%
  filter(!is.na(n_rep1) &
           !is.na(n_rep2) & n_rep1 > 0 & n_rep2 > 0) %>%
  mutate(log_n_rep1 = log2(n_rep1),
         log_n_rep2 = log2(n_rep2))
correlation = cor(joined$log_n_rep1, joined$log_n_rep2, method = "pearson")
plot(
  joined$log_n_rep1,
  joined$log_n_rep2,
  xlab = "log2(n) Technical Replicate 1",
  ylab = "log2(n) Technical Replicate 2",
  main = paste("Log2 Correlation:", round(correlation, 3))
)
abline(lm(log_n_rep2 ~ log_n_rep1, data = joined), col = "red")
ggplot(joined, aes(x = log_n_rep1, y = log_n_rep2)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  labs(
    x = "log2(n) Technical Replicate 1",
    y = "log2(n) Technical Replicate 2",
    title = "P1") +
  theme_minimal()
# Correlation histogram per sample
P_wide = P_raw %>%
  filter(Replicate %in% c("1", "2")) %>%
  select(protein, Sample, Replicate, n) %>%
  pivot_wider(names_from = Replicate,
              values_from = n,
              names_prefix = "rep") %>%
  filter(!is.na(rep1) & !is.na(rep2) & rep1 > 0 & rep2 > 0) %>%
  mutate(log_rep1 = log2(rep1), log_rep2 = log2(rep2))
SampleCorr = P_wide %>%
  group_by(Sample) %>%
  summarise(
    n_points = n(),
    pearson_corr = cor(log_rep1, log_rep2, method = "pearson"),
    spearman_corr = cor(log_rep1, log_rep2, method = "spearman")
  )
hist(SampleCorr$pearson_corr, 20, main = "Replicate Correlation per Sample")
# Average CoV by protein (replicate level)
P_CoV = P_raw %>%
  group_by(protein, Sample, Condition, Region) %>%
  summarise(
    mean_n = mean(n, na.rm = TRUE),
    sd_n   = ifelse(sum(!is.na(n)) > 1, sd(n, na.rm = TRUE), NA_real_),
    cv_n   = ifelse(sum(!is.na(n)) > 1, sd(n, na.rm = TRUE) / mean(n, na.rm = TRUE), NA_real_),
    .groups = "drop"
  )  %>%
  group_by(protein) %>%
  summarise(
    mean_cv = mean(cv_n, na.rm = TRUE),
    sd_cv   = sd(cv_n, na.rm = TRUE),
    n_obs   = sum(!is.na(cv_n)),
    .groups = "drop"
  )
hist(P_CoV$mean_cv, 100)
abline(v = 0.4, col = "red", lwd = 2, lty = 2)
ggplot(P_CoV, aes(x = mean_cv)) +
  geom_histogram(bins = 80) +
  geom_vline(xintercept = 0.4, color = "red", linewidth = 1, linetype = "dashed") +
  labs(
    x = "Mean CoV",
    y = "Count",
    title = "P1"
  ) +
  theme_minimal()
plot(P_CoV$mean_cv, P_CoV$sd_cv)
# ---- Filter proteins that exhibit high measurement variation between technical replicates ----
to_remove = P_CoV %>% filter(mean_cv > 0.4) %>% pull(protein)
P_raw = P_raw %>% filter(!protein %in% to_remove)
message("Removed ", length(to_remove), " proteins due to high replicate CoV.")
# ---- Calculate mean quantification between technical replicates ----
#(in the case of one numeric value and one NA value, keep the numeric value)
P = P_raw %>%
  group_by(protein, Sample, Region, Condition) %>%
  summarise(n = if (all(is.na(n))) {
    NA_real_
  } else if (any(is.na(n))) {
    n[!is.na(n)][1]  # keep the non-NA value
  } else {
    mean(n)  # take mean of two values
  }, .groups = "drop")
# ---- Pivoting dataset ----
P_spread = P %>% pivot_wider(names_from = protein, values_from = n)
protein_idxs = c(4:dim(P_spread)[2])
metadata_idxs = c(1:3)
# ---- Removing entries that exhibit too many NAs ----
#(within one region, for each protein, if we observe that in all conditions
#there are >30% missing values, we remove those quantification)
protein_cols = colnames(P_spread)[4:ncol(P_spread)]
clean_region = function(df_region) {
  n_before = length(protein_cols)
  na_props = df_region %>%
    group_by(Condition) %>%
    summarise(across(all_of(protein_cols), ~ mean(is.na(.))), .groups = "drop") %>%
    column_to_rownames("Condition")
  n_already_NA = length(which(colSums(na_props==1) == nrow(na_props)))
  cols_to_na = names(which(colSums(na_props > 0.3) == nrow(na_props)))
  df_region = df_region %>%
    mutate(across(all_of(cols_to_na), ~ NA_real_))
  all_na_cols = df_region %>%
    select(all_of(protein_cols)) %>%
    select(where(~ all(is.na(.)))) %>%
    colnames()
  n_removed = length(all_na_cols) - n_already_NA
  message("Removed ", n_removed, " proteins (",
          round(100 * n_removed / n_before, 1), "%) due to >30% missing values in all conditions")
  df_region = df_region %>%
    select(where(~ !all(is.na(.))))
  return(df_region)
}
P_spread = P_spread %>%
  group_split(Region) %>%
  map_df(clean_region) %>%
  select(where(~!all(is.na(.))))
protein_idxs = c(4:dim(P_spread)[2])
# ---- Inspect PCA ----
abund_matrix =  P_spread[, 4:ncol(P_spread)] %>%
                as.matrix() %>%
                t() %>%
                log2() %>%
                median_normalization()
pca_res = pcaMethods::pca(t(abund_matrix), scale = "vector", center = TRUE)
P_pca = cbind(P_spread[, 1:3], pca_res@scores[, 1], pca_res@scores[, 2])
colnames(P_pca)[4:5] = c("PCA1", "PCA2")
ggplot(P_pca, aes(x = PCA1, y = PCA2, color = Condition, label = Sample)) +
  geom_point(size = 3) +
  #geom_text(vjust = -0.5, size = 3) +  # Add sample labels
  facet_wrap(~ Region, scales = "free") +              # Facet by Region
  theme_minimal() +
  labs(
    title = "P1",
    x = "PC1",
    y = "PC2"
  ) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold"),
    axis.text.x = element_blank(),
    axis.text.y = element_blank()
  )
# ---- Inspect sample dendrogram ----
sampleTree = hclust(dist(P_spread[, protein_idxs]), method = "average")
sampleTree$labels = paste(P_spread$Sample, P_spread$Region, P_spread$Condition)
plot(
  sampleTree,
  main = "Sample clustering",
  sub = "",
  xlab = "",
  cex.lab = 1.5,
  cex.axis = 1.5,
  cex.main = 2
)
# ---- Filtering samples ----
# No sample to filter
# ---- Check correlation between number of NAs and mean expression
P_spread %>%
  mutate(mean_protein = rowMeans(select(., all_of(protein_idxs)), na.rm = TRUE),
         na_count = rowSums(is.na(select(
           ., all_of(protein_idxs)
         )))) %>%
  ggplot(aes(x = mean_protein, y = na_count, color = Region)) +
  geom_point() +
  labs(x = "Mean protein abundance",
       y = "Count of NA values",
       title = "Mean Protein Abundance vs. Number of NA by Region") +
  theme_minimal()
# ---- Export matrix for Figure 2 Heatmap ----
out_dir = "./Results/Fig2/"
genes_to_print = list(
  LU = c(
    # Immune / macrophage / acute-phase
    "Marco",
    "Saa3",
    "Il4i1",
    "Retnla",
    "Chi3l1",
    "Cd68",
    "Itgax",
    # Redox / iron handling
    "Hmox1",
    "Steap4",
    "Fth1",
    # Epithelial / barrier / surfactant
    "Ces1e",
    "Sftpb",
    "Scgb1a1",
    "Cyp2f2",
    "Pigr",
    "Sftpd",
    # Endothelial / vascular / ECM
    "Cav1",
    "Thbd",
    "Fbn1",
    "Col4a1"
  )
)
figure2_marker_annotation_lung = tibble::tribble(
  ~Region, ~Module, ~GeneSymbol,
  # Immune / macrophage / acute-phase
  "LU", "Immune / macrophage / acute-phase", "Marco",
  "LU", "Immune / macrophage / acute-phase", "Saa3",
  "LU", "Immune / macrophage / acute-phase", "Il4i1",
  "LU", "Immune / macrophage / acute-phase", "Retnla",
  "LU", "Immune / macrophage / acute-phase", "Chi3l1",
  "LU", "Immune / macrophage / acute-phase", "Cd68",
  "LU", "Immune / macrophage / acute-phase", "Itgax",
  # Redox / iron handling
  "LU", "Redox / iron handling", "Hmox1",
  "LU", "Redox / iron handling", "Steap4",
  "LU", "Redox / iron handling", "Fth1",
  # Epithelial / barrier / surfactant
  "LU", "Epithelial / barrier / surfactant", "Ces1e",
  "LU", "Epithelial / barrier / surfactant", "Sftpb",
  "LU", "Epithelial / barrier / surfactant", "Scgb1a1",
  "LU", "Epithelial / barrier / surfactant", "Cyp2f2",
  "LU", "Epithelial / barrier / surfactant", "Pigr",
  "LU", "Epithelial / barrier / surfactant", "Sftpd",
  # Endothelial / vascular / ECM
  "LU", "Endothelial / vascular / ECM", "Cav1",
  "LU", "Endothelial / vascular / ECM", "Thbd",
  "LU", "Endothelial / vascular / ECM", "Fbn1",
  "LU", "Endothelial / vascular / ECM", "Col4a1"
)
figure2_markers_lung = list(
  LU = figure2_marker_annotation_lung$GeneSymbol
)
meta_cols = c("Sample", "Region", "Condition")
protein_cols = setdiff(names(P_spread), meta_cols)
protein_key = tibble(protein_col = protein_cols) %>%
  mutate(From = strsplit(protein_col, ";")) %>%
  unnest(From) %>%
  left_join(ProteinToGene, by = "From") %>%
  filter(!is.na(To)) %>%
  transmute(
    protein_col,
    accession = From,
    GeneSymbol = To
  ) %>%
  distinct()
heatmap_original_values = imap(genes_to_print, function(gene_vec, reg) {
  key_reg = protein_key %>%
    filter(GeneSymbol %in% gene_vec) %>%
    mutate(GeneSymbol = factor(GeneSymbol, levels = gene_vec)) %>%
    arrange(GeneSymbol, protein_col)
  cols_reg = unique(key_reg$protein_col)
  values_reg = P_spread %>%
    filter(Region == reg) %>%
    select(all_of(meta_cols), all_of(cols_reg))
  marker_annotation_reg = figure2_marker_annotation_lung %>%
    filter(Region == reg, GeneSymbol %in% gene_vec)
  list(
    values = values_reg,
    annotation = key_reg,
    marker_annotation = marker_annotation_reg,
    missing_genes = setdiff(gene_vec, key_reg$GeneSymbol)
  )
})
saveRDS(
  heatmap_original_values,
  file = file.path(out_dir, "P_spread_subset_lung_Fig2.rds")
)
# ---- Differential-abundance analysis with proDA (factorial model: PM * NOISE) ----
DEGS = data.frame()
for (region in c("LU", "BR", "HR")) {
  cat("\n")
  cat(region)
  cat("\n")
  expr_matrix = P_spread %>%
    filter(Region == region) %>%
    select(-Condition, -Region, -Sample) %>%
    as.matrix() %>%
    t() %>%
    log2()
  expr_matrix = expr_matrix[which(rowSums(is.na(expr_matrix))!=dim(expr_matrix)[2]), ]
  expr_matrix = proDA::median_normalization(expr_matrix)
  sample_info = P_spread %>%
    filter(Region == region) %>%
    select(Condition)
  sample_info$Condition = factor(sample_info$Condition,
                                 levels = c("CTRL", "NOISE", "NIST", "NN"))
  # 2x2 factorial encoding:
  # PM = 1 for NIST and NN
  # NOISE = 1 for NOISE and NN
  sample_info$PM    = factor(ifelse(sample_info$Condition %in% c("NIST", "NN"), 1, 0))
  sample_info$NOISE = factor(ifelse(sample_info$Condition %in% c("NOISE", "NN"), 1, 0))
  fit = proDA(
    data = expr_matrix,
    design = ~ PM * NOISE,
    col_data = sample_info,
    verbose = TRUE
  )
  DEG_PM = test_diff(fit, "PM1")
  DEG_PM$ident.1 = "PM"
  DEG_PM$ident.2 = "CTRL"
  DEG_NOISE = test_diff(fit, "NOISE1")
  DEG_NOISE$ident.1 = "NOISE"
  DEG_NOISE$ident.2 = "CTRL"
  DEG_INTERACTION = proDA::test_diff(fit, "`PM1:NOISE1`")
  DEG_INTERACTION$ident.1 = "PM:NOISE"
  DEG_INTERACTION$ident.2 = "additive"
  temp = rbind(DEG_PM, DEG_NOISE, DEG_INTERACTION)
  temp$Region = region
  DEGS = rbind(DEGS, temp)
}
DEGS = DEGS %>%
  mutate(GeneSymbol = map_chr(name, map_protein_to_gene)) %>%
  arrange(adj_pval)
write.xlsx(DEGS, file = paste0(res_folder, "DiffExpr_ProDA_Factorial.xlsx"))
# ---- Differential-abundance analysis with proDA (marginal contrasts) ----
DEGS = data.frame()
for (region in c("LU", "BR", "HR")) {
  cat("\n")
  cat(region)
  cat("\n")
  expr_matrix = P_spread %>%
    filter(Region == region) %>%
    select(-Condition, -Region, -Sample) %>%
    as.matrix() %>%
    t() %>%
    log2()
  expr_matrix = expr_matrix[which(rowSums(is.na(expr_matrix))!=dim(expr_matrix)[2]), ]
  expr_matrix = proDA::median_normalization(expr_matrix)
  sample_info = P_spread %>%
    filter(Region == region) %>%
    select(Condition)
  sample_info$Condition = factor(sample_info$Condition,
                                 levels = c("CTRL", "NOISE", "NIST", "NN"))
  fit = proDA(
    data = expr_matrix,
    design = ~ Condition,
    col_data = sample_info,
    verbose = TRUE,
    reference_level = "CTRL"
  )
  DEG_NOISE = test_diff(fit, "ConditionNOISE")
  DEG_NOISE$ident.1 = "NOISE"
  DEG_NOISE$ident.2 = "CTRL"
  DEG_NIST = test_diff(fit, "ConditionNIST")
  DEG_NIST$ident.1 = "NIST"
  DEG_NIST$ident.2 = "CTRL"
  DEG_NN = test_diff(fit, "ConditionNN")
  DEG_NN$ident.1 = "NN"
  DEG_NN$ident.2 = "CTRL"
  temp = rbind(DEG_NOISE, DEG_NIST, DEG_NN)
  temp$Region = region
  DEGS = rbind(DEGS, temp)
}
DEGS = DEGS %>%
  mutate(GeneSymbol = map_chr(name, map_protein_to_gene)) %>%
  arrange(adj_pval)
write.xlsx(DEGS, file = paste0(res_folder, "DiffExpr_ProDA.xlsx"))
# ---- Differential-abundance analysis with proDA (marginal contrasts NN vs stressors) ----
DEGS_NN = data.frame()
for (region in c("LU", "BR", "HR")) {
  cat("\n")
  cat(region)
  cat("\n")
  expr_matrix = P_spread %>%
    filter(Region == region) %>%
    select(-Condition, -Region, -Sample) %>%
    as.matrix() %>%
    t() %>%
    log2()
  expr_matrix = expr_matrix[which(rowSums(is.na(expr_matrix))!=dim(expr_matrix)[2]), ]
  expr_matrix = proDA::median_normalization(expr_matrix)
  sample_info = P_spread %>%
    filter(Region == region) %>%
    select(Condition)
  sample_info$Condition = factor(sample_info$Condition,
                                 levels = c("CTRL", "NOISE", "NIST", "NN"))
  fit = proDA(
    data = expr_matrix,
    design = ~ Condition,
    col_data = sample_info,
    verbose = TRUE,
    reference_level = "NN"
  )
  DEG_NOISE = test_diff(fit, "ConditionNOISE")
  DEG_NOISE$ident.1 = "NOISE"
  DEG_NOISE$ident.2 = "NN"
  DEG_NIST = test_diff(fit, "ConditionNIST")
  DEG_NIST$ident.1 = "NIST"
  DEG_NIST$ident.2 = "NN"
  temp = rbind(DEG_NOISE, DEG_NIST)
  temp$Region = region
  DEGS_NN = rbind(DEGS_NN, temp)
}
DEGS_NN = DEGS_NN %>%
  mutate(GeneSymbol = map_chr(name, map_protein_to_gene)) %>%
  arrange(adj_pval)
write.xlsx(DEGS_NN, file = paste0(res_folder, "DiffExpr_ProDA_NOISE_NIST_VS_NN.xlsx"))
# ---- Boxplot of a specific protein/protein group ----
P_spread %>%
  mutate(ConditionRegion = paste(Region, Condition)) %>%
  ggplot(aes(x = ConditionRegion, y = log2(`Q64314`))) +
  geom_boxplot(fill = "steelblue", alpha = 0.7) +
  geom_jitter(width = 0.2,
              alpha = 0.4,
              color = "black") +
  labs(title = "Protein Intensity by Condition", y = "Intensity", x = "Condition") +
  theme_minimal() +
  theme(text = element_text(size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1))
P_spread %>%
  group_by(Region) %>%
  mutate(Q64314_scaled = `Q61093` / mean(`Q61093`[Condition == "CTRL"], na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(ConditionRegion = paste(Region, Condition)) %>%
  ggplot(aes(x = ConditionRegion, y = log2(Q64314_scaled))) +
  geom_boxplot(fill = "steelblue", alpha = 0.7) +
  geom_jitter(width = 0.2,
              alpha = 0.4,
              color = "black") +
  labs(title = "Protein Intensity by Condition (scaled to CTRL mean per Region)",
       y = "log2(Intensity scaled to CTRL mean)",
       x = "Region & Condition") +
  theme_minimal() +
  theme(text = element_text(size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1))
# ---- Boxplot of top 12 proDA results ----
to_plot = DEGS %>%
  mutate(rank= -log10(pval)*diff) %>%
  group_by(ident.1, Region) %>%
  arrange(-rank, .by_group = TRUE) %>%
  slice_head(n = 12) %>%
  ungroup()
for (region in c("LU", "BR", "HR")) {
  for (cond in c("NOISE", "NIST", "NN")) {
    to_plot_sub = to_plot %>%
      filter(Region == region) %>%
      filter(ident.1 == cond) %>%
      pull(name)
    p = P_spread %>%
      pivot_longer(
        cols = protein_idxs,
        names_to = "protein",
        values_to = "n") %>%
      filter(protein %in% to_plot_sub)  %>%
      filter(Region == region) %>%
      filter(Condition == cond | Condition == "CTRL" ) %>%
      mutate(protein = factor(protein, levels = to_plot_sub)) %>%
      mutate(ConditionRegion = paste(Region, Condition)) %>%
      ggplot(aes(x = ConditionRegion, y = log2(n))) +
      geom_boxplot(fill = "steelblue", alpha = 0.7, outlier.color=NA) +
      geom_jitter(width = 0.2,
                  alpha = 0.4,
                  color = "black") +
      facet_wrap(~ protein, scales = "free_y") +
      labs(title = paste0("Protein Intensity by Condition ", region, " ", cond), y = "Intensity", x = "Condition") +
      theme_minimal() +
      theme(text = element_text(size = 12),
            axis.text.x = element_text(angle = 45, hjust = 1))
    ggsave(
      filename = paste0(res_folder, "Top12_boxplot_", region, "_", cond, ".pdf"),
      plot = p,
      width = 8, height = 8
    )
  }
}
# ---- Pull databases for enrichment ----
get_na_pos_neg = function(core_genes_str, gene_list) {
  genes = unlist(strsplit(core_genes_str, "/"))
  genes_in_rank = gene_list[genes]
  na_pos = sum(genes_in_rank > 0)
  na_neg = sum(genes_in_rank < 0)
  return(c(na_pos = na_pos, na_neg = na_neg))
}
if(recompute_databases){
  for(collect in c("M1", "M2", "M3", "M5", "MH")) {
    cat("\n")
    cat(collect)
    cat("\n")
    msigdb_sets = msigdbr(species = "mouse", collection = collect)
    write.xlsx(msigdb_sets, file = paste0("./Current/Extra/MSIGDB/", collect, "_db.xlsx"))
  }
}
# ---- Enrichment with fgsea on MSigDB ----
DAPS = read.xlsx("./Results/DiffExpr_ProDA.xlsx")
for (reg in c("LU", "BR", "HR")) {
  for (cond in c("NIST", "NOISE", "NN")){
    DEGS_ranked = DAPS %>%
      filter(ident.1==cond) %>%
      filter(Region==reg) %>%
      mutate(rank = -log10(pval)*diff) %>%
      dplyr::select(GeneSymbol, rank) %>%
      distinct(GeneSymbol, .keep_all = TRUE) %>%
      arrange(desc(rank))
    ranking_vector = DEGS_ranked$rank
    names(ranking_vector) = DEGS_ranked$GeneSymbol
    ranking_vector = sort(ranking_vector, decreasing = TRUE)
    ranking_vector = ranking_vector[!is.na(names(ranking_vector)) & !is.na(ranking_vector)]
    for(collect in c("M1", "M2", "M3", "M5", "MH")) {
      cat("\n")
      cat(paste(reg, cond, collect))
      cat("\n")
      msigdb_sets = read.xlsx(paste0("./Current/Extra/MSIGDB/", collect, "_db.xlsx"))
      gene_sets = msigdb_sets %>%
        dplyr::select(gs_name, gene_symbol) %>%
        group_by(gs_name) %>%
        summarise(genes = list(unique(gene_symbol))) %>%
        deframe()
      gsea_results = fgsea::fgsea(pathways = gene_sets,
                                  stats = ranking_vector,
                                  minSize = 15,
                                  maxSize = 250,
                                  nproc = 1)
      gsea_results = gsea_results[order(gsea_results$padj), ]
      #gsea_results$go_id = gmt_file$gs_exact_source[match(gsea_results$pathway, gmt_file$gs_name)]
      res_df = as.data.frame(gsea_results)
      write.xlsx(res_df, file = paste0(res_folder, "fGSEA/", collect, "_", reg, "_", cond, ".xlsx"))
    }
  }
}
# ---- Enrichment with fgsea (KEGG, Mouse) ----
DAPS = read.xlsx("./Results/DiffExpr_ProDA.xlsx")
res_folder = "./Results/"
symbols = unique(DAPS$GeneSymbol)
gene_map = AnnotationDbi::select(
  org.Mm.eg.db,
  keys = symbols,
  columns = "ENTREZID",
  keytype = "SYMBOL"
) %>%
  filter(!is.na(ENTREZID)) %>%
  distinct(SYMBOL, .keep_all = TRUE) %>%
  rename(GeneSymbol = SYMBOL)
DAPS_mapped = DAPS %>%
  left_join(gene_map, by = "GeneSymbol")
kegg_sets = {
  ks = clusterProfiler::download_KEGG("mmu")
  merged = merge(ks$KEGGPATHID2EXTID, ks$KEGGPATHID2NAME, by = "from")
  colnames(merged) = c("pathway_id", "ENTREZID", "pathway_name")
  split(merged$ENTREZID, merged$pathway_name)
}
kegg_list = list()
for (reg in c("LU", "BR", "HR")) {
  for (cond in c("NIST", "NOISE", "NN")) {
    DEGS_ranked = DAPS_mapped %>%
      filter(ident.1 == cond, Region == reg) %>%
      mutate(rank = -log10(pval) * diff) %>%
      select(GeneSymbol, ENTREZID, rank) %>%
      distinct(GeneSymbol, .keep_all = TRUE) %>%
      filter(!is.na(ENTREZID)) %>%
      arrange(desc(rank))
    ranking_vector = DEGS_ranked$rank
    names(ranking_vector) = DEGS_ranked$ENTREZID
    ranking_vector = ranking_vector[!is.na(names(ranking_vector)) & !is.na(ranking_vector)]
    ranking_vector = sort(ranking_vector, decreasing = TRUE)
    cat("\nRunning KEGG fgsea for:", reg, cond, "\n")
    gsea_results = fgsea::fgsea(
      pathways = kegg_sets,
      stats = ranking_vector,
      minSize = 15,
      maxSize = 250,
      nproc = 1) %>%
      arrange(padj) %>%
      mutate(leadingEdge = map_chr(leadingEdge, ~ paste(unique(na.omit(
        gene_map$GeneSymbol[match(.x, gene_map$ENTREZID)]
      )), collapse=", ")))
    write.xlsx(as.data.frame(gsea_results),
               file = paste0(res_folder, "/fGSEA/KEGG_", reg, "_", cond, ".xlsx"))
    kegg_list[[paste(reg, cond, sep="_")]] = gsea_results
  }
}
# ---- Output genes and foldchanges for KEGG Pathway visualization ----
kegg_map = bitr(ProteinToGene$To, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
ProteinToGene$ToKegg = paste0("mmu:", kegg_map$ENTREZID[match(ProteinToGene$To, kegg_map$SYMBOL)])
path_id = names(keggList("pathway", "mmu"))[grep("AGE-RAGE signaling pathway in diabetic complications", keggList("pathway", "mmu"), ignore.case=TRUE)[1]]
genes = keggLink("mmu", path_id)
DAPS$ToKegg = ProteinToGene$ToKegg[match(DAPS$GeneSymbol, ProteinToGene$To)]
kegg_output = DAPS %>%
  filter(ToKegg %in% genes) %>%
  filter(ident.1=="NIST") %>%
  filter(Region == "LU") %>%
  dplyr::select(ToKegg, diff) %>%
  mutate(diff = round(diff, 2))
cat(apply(kegg_output, 1, \(x) paste(x[1], sprintf("%.3f", as.numeric(x[2])))), sep="\n")
# ------------------------------------------------------------------------------
# MANUSCRIPT PLOTS
# ------------------------------------------------------------------------------
# ---- Plot volcano plots vs CTRL----
DEGS = read.xlsx("./Results/DiffExpr_ProDA.xlsx")
# LUNGS
region = "LU"
DEGS_sub = DEGS %>%
  filter(Region == region) %>%
  mutate(
    significant = adj_pval < 0.05 & abs(diff) >= 0.5,
    diff_capped = pmax(pmin(diff, 3), -3),
    score = -log10(pval) * abs(diff)
  )
top_labels = DEGS_sub %>%
  filter(abs(diff) > 0.5, adj_pval < 0.05) %>%   # <- add here
  group_by(ident.1) %>%
  filter(diff > 0) %>%
  slice_max(order_by = score, n = 8, with_ties = FALSE) %>%
  ungroup() %>%
  bind_rows(
    DEGS_sub %>%
      filter(abs(diff) > 0.5, adj_pval < 0.05) %>%  # <- and here
      group_by(ident.1) %>%
      filter(diff < 0) %>%
      slice_max(order_by = score, n = 8, with_ties = FALSE) %>%
      ungroup()
  )
ggplot(DEGS_sub, aes(
  x = diff_capped,
  y = -log10(pval),
  color = significant)) +
  geom_point(alpha = 0.8) +
  geom_vline(xintercept = c(-1, 1),
             colour = "grey70",
             linetype = "dotted") +
  geom_hline(yintercept = 2,
             colour = "grey70",
             linetype = "dotted") +
  geom_text_repel(
    data = top_labels,
    aes(label = GeneSymbol),
    color = "black",
    size = 3.5,
    max.overlaps = Inf,
    box.padding = 0.4,
    point.padding = 0.4,
    segment.color = "grey70",
    segment.size = 0.5) +
  facet_wrap(~ident.1) +
  scale_color_manual(values = c("grey80", "#74a892")) +
  scale_x_continuous(
    limits = c(-3, 3),
    breaks = c(-1, 0, 1)) +
  scale_y_continuous(
    limits = c(0, 10),
    breaks = c(0, 2, 10)) +
  theme_minimal(base_size = 16) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_blank(),
    legend.position = "none") +
  labs(
    x = expression("fold change ("*log[2]*")"),
    y = expression("P value (- "*log[10]*")"))
# BRAIN
region = "BR"
DEGS_sub = DEGS %>%
  filter(Region == region) %>%
  mutate(
    significant = pval < 0.01 & abs(diff) >= 0.5,
    diff_capped = pmax(pmin(diff, 3), -3),
    score = -log10(pval) * abs(diff)
  )
top_labels = DEGS_sub %>%
  filter(abs(diff) > 0.5) %>%   # FC threshold only
  group_by(ident.1) %>%
  # top positive side
  filter(diff > 0) %>%
  slice_max(order_by = score, n = 8, with_ties = FALSE) %>%
  ungroup() %>%
  bind_rows(
    DEGS_sub %>%
      filter(abs(diff) > 0.5) %>%   # FC threshold only
      group_by(ident.1) %>%
      # top negative side
      filter(diff < 0) %>%
      slice_max(order_by = score, n = 8, with_ties = FALSE) %>%
      ungroup()
  )
top_labels$GeneSymbol[top_labels$GeneSymbol=="Atp5mc2;Atp5mc3;Atp5mc1"] = "Atp5"
top_labels = top_labels %>% filter(pval < 0.01)
ggplot(DEGS_sub, aes(
  x = diff_capped,
  y = -log10(pval),
  color = significant)) +
  geom_point(alpha = 0.8) +
  geom_vline(xintercept = c(-1, 1),
             colour = "grey70",
             linetype = "dotted") +
  geom_hline(yintercept = 2,
             colour = "grey70",
             linetype = "dotted") +
  geom_text_repel(
    data = top_labels,
    aes(label = GeneSymbol),
    color = "black",
    size = 3.5,
    max.overlaps = Inf,
    box.padding = 0.4,
    point.padding = 0.4,
    segment.color = "grey70",
    segment.size = 0.5) +
  facet_wrap(~ident.1) +
  scale_color_manual(values = c("grey80", "#e5c185")) +
  scale_x_continuous(
    limits = c(-3, 3),
    breaks = c(-1, 0, 1)) +
  scale_y_continuous(
    limits = c(0, 6),
    breaks = c(0, 2, 6)) +
  theme_minimal(base_size = 16) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_blank(),
    legend.position = "none") +
  labs(
    x = expression("fold change ("*log[2]*")"),
    y = expression("P value (- "*log[10]*")"))
# HEART
region = "HR"
DEGS_sub = DEGS %>%
  filter(Region == region) %>%
  mutate(
    significant = pval < 0.01 & abs(diff) >= 0.5,
    diff_capped = pmax(pmin(diff, 3), -3),
    score = -log10(pval) * abs(diff)
  )
top_labels = DEGS_sub %>%
  filter(abs(diff) > 0.5) %>%   # FC threshold only
  group_by(ident.1) %>%
  # top positive side
  filter(diff > 0) %>%
  slice_max(order_by = score, n = 8, with_ties = FALSE) %>%
  ungroup() %>%
  bind_rows(
    DEGS_sub %>%
      filter(abs(diff) > 0.5) %>%   # FC threshold only
      group_by(ident.1) %>%
      # top negative side
      filter(diff < 0) %>%
      slice_max(order_by = score, n = 8, with_ties = FALSE) %>%
      ungroup()
  )
top_labels$GeneSymbol[top_labels$GeneSymbol=="Igkv17-127;Igkv17-121"] = "Igkv17"
top_labels = top_labels %>% filter(pval < 0.01)
ggplot(DEGS_sub, aes(
  x = diff_capped,
  y = -log10(pval),
  color = significant)) +
  geom_point(alpha = 0.8) +
  geom_vline(xintercept = c(-1, 1),
             colour = "grey70",
             linetype = "dotted") +
  geom_hline(yintercept = 2,
             colour = "grey70",
             linetype = "dotted") +
  geom_text_repel(
    data = top_labels,
    aes(label = GeneSymbol),
    color = "black",
    size = 3.5,
    max.overlaps = Inf,
    box.padding = 0.4,
    point.padding = 0.4,
    segment.color = "grey70",
    segment.size = 0.5) +
  facet_wrap(~ident.1) +
  scale_color_manual(values = c("grey80", "#c7522a")) +
  scale_x_continuous(
    limits = c(-3, 3),
    breaks = c(-1, 0, 1)) +
  scale_y_continuous(
    limits = c(0, 6),
    breaks = c(0, 2, 6)) +
  theme_minimal(base_size = 16) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_blank(),
    legend.position = "none") +
  labs(
    x = expression("fold change ("*log[2]*")"),
    y = expression("P value (- "*log[10]*")"))
# ---- Plot volcano plots vs NN ----
DEGS_NN = read.xlsx("./Results/DiffExpr_ProDA_NOISE_NIST_VS_NN.xlsx")
region = "HR"
DEGS_sub = DEGS_NN %>%
  filter(Region == region) %>%
  mutate(
    significant = adj_pval < 0.05 & abs(diff) >= 0.5,
    diff_capped = pmax(pmin(diff, 3), -3),
    score = -log10(pval) * abs(diff)
  )
top_labels = DEGS_sub %>%
  filter(abs(diff) > 0.5, adj_pval < 0.05) %>%   # <- add here
  group_by(ident.1) %>%
  filter(diff > 0) %>%
  slice_max(order_by = score, n = 8, with_ties = FALSE) %>%
  ungroup() %>%
  bind_rows(
    DEGS_sub %>%
      filter(abs(diff) > 0.5, adj_pval < 0.05) %>%  # <- and here
      group_by(ident.1) %>%
      filter(diff < 0) %>%
      slice_max(order_by = score, n = 8, with_ties = FALSE) %>%
      ungroup()
  )
ggplot(DEGS_sub, aes(
  x = diff_capped,
  y = -log10(pval),
  color = significant)) +
  geom_point(alpha = 0.8) +
  geom_vline(xintercept = c(-1, 1),
             colour = "grey70",
             linetype = "dotted") +
  geom_hline(yintercept = 2,
             colour = "grey70",
             linetype = "dotted") +
  geom_text_repel(
    data = top_labels,
    aes(label = GeneSymbol),
    color = "black",
    size = 3.5,
    max.overlaps = Inf,
    box.padding = 0.4,
    point.padding = 0.4,
    segment.color = "grey70",
    segment.size = 0.5) +
  facet_wrap(~ident.1) +
  scale_color_manual(values = c("grey80", "#74a892")) +
  scale_x_continuous(
    limits = c(-3, 3),
    breaks = c(-1, 0, 1)) +
  scale_y_continuous(
    limits = c(0, 10),
    breaks = c(0, 2, 10)) +
  theme_minimal(base_size = 16) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_blank(),
    legend.position = "none") +
  labs(
    x = expression("fold change ("*log[2]*")"),
    y = expression("P value (- "*log[10]*")"))
# SUPPL 1A
ggplot(P_dev_filtered, aes(x = Deviation_clamped, y = GeneSymbol, fill = Broad_Category)) +
  geom_boxplot(outlier.shape = NA) +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.5) +
  geom_hline(data = cat_divs, aes(yintercept = cut_pos), color = "black", linewidth = 0.4) +
  facet_grid(. ~ Condition) +
  coord_cartesian(xlim = c(-4, 4)) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    strip.text.x = element_text(size = 12, face = "bold"),
    axis.text.y = element_text(size = 7),
    panel.spacing = unit(0.5, "lines"),
    panel.grid.major = element_line(color = "grey85", linewidth = 0.3),
    panel.grid.minor = element_blank()
  ) +
  labs(
    x = "Deviation from Control",
    y = "Gene Symbol"
  )
# SUPPL 2A
DAPS = read.xlsx("./Results/DiffExpr_ProDA_Brain_Proteomics_Categorized_FINAL_FULL.xlsx")
DAPS = DAPS %>% filter(!Broad_Category %in% c("Unknown", "Endothelial (Angiogenesis)"))
df = P_spread %>% filter(Region == "BR")
prot_cols = setdiff(colnames(df), c("Sample", "Region", "Condition"))
mat = as.matrix(df[, prot_cols])
mat_norm = proDA::median_normalization(log2(mat))
df_norm = df
df_norm[, prot_cols] = mat_norm
genes_to_plot = DAPS %>%
  distinct(name, GeneSymbol, Broad_Category) %>%
  group_by(GeneSymbol) %>%
  slice(1) %>%
  ungroup()
P_long = df_norm %>%
  pivot_longer(all_of(prot_cols), names_to = "Protein", values_to = "Abundance") %>%
  separate_rows(Protein, sep = ";") %>%
  inner_join(genes_to_plot, by = c("Protein" = "name"))
ctrl_means = P_long %>%
  filter(Condition == "CTRL") %>%
  group_by(GeneSymbol) %>%
  summarise(ctrl_mean = mean(Abundance, na.rm = TRUE), .groups = "drop")
P_dev = P_long %>%
  inner_join(ctrl_means, by = "GeneSymbol") %>%
  mutate(Deviation = Abundance - ctrl_mean)
exposures = unique(DAPS$ident.1)
P_dev_filtered = P_dev %>%
  filter(Condition %in% exposures) %>%
  mutate(Deviation_clamped = pmax(pmin(Deviation, 4), -4))
nn_order = P_dev_filtered %>%
  filter(Condition == "NN") %>%
  group_by(GeneSymbol, Broad_Category) %>%
  summarise(mean_dev_NN = mean(Deviation, na.rm = TRUE), .groups = "drop")
gene_order = nn_order %>%
  arrange(Broad_Category, mean_dev_NN) %>%
  pull(GeneSymbol)
P_dev_filtered = P_dev_filtered %>%
  mutate(GeneSymbol = factor(GeneSymbol, levels = gene_order))
cat_divs = nn_order %>%
  mutate(GeneSymbol = factor(GeneSymbol, levels = gene_order),
         y = as.numeric(GeneSymbol)) %>%
  group_by(Broad_Category) %>%
  summarise(cut_pos = max(y) + 0.5, .groups = "drop")
gene_means = P_dev_filtered %>%
  group_by(Condition, GeneSymbol, Broad_Category) %>%
  summarise(mean_dev = mean(Deviation_clamped, na.rm = TRUE), .groups = "drop")
cat_order = gene_means %>%
  filter(Condition == "NN") %>%
  group_by(Broad_Category) %>%
  summarise(mean_NN = mean(mean_dev, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(mean_NN)) %>%
  pull(Broad_Category)
cat_level = gene_means %>%
  group_by(Condition, Broad_Category) %>%
  summarise(cat_values = list(mean_dev), .groups = "drop") %>%
  unnest(cat_values) %>%
  mutate(Broad_Category = factor(Broad_Category, levels = cat_order))
ggplot(cat_level, aes(y = Broad_Category, x = cat_values, fill = Broad_Category)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  #geom_jitter(height = 0.15, alpha = 0.6, size = 1.8) +   # each dot = gene mean
  geom_vline(xintercept = 0, color = "black", linewidth = 0.4) +
  facet_wrap(~Condition, nrow = 1) +
  coord_cartesian(xlim = c(-2.5, 2.5)) +
  theme_minimal() +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 12, face = "bold"),
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 10)
  ) +
  labs(
    y = "",
    x = "Mean Protein Deviation from CTRL"
  )
ggplot(P_dev_filtered, aes(x = Deviation_clamped, y = GeneSymbol, fill = Broad_Category)) +
  geom_boxplot(outlier.shape = NA) +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.5) +
  geom_hline(data = cat_divs, aes(yintercept = cut_pos), color = "black", linewidth = 0.4) +
  facet_grid(. ~ Condition) +
  coord_cartesian(xlim = c(-4, 4)) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    strip.text.x = element_text(size = 12, face = "bold"),
    axis.text.y = element_text(size = 7),
    panel.spacing = unit(0.5, "lines"),
    panel.grid.major = element_line(color = "grey85", linewidth = 0.3),
    panel.grid.minor = element_blank()
  ) +
  labs(
    x = "Deviation from Control",
    y = "Gene Symbol"
  )
# SUPPL 3A
DAPS = read.xlsx("./Results/DiffExpr_ProDA_Heart_Proteomics_Categorized_FINAL_FULL.xlsx")
df = P_spread %>% filter(Region == "HR")
prot_cols = setdiff(colnames(df), c("Sample", "Region", "Condition"))
mat = as.matrix(df[, prot_cols])
mat_norm = proDA::median_normalization(log2(mat))
df_norm = df
df_norm[, prot_cols] = mat_norm
genes_to_plot = DAPS %>%
  distinct(name, GeneSymbol, Broad_Category) %>%
  group_by(GeneSymbol) %>%
  slice(1) %>%
  ungroup()
P_long = df_norm %>%
  pivot_longer(all_of(prot_cols), names_to = "Protein", values_to = "Abundance") %>%
  separate_rows(Protein, sep = ";") %>%
  inner_join(genes_to_plot, by = c("Protein" = "name"))
ctrl_means = P_long %>%
  filter(Condition == "CTRL") %>%
  group_by(GeneSymbol) %>%
  summarise(ctrl_mean = mean(Abundance, na.rm = TRUE), .groups = "drop")
P_dev = P_long %>%
  inner_join(ctrl_means, by = "GeneSymbol") %>%
  mutate(Deviation = Abundance - ctrl_mean)
exposures = unique(DAPS$ident.1)
P_dev_filtered = P_dev %>%
  filter(Condition %in% exposures) %>%
  mutate(Deviation_clamped = pmax(pmin(Deviation, 4), -4))
nn_order = P_dev_filtered %>%
  filter(Condition == "NN") %>%
  group_by(GeneSymbol, Broad_Category) %>%
  summarise(mean_dev_NN = mean(Deviation, na.rm = TRUE), .groups = "drop")
gene_order = nn_order %>%
  arrange(Broad_Category, mean_dev_NN) %>%
  pull(GeneSymbol)
P_dev_filtered = P_dev_filtered %>%
  mutate(GeneSymbol = factor(GeneSymbol, levels = gene_order))
cat_divs = nn_order %>%
  mutate(GeneSymbol = factor(GeneSymbol, levels = gene_order),
         y = as.numeric(GeneSymbol)) %>%
  group_by(Broad_Category) %>%
  summarise(cut_pos = max(y) + 0.5, .groups = "drop")
gene_means = P_dev_filtered %>%
  group_by(Condition, GeneSymbol, Broad_Category) %>%
  summarise(mean_dev = mean(Deviation_clamped, na.rm = TRUE), .groups = "drop")
cat_order = gene_means %>%
  filter(Condition == "NN") %>%
  group_by(Broad_Category) %>%
  summarise(mean_NN = mean(mean_dev, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(mean_NN)) %>%
  pull(Broad_Category)
cat_level = gene_means %>%
  group_by(Condition, Broad_Category) %>%
  summarise(cat_values = list(mean_dev), .groups = "drop") %>%
  unnest(cat_values) %>%
  mutate(Broad_Category = factor(Broad_Category, levels = cat_order))
ggplot(cat_level, aes(y = Broad_Category, x = cat_values, fill = Broad_Category)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  #geom_jitter(height = 0.15, alpha = 0.6, size = 1.8) +   # each dot = gene mean
  geom_vline(xintercept = 0, color = "black", linewidth = 0.4) +
  facet_wrap(~Condition, nrow = 1) +
  coord_cartesian(xlim = c(-2, 2)) +
  theme_minimal() +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 12, face = "bold"),
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 10)
  ) +
  labs(
    y = "",
    x = "Mean Protein Deviation from CTRL"
  )
ggplot(P_dev_filtered, aes(x = Deviation_clamped, y = GeneSymbol, fill = Broad_Category)) +
  geom_boxplot(outlier.shape = NA) +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.5) +
  geom_hline(data = cat_divs, aes(yintercept = cut_pos), color = "black", linewidth = 0.4) +
  facet_grid(. ~ Condition) +
  coord_cartesian(xlim = c(-4, 4)) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    strip.text.x = element_text(size = 12, face = "bold"),
    axis.text.y = element_text(size = 7),
    panel.spacing = unit(0.5, "lines"),
    panel.grid.major = element_line(color = "grey85", linewidth = 0.3),
    panel.grid.minor = element_blank()
  ) +
  labs(
    x = "Deviation from Control",
    y = "Gene Symbol"
  )
# ---- Enrichment Plots ----
kegg_unif = list.files(paste0(res_folder, "/fGSEA/"), "^KEGG_", full.names = TRUE) %>%
  map_dfr(~ read.xlsx(.x) %>%
            mutate(Region=str_match(basename(.x), "KEGG_(.*)_(.*)\\.xlsx$")[, 2],
                   Condition=str_match(basename(.x), "KEGG_(.*)_(.*)\\.xlsx$")[, 3]))
NES_table = kegg_unif %>%
  filter(padj < 0.25) %>%
  mutate(RegCond = paste(Region, Condition, sep = "_")) %>%
  dplyr::select(pathway, RegCond, NES) %>%
  pivot_wider(names_from = RegCond, values_from = NES) %>%
  mutate(n = rowSums(!is.na(across(-pathway))),
         NES_SUM = rowSums(abs(across(-pathway)), na.rm = TRUE)) %>%
  arrange(desc(n), desc(NES_SUM))
# LUNGS
pathways_not_to_consider = c("Coronavirus disease - COVID-19",
                             "Human papillomavirus infection",
                             "Small cell lung cancer",
                             "Pathways in cancer",
                             "Dilated cardiomyopathy",
                             "Hypertrophic cardiomyopathy",
                             "Amoebiasis",
                             "Cornified envelope formation",
                             "Arrhythmogenic right ventricular cardiomyopathy",
                             "Metabolism of xenobiotics by cytochrome P450",
                             "Tuberculosis",
                             "Cardial muscle contraction",
                             "Ribosome biogenesis in eukaryotes",
                             "Cardiac muscle contraction",
                             "Huntington disease",
                             "Circadian entrainment"
)
df_filtered = kegg_unif %>%
  filter(Region == "LU") %>%
  filter(!pathway %in% pathways_not_to_consider) %>%
  mutate(
    negLogFDR = -log10(padj),   # compute -log10(FDR) for scoring
    negLogFDR = ifelse(is.infinite(negLogFDR), max(negLogFDR[is.finite(negLogFDR)]), negLogFDR)
  )
df_filtered = df_filtered %>%
  group_by(pathway) %>%
  filter(any(padj < 0.25)) %>%
  ungroup()
top_pathways = df_filtered %>%
  group_by(pathway) %>%
  summarize(score = mean(negLogFDR * NES, na.rm = TRUE)) %>%
  arrange(desc(abs(score))) %>%
  slice_head(n = 15) %>%
  pull(pathway)
df_top = df_filtered %>%
  filter(pathway %in% top_pathways) %>%
  mutate(pathway = factor(pathway, levels = top_pathways))   # preserve order
ggplot(df_top, aes(x = Condition, y = pathway, color = NES, size = negLogFDR)) +
  geom_point() +
  scale_color_gradient2(
    low = "#44749D",
    mid = "#EBE7E0",
    high = "#9d6d44",
    midpoint = 0,
    breaks = c(-1, 0, 1)
  ) +
  scale_size_continuous(range = c(3, 7)) +
  theme_minimal(base_size = 16) +
  labs(
    x = "",
    y = "",
    color = "NES",
    size = expression("FDR (- "*log[10]*")")
  ) +
  theme(
    axis.text.y = element_text(size = 10),
    panel.grid.major.y = element_blank()
  )
# HEART
pathways_not_to_consider = c("Coronavirus disease - COVID-19",
                             "Human papillomavirus infection",
                             "Small cell lung cancer",
                             "Pathways in cancer",
                             "Dilated cardiomyopathy",
                             "Hypertrophic cardiomyopathy",
                             "Amoebiasis",
                             "Cornified envelope formation",
                             "Arrhythmogenic right ventricular cardiomyopathy",
                             "Metabolism of xenobiotics by cytochrome P450",
                             "Tuberculosis",
                             "Cardial muscle contraction",
                             "Ribosome biogenesis in eukaryotes",
                             "Cardiac muscle contraction",
                             "Huntington disease",
                             "Circadian entrainment",
                             "Non-alcoholic fatty liver disease",
                             "Thermogenesis",
                             "Salmonella infection",
                             "Serotonergic synapse",
                             "Morphine addiction",
                             "Alzheimer disease",
                             "Parkinson disease"
)
df_filtered = kegg_unif %>%
  filter(Region == "HR") %>%
  filter(!pathway %in% pathways_not_to_consider) %>%
  mutate(
    negLogFDR = -log10(padj),   # compute -log10(FDR) for scoring
    negLogFDR = ifelse(is.infinite(negLogFDR), max(negLogFDR[is.finite(negLogFDR)]), negLogFDR)
  )
df_filtered = df_filtered %>%
  group_by(pathway) %>%
  filter(any(padj < 0.25)) %>%
  ungroup()
top_pathways = df_filtered %>%
  group_by(pathway) %>%
  summarize(score = mean(negLogFDR * NES, na.rm = TRUE)) %>%
  arrange(desc(abs(score))) %>%
  slice_head(n = 10) %>%
  pull(pathway)
df_top = df_filtered %>%
  filter(pathway %in% top_pathways) %>%
  mutate(pathway = factor(pathway, levels = top_pathways))   # preserve order
ggplot(df_top, aes(x = Condition, y = pathway, color = NES, size = negLogFDR)) +
  geom_point() +
  scale_color_gradient2(
    low = "#44749D",
    mid = "#EBE7E0",
    high = "#9d6d44",
    midpoint = 0,
    breaks = c(-1, 0, 1)
  ) +
  scale_size_continuous(range = c(3, 7)) +
  theme_minimal(base_size = 16) +
  labs(
    x = "",
    y = "",
    color = "NES",
    size = expression("FDR (- "*log[10]*")")
  ) +
  theme(
    axis.text.y = element_text(size = 10),
    panel.grid.major.y = element_blank()
  )
