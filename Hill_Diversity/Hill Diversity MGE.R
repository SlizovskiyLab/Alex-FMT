library(stringr)
library(microViz)
library(dplyr)
library(ggplot2)
library(BiocManager)
library(phyloseq)
library(ggpubr)
library(vegan)
library(tidyr)
library(hilldiv2)

###MOBILOME###

# Load datasets
metadata <- read.csv("FMT_full_dataset_03_21.csv", header = TRUE)
mge_matrix <- read.csv("Telcomb_MGE_analytical_matrix.csv", header = TRUE, row.names = 1)
mge_classification <- read.csv("MGE_total_classification.csv", header = TRUE, stringsAsFactors = FALSE)

metadata_table_MGE <- metadata[1:148,]

# --- Format metadata for phyloseq ---
metadata_data_frame <- metadata_table_MGE %>%
  select(run_accession, donor_or_recipient, donor_pre_post_timepoint, donor_pre_post, timepoint, timepoint_type, timepoint_group) %>%
  mutate(timepoint = ifelse(donor_or_recipient == "Donor", -7, timepoint))

metadata_sample_data <- sample_data(metadata_data_frame)
rownames(metadata_sample_data) <- metadata_table_MGE$run_accession

# Prepare taxonomy table
colnames(mge_classification)[1:2] <- c("MGE_ID", "Classification")
mge_classification <- tibble::column_to_rownames(mge_classification, var = "MGE_ID")
tax_mat <- as.matrix(mge_classification)
tax_table_obj <- tax_table(tax_mat)

# Prepare OTU table
otu_table_obj <- otu_table(as.matrix(mge_matrix), taxa_are_rows = TRUE)

# Combine phyloseq object
physeq_obj <- phyloseq(otu_table_obj, metadata_sample_data, tax_table_obj)

# Calculate Hill diversity
compute_hill <- function(ps_obj, qs) {
  # 1. Extract OTU matrix
  otu_mat <- as(otu_table(ps_obj), "matrix")
  if (!taxa_are_rows(ps_obj)) otu_mat <- t(otu_mat)
  
  # 2. Compute Hill numbers
  hillmat <- hilldiv(data = otu_mat, q = qs)
  
  # 3. Turn into tibble, *binding* qs as the first column
  hill_df <- tibble(
    q = qs
  ) %>%
    bind_cols(
      as_tibble(hillmat)   # columns are sample names
    ) %>%
    pivot_longer(
      cols      = -q,
      names_to  = "Sample",
      values_to = "Diversity"
    )
  
  return(hill_df)
}

# q range for Hill diversity
q_vals <- seq(0, 2, 0.1)

# Run calculation
hill_results <- compute_hill(physeq_obj, q_vals)


hill_results_with_meta <- hill_results %>%
  left_join(
    metadata_data_frame %>% select(run_accession, donor_pre_post_timepoint),
    by = c("Sample" = "run_accession")
  )

hill_results_with_meta <- hill_results_with_meta %>%
  mutate(
    Timepoint = dplyr::recode(
      donor_pre_post_timepoint,
      "PostFMT"  = "PostFMT 1–30 days",
      "PostFMT2" = "PostFMT 31–60 days",
      "PostFMT3" = "PostFMT 60+ days",
      .default   = donor_pre_post_timepoint  # keeps "PreFMT" and "Donor"
    ),
    Timepoint = factor(
      Timepoint,
      levels = c("PreFMT", "Donor", "PostFMT 1–30 days", "PostFMT 31–60 days", "PostFMT 60+ days")
    )
  )

# Compute mean ± SE by q × Timepoint 
mge_summary <- hill_results_with_meta %>%
  group_by(Timepoint, q) %>%
  summarise(
    mean_div = mean(Diversity, na.rm = TRUE),
    se       = sd(Diversity, na.rm = TRUE)/sqrt(dplyr::n()),
    .groups  = "drop"
  )

# Palette similar to ARG plot
pal_col <- c(
  "PreFMT"             = "#CC79A7",  # magenta
  "Donor"              = "#E69F00",  # orange
  "PostFMT 1–30 days"  = "#F0E442",  # yellow
  "PostFMT 31–60 days" = "#009E73",  # green
  "PostFMT 60+ days"   = "#56B4E9"   # light blue
)



q_values <- c(0, 1, 2)

kw_results <- map_dfr(q_values, function(qv) {
  dfq <- hill_results_with_meta %>% filter(q == qv)
  # Need ≥2 groups to run KW
  if (n_distinct(dfq$Timepoint, na.rm = TRUE) < 2) {
    return(tibble(q = qv, p_value = NA_real_, label = "NA"))
  }
  kw <- kruskal.test(Diversity ~ Timepoint, data = dfq)
  tibble(
    q = qv,
    p_value = kw$p.value,
    label = case_when(
      kw$p.value < 0.001 ~ "***",
      kw$p.value < 0.01  ~ "**",
      kw$p.value < 0.05  ~ "*",
      TRUE               ~ "ns"
    )
  )
})


# Where to place the stars
annot_y <- mge_summary %>%
  filter(q %in% q_values) %>%
  group_by(q) %>%
  summarise(y_pos = max(mean_div + se, na.rm = TRUE) * 1.05, .groups = "drop")

kw_ann <- left_join(kw_results, annot_y, by = "q")

p_mge <- ggplot(mge_summary,
                aes(x = q, y = mean_div,
                    color = Timepoint, shape = Timepoint,
                    linetype = Timepoint, group = Timepoint)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = mean_div - se, ymax = mean_div + se, fill = Timepoint),
              alpha = 0.2, color = NA) +
  scale_x_continuous(breaks = c(0, 1, 2)) +
  scale_color_manual(values = pal_col) +
  scale_fill_manual(values = pal_col) +
  labs(
    x = "Hill order (q)",
    y = "Mean Effective number of features",
    color = "Timepoint", fill = "Timepoint",
    shape = "Timepoint", linetype = "Timepoint"
  ) +
  theme_minimal(base_size = 14) +
  theme(panel.grid.minor = element_blank()) +
  # ← asterisks
  geom_text(data = kw_ann,
            aes(x = q, y = y_pos, label = label),
            inherit.aes = FALSE, size = 4, fontface = "bold", nudge_x = 0.05)

p_mge





#### PUTTING IT ALL TOGETHER ####

library(tidyverse)
library(phyloseq)
library(hilldiv2)

# Load and prepare metadata
metadata <- read.csv("FMT_full_dataset_03_21.csv", header = TRUE)
mge_matrix <- read.csv("Telcomb_MGE_analytical_matrix.csv", header = TRUE, row.names = 1)
mge_classification <- read.csv("MGE_total_classification.csv", header = TRUE, stringsAsFactors = FALSE)


metadata_df <- metadata[1:148, ] %>%
  select(run_accession, donor_or_recipient, donor_pre_post_timepoint,
         donor_pre_post, timepoint, timepoint_type, timepoint_group) %>%
  mutate(timepoint = ifelse(donor_or_recipient == "Donor", -7, timepoint))

# clean sample_data for joining
meta_keep <- metadata_df %>%
  select(run_accession, donor_pre_post_timepoint)

# classification table (IDs must match rownames in mge_matrix)
classification_tbl <- mge_classification %>%
  select(IDs, final_classification) %>%
  column_to_rownames("IDs")

## ---------- hill function ----------
compute_hill_tbl <- function(mat, qs) {
  mat <- as.matrix(mat)
  storage.mode(mat) <- "double"
  
  # drop all-zero features
  if (nrow(mat) > 0) mat <- mat[rowSums(mat) > 0, , drop = FALSE]
  
  # drop samples with too few counts (and then TSS-normalize)
  keep <- colSums(mat) > 0
  if (any(!keep)) mat <- mat[, keep, drop = FALSE]
  if (ncol(mat) == 0L || nrow(mat) == 0L) {
    return(tibble(q = numeric(0), Sample = character(0), Diversity = numeric(0)))
  }
  mat <- sweep(mat, 2, colSums(mat), "/")  # TSS to relative abundances
  
  hm <- hilldiv2::hilldiv(mat, q = qs) |> as.matrix()
  colnames(hm) <- colnames(mat)
  
  tibble(
    q         = rep(qs, times = ncol(hm)),
    Sample    = rep(colnames(hm), each  = length(qs)),
    Diversity = as.vector(hm)
  ) %>% mutate(Diversity = if_else(is.finite(Diversity), Diversity, NA_real_))
}

## ---------- per-MGE type hill (returns long table) ----------
compute_hill_with_label <- function(mge_type, q_vals = seq(0, 2, 0.1)) {
  mge_ids <- rownames(classification_tbl)[classification_tbl$final_classification == mge_type]
  if (length(mge_ids) == 0) return(NULL)
  
  sub_mat <- mge_matrix[rownames(mge_matrix) %in% mge_ids, , drop = FALSE]
  if (nrow(sub_mat) == 0) return(NULL)
  
  hill_long <- compute_hill_tbl(sub_mat, q_vals) %>%
    left_join(meta_keep, by = c("Sample" = "run_accession")) %>%
    mutate(MGE_type = mge_type)
  
  hill_long
}

## ---------- run for all MGE types ----------
mge_types <- c("plasmid", "ICE", "prophage", "virus", "likely IS/TE")
q_vals <- seq(0, 2, 0.1)

hill_all_raw <- purrr::map_dfr(mge_types, compute_hill_with_label)

# Tidy labels and 5-level timepoint like ARG plot
hill_all <- hill_all_raw %>%
  # keep only the 5 levels you plot; drop missing labels
  dplyr::filter(!is.na(donor_pre_post_timepoint)) %>%
  dplyr::mutate(
    Timepoint = dplyr::recode(
      donor_pre_post_timepoint,
      "PostFMT"  = "PostFMT 1–30 days",
      "PostFMT2" = "PostFMT 31–60 days",
      "PostFMT3" = "PostFMT 60+ days",
      .default   = donor_pre_post_timepoint
    ),
    Timepoint = factor(Timepoint,
                       levels = c("PreFMT","Donor","PostFMT 1–30 days","PostFMT 31–60 days","PostFMT 60+ days")
    )
  ) %>%
  # drop non-finite values just in case
  dplyr::filter(is.finite(Diversity))

# Rename MGE type Labels
hill_all <- hill_all %>%
  mutate(MGE_type = recode(MGE_type,
                           "plasmid" = "Plasmid",
                           "ICE" = "ICE",
                           "prophage" = "Prophage",
                           "virus" = "Virus",
                           "likely IS/TE" = "IS/TE"
  ))

## ---------- summary (mean ± SE) for ribbons/lines ----------
mge_summary <- hill_all %>%
  group_by(MGE_type, Timepoint, q) %>%
  summarise(
    mean_div = mean(Diversity, na.rm = TRUE),
    se       = sd(Diversity,   na.rm = TRUE) / sqrt(n()),
    .groups  = "drop"
  )

## ---------- Kruskal–Wallis stars at q = 0,1,2 per facet ----------
q_values <- c(0,1,2)

kw_ann <- hill_all %>%
  filter(q %in% q_values) %>%
  group_by(MGE_type, q) %>%
  summarise(
    p_value = tryCatch(kruskal.test(Diversity ~ Timepoint)$p.value, error = function(e) NA_real_),
    .groups = "drop"
  ) %>%
  mutate(label = case_when(
    p_value < 0.001 ~ "***",
    p_value < 0.01  ~ "**",
    p_value < 0.05  ~ "*",
    TRUE            ~ "ns"
  )) %>%
  left_join(
    mge_summary %>%
      filter(q %in% q_values) %>%
      group_by(MGE_type, q) %>%
      summarise(y_pos = max(mean_div + se, na.rm = TRUE) * 1.05, .groups = "drop"),
    by = c("MGE_type","q")
  )

## ---------- palettes to match ARG style ----------
pal_col <- c(
  "PreFMT"             = "#CC79A7",
  "Donor"              = "#E69F00",
  "PostFMT 1–30 days"  = "#F0E442",
  "PostFMT 31–60 days" = "#009E73",
  "PostFMT 60+ days"   = "#56B4E9"
)

## ---------- plot ----------
p_mge_byType <- ggplot(mge_summary,
       aes(x = q, y = mean_div,
           color = Timepoint, shape = Timepoint,
           linetype = Timepoint, group = Timepoint)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = mean_div - se, ymax = mean_div + se, fill = Timepoint),
              alpha = 0.2, color = NA, show.legend = FALSE) +   # <- hide fill legend
  facet_wrap(~ MGE_type, scales = "free_y") +
  scale_x_continuous(breaks = c(0,1,2)) +
  scale_color_manual(values = pal_col) +
  scale_shape_manual(values = c(16, 16, 16, 16, 16)) +          # optional: all circles
  scale_linetype_manual(values = c("solid","dashed","dotted","dotdash","longdash")) +
  labs(
    x = "Hill order (q)",
    y = "Mean Effective number of features",
    color = "Timepoint", shape = "Timepoint", linetype = "Timepoint"
  ) +
  theme_minimal(base_size = 13) +
  theme(panel.grid.minor = element_blank()) +
  geom_text(data = kw_ann,
            aes(x = q, y = y_pos, label = label),
            inherit.aes = FALSE, size = 4, fontface = "bold", nudge_x = 0.05)

p_mge_byType