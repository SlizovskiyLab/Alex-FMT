library(stringr)
library(microViz)
library(ggplot2)
library(phyloseq)
library(tidyverse)
library(hilldiv2)
library(purrr)

conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")

# Load previously saved phyloseq objects
physeq_metadata_rCDI   <- readRDS("physeq_metadata_rCDI.rds")
physeq_metadata_Melanoma <- readRDS("physeq_metadata_Melanoma.rds")
physeq_metadata_MDRB   <- readRDS("physeq_metadata_MDRB.rds")

# ---- 1) One robust helper (defined once) ----
compute_hill <- function(ps_obj, qs) {
  m <- as(otu_table(ps_obj), "matrix")
  if (!taxa_are_rows(ps_obj)) m <- t(m)
  
  hm <- hilldiv(data = m, q = qs)
  hm <- as.matrix(hm)
  # keep sample names
  colnames(hm) <- colnames(m)
  
  tibble(q = qs) |>
    bind_cols(as_tibble(hm, .name_repair = "minimal")) |>
    pivot_longer(-q, names_to = "Sample", values_to = "Diversity")
}

# ---- 2) Full pipeline for one phyloseq object ----
hill_plot_for_phyloseq <- function(ps, cohort_label,
                                   qs = seq(0, 2, by = 0.1),
                                   prefer_time_col = c("donor_pre_post_timepoint","donor_pre_post")) {
  
  # Aggregate at your three ranks
  glom_list <- list(
    Class     = tax_glom(ps, "class"),
    Mechanism = tax_glom(ps, "mechanism"),
    Group     = tax_glom(ps, "group")
  )
  
  # TSS normalize
  glom_list_norm <- imap(glom_list, ~ transform_sample_counts(.x, function(z) z / sum(z)))
  
  # Hill numbers for each rank
  hill_all <- imap_dfr(glom_list_norm, ~ compute_hill(.x, qs) %>% mutate(FeatureType = .y))
  
  # ---- metadata from the SAME object ----
  meta_df <- as(sample_data(glom_list_norm[[1]]), "data.frame") |>
    rownames_to_column("Sample")
  
  # Clean the weird suffixes and remove dup names
  names(meta_df) <- sub("\\.\\.\\.\\..*$", "", names(meta_df))
  meta_df <- meta_df[, !duplicated(names(meta_df))]
  
  # Find a usable time column
  have_cols <- intersect(prefer_time_col, names(meta_df))
  if (length(have_cols) == 0) {
    stop("No timepoint-like column found. Available columns: ",
         paste(names(meta_df), collapse = ", "))
  }
  time_col <- have_cols[1]
  
  # Build a clean Timepoint
  meta_df <- meta_df |>
    mutate(Timepoint = .data[[time_col]]) |>
    mutate(Timepoint = dplyr::recode(
      Timepoint,
      "PostFMT"  = "PostFMT 1–30 days",
      "PostFMT2" = "PostFMT 31–60 days",
      "PostFMT3" = "PostFMT 61+ days",
      .default = Timepoint
    )) |>
    select(Sample, Timepoint)
  
  # Join hill numbers with metadata
  hill_plot_df <- hill_all |>
    left_join(meta_df, by = "Sample")
  
  # quick sanity check
  # stopifnot(all(c("q","Sample","Diversity","FeatureType","Timepoint") %in% names(hill_plot_df)))
  
  # ---- Summary (mean ± SE) ----
  hill_summary <- hill_plot_df |>
    group_by(FeatureType, Timepoint, q) |>
    summarise(
      mean_div = mean(Diversity, na.rm = TRUE),
      se       = sd(Diversity,  na.rm = TRUE) / sqrt(dplyr::n()),
      .groups  = "drop"
    )
  
  # ---- Kruskal–Wallis at q = 0,1,2 ----
  q_values <- c(0, 1, 2)
  feature_types <- c("Class", "Group", "Mechanism")
  
  kw_results_all <- expand.grid(q = q_values, FeatureType = feature_types) |>
    pmap_dfr(function(q, FeatureType) {
      df_q <- hill_plot_df |>
        filter(q == !!q, FeatureType == !!FeatureType)
      
      # Need at least 2 groups and a Diversity column
      if (!all(c("Diversity","Timepoint") %in% names(df_q)) ||
          nrow(df_q) == 0 ||
          df_q$Timepoint |> n_distinct(na.rm = TRUE) < 2) {
        return(tibble(q = q, FeatureType = FeatureType, p_value = NA_real_, label = "NA"))
      }
      
      # Build formula: Diversity ~ Timepoint
      fml <- reformulate("Timepoint", response = "Diversity")
      test <- kruskal.test(fml, data = df_q)
      
      tibble(
        q = q,
        FeatureType = FeatureType,
        p_value = test$p.value,
        label = case_when(
          test$p.value < 0.001 ~ "***",
          test$p.value < 0.01  ~ "**",
          test$p.value < 0.05  ~ "*",
          TRUE                 ~ "ns"
        )
      )
    })
  
  # y-position for annotation stars
  max_y <- hill_summary |>
    filter(q %in% q_values) |>
    group_by(q, FeatureType) |>
    summarise(y_pos = max(mean_div + se, na.rm = TRUE) * 1.05, .groups = "drop")
  
  kw_results_all <- left_join(kw_results_all, max_y, by = c("q","FeatureType"))
  
  # ---- Plot ----
  ggplot(hill_summary,
         aes(x = q, y = mean_div,
             color = Timepoint, shape = Timepoint, linetype = Timepoint, group = Timepoint)) +
    geom_line(linewidth = 1) +
    geom_ribbon(aes(ymin = mean_div - se, ymax = mean_div + se, fill = Timepoint),
                alpha = 0.2, color = NA) +
    facet_wrap(~FeatureType, scales = "free_y") +
    scale_x_continuous(breaks = c(0, 1, 2)) +
    labs(
      title = cohort_label,
      x = "Hill order (q)",
      y = "Mean Hill diversity (Effective number of features)",
      color = "Timepoint", fill = "Timepoint", shape = "Timepoint", linetype = "Timepoint"
    ) +
    theme_minimal() +
    geom_text(
      data = kw_results_all,
      aes(x = q, y = y_pos, label = label),
      inherit.aes = FALSE, size = 3, fontface = "bold", nudge_x = 0.1
    )
}

# ---- 3) Run once per cohort ----
res_rCDI     <- hill_plot_for_phyloseq(physeq_metadata_rCDI,     "rCDI")
res_Melanoma <- hill_plot_for_phyloseq(physeq_metadata_Melanoma, "Melanoma")
res_MDRB     <- hill_plot_for_phyloseq(physeq_metadata_MDRB,     "MDRB")


res_rCDI
res_Melanoma
res_MDRB
