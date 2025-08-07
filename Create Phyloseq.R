library(stringr)
library(microViz)
library(dplyr)
library(ggplot2)
library(BiocManager)
library(phyloseq)
library(ggpubr)

# Load data
full_dataset_metadata_03_21 <- read.csv("FMT_full_dataset_03_21.csv")
dedup_AMR_analytic_matrix <- read.csv("dedup_AMR_analytic_matrix.csv")
megares_annotations_v3_00 <- read.csv("megares_annotations_v3.00.csv")


create_phyloseq_object <- function(
    SRR_list, metadata_rows, dedup_matrix, metadata_df, annotation_df
) {
  library(dplyr)
  library(stringr)
  library(phyloseq)
  
  SRR_vector <- unlist(SRR_list, use.names=FALSE)
  
  # Metadata
  metadata_table <- metadata_df[metadata_rows, ]
  metadata_data_frame <- data.frame(
    timepoint = metadata_table$timepoint,
    timepoint_type = metadata_table$timepoint_type,
    timepoint_group = metadata_table$timepoint_group,
    sequencer = metadata_table$sequencer,
    lib_prep = metadata_table$library_preparation,
    fmt_route = metadata_table$fmt_route,
    fmt_prep = metadata_table$fmt_prep,
    bowel_preparation = metadata_table$bowel_preparation,
    abx_use = metadata_table$abx_use,
    donor_pre_post_timepoint = metadata_table$donor_pre_post_timepoint,
    ID = metadata_table$ID,
    donor_recipient = metadata_table$donor_or_recipient,
    donor_pre_post = metadata_table$donor_pre_post,
    row.names = metadata_table$run_accession
  )
  metadata_data_frame <- metadata_data_frame %>%
    mutate(timepoint = ifelse(donor_recipient == "Donor", -7, timepoint))
  metadata_sample_data <- sample_data(metadata_data_frame)
  
  # OTU Table
  sample_count_import <- dedup_matrix[
    , colnames(dedup_matrix) %in% SRR_vector | colnames(dedup_matrix) == "gene_accession"
  ]
  sample_count_matrix <- data.matrix(sample_count_import[, -1])
  sample_otu_table <- otu_table(sample_count_matrix, taxa_are_rows = TRUE)
  rownames(sample_otu_table) <- sample_count_import[, 1]
  
  # Prune and remove SNPs
  pruned_initial <- prune_taxa(taxa_sums(sample_otu_table) > 0, sample_otu_table)
  for (k in 1:nrow(pruned_initial)) {
    if (str_sub(rownames(pruned_initial)[k], -23, -1) == "RequiresSNPConfirmation") {
      pruned_initial[k, ] <- 0
    }
  }
  pruned_no_snp <- prune_taxa(taxa_sums(pruned_initial) > 0, pruned_initial)
  
  # Tax Table
  ARG_table <- tax_table(as.matrix(annotation_df[, -1]))
  colnames(ARG_table)[1:4] <- colnames(annotation_df)[2:5]
  rownames(ARG_table) <- annotation_df[, 1]
  
  # Create phyloseq
  physeq <- phyloseq(pruned_no_snp, ARG_table)
  physeq <- ps_reorder(physeq, SRR_vector)
  physeq_merged <- merge_phyloseq(physeq, metadata_sample_data)
  
  return(physeq_merged)
}

# Create phyloseq objects
physeq_metadata_rCDI <- create_phyloseq_object(SRR_list_R, 1:148,
                                               dedup_AMR_analytic_matrix, full_dataset_metadata_03_21, megares_annotations_v3_00)

physeq_metadata_Melanoma <- create_phyloseq_object(SRR_list_M, 149:203,
                                                   dedup_AMR_analytic_matrix, full_dataset_metadata_03_21, megares_annotations_v3_00)

physeq_metadata_MDRB <- create_phyloseq_object(SRR_list_MD, 204:263,
                                               dedup_AMR_analytic_matrix, full_dataset_metadata_03_21, megares_annotations_v3_00)


# Save RDS files
saveRDS(physeq_metadata_rCDI, "physeq_metadata_rCDI.rds")
saveRDS(physeq_metadata_Melanoma, "physeq_metadata_Melanoma.rds")
saveRDS(physeq_metadata_MDRB, "physeq_metadata_MDRB.rds")