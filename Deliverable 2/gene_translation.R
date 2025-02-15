# This file creates a count data matrix with translated gene names as the row names

library(SummarizedExperiment)
library(tidyverse)
library(biomaRt)

# load data
load(file="readMatrix_PBMC_BulkRNAseq_CVID.Rdata")
european_gene_names <- rownames(mt.counts)

# Extract Ensembl gene names from count matrix
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

# Retrieve gene conversions
gene_conversion <- getBM(
    attributes = c("ensembl_gene_id", "hgnc_symbol"), 
    filters = "ensembl_gene_id", 
    values = european_gene_names, 
    mart = ensembl
)

# Clean translation table by keeping only matching and non-NA entries
gene_conversion <- gene_conversion %>%
    filter(ensembl_gene_id %in% european_gene_names, hgnc_symbol != "") %>%
    drop_na()

# Filter count matrix to retain only valid genes
valid_genes <- row.names(mt.counts) %in% gene_conversion$ensembl_gene_id
mt.counts <- mt.counts[valid_genes, , drop = FALSE]  # Preserve matrix structure

#### Handle Duplicated Genes

# Identify Duplicated Genes
gene_counts <- table(gene_conversion$hgnc_symbol)
duplicated_hgnc_symbols <- names(gene_counts)[gene_counts > 1]
duplicated_genes <- gene_conversion$ensembl_gene_id[gene_conversion$hgnc_symbol %in% duplicated_hgnc_symbols]

# Split duplicated and unduplicated Genes into their own matrices
duplicated_mask <- row.names(mt.counts) %in% duplicated_genes
duplicated_matrix <- mt.counts[duplicated_mask, , drop = FALSE]
unique_matrix <- mt.counts[!duplicated_mask, , drop = FALSE]


# translate unduplicated genes - easy
if(nrow(unique_matrix) > 0) {
    row.names(unique_matrix) <- gene_conversion$hgnc_symbol[
        match(row.names(unique_matrix), gene_conversion$ensembl_gene_id)
    ]
}

# translate unduplicated genes - group_by %>% sum
if(nrow(duplicated_matrix) > 0) {
    # prepare data frame
    gene_df <- as.data.frame(duplicated_matrix)
    hgnc_genes <- gene_conversion$hgnc_symbol[
        match(row.names(duplicated_matrix), gene_conversion$ensembl_gene_id)
    ]
    
    # aggregate gene counts
    gene_df <- gene_df %>%
        mutate(gene = hgnc_genes) %>%
        group_by(gene) %>%
        summarise(across(everything(), sum, na.rm = TRUE), .groups = "drop")
    
    # reformat into matrix, delete gene column
    names <- gene_df$gene
    gene_df$gene <- NULL
    duplicated_matrix <- as.matrix(gene_df)
    rownames(duplicated_matrix) <- names
}

# combine translated results
mt.counts <- rbind(unique_matrix, duplicated_matrix)

# Ensure no NA row names
mt.counts <- mt.counts[
    !is.na(row.names(mt.counts)) & row.names(mt.counts) != "", , drop = FALSE
]

# save matrix
write.csv(mt.counts, "translated_genes.txt")