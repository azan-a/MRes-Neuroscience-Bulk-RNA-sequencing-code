library(pathview)
library(readxl)
library(dplyr)

setwd("path/to/your/data")
dds4 <- read_excel("SIGNIF_Control_vs_Drp-1_CA.xlsx", col_names = TRUE)
f_dds4 <- dds4 %>%
  filter(abs(log2FoldChange) > 0.26)

gene_data <- dds4$log2FoldChange
names(gene_data) <- dds4$GeneID

pathway_id <- "rno04930"  # Replace with KEGG pathway ID of interest

# Visualize the KEGG pathway
pathview(gene.data = gene_data, pathway.id = pathway_id, species = "rno",
         node.sum = "mean", low = list(gene = "#4A6FE3"), 
         high = list(gene = "#D33F6A"),
         limit = list(gene = c(-1, 1)),
         bins = 20)
head(gene_data)
