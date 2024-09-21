library(clusterProfiler)
library(org.Rn.eg.db)
library(DOSE)
library(enrichplot)
library(readxl)
library(dplyr)
library(ggplot2)

setwd("path/to/your/data")

dds1 <- read_excel("SIGNIF_Control_vs_Insulin.xlsx", col_names = TRUE)
dds2 <- read_excel("SIGNIF_Control_vs_Drp-1_CA.xlsx", col_names = TRUE)
dds4 <- read_excel("SIGNIF_Drp-1_CA_vs_Drp-1_CA_Insulin.xlsx", col_names = TRUE)
gene_universe <- readRDS("gene_universe.rds")

f_dds1 <- dds1 %>%
  dplyr::filter(abs(log2FoldChange) > 0.26)

f_dds2 <- dds2 %>%
  dplyr::filter(abs(log2FoldChange) > 0.26)

f_dds4 <- dds4 %>%
  dplyr::filter(abs(log2FoldChange) > 0.26)

dds1_up <- f_dds1 %>%
  filter(log2FoldChange > 0) %>%
  dplyr::pull(GeneID)

dds1_down <- f_dds1 %>%
  filter(log2FoldChange < 0) %>%
  dplyr::pull(GeneID)

dds2_up <- f_dds2 %>%
  filter(log2FoldChange > 0) %>%
  dplyr::pull(GeneID)

dds2_down <- f_dds2 %>%
  filter(log2FoldChange < 0) %>%
  dplyr::pull(GeneID)

dds4_up <- f_dds4 %>%
  filter(log2FoldChange > 0) %>%
  dplyr::pull(GeneID)

dds4_down <- f_dds4 %>%
  filter(log2FoldChange < 0) %>%
  dplyr::pull(GeneID)

options(max.print = .Machine$integer.max, scipen = 999, stringsAsFactors = F, 
        dplyr.summarise.inform = F) 

up_input <- list(GFP = dds1_up, DRP1 = dds4_up)
down_input <- list(GFP = dds1_down, DRP1 = dds4_down)


compare_results_kegg <- compareCluster(
  geneClusters = up_input,
  fun = "enrichKEGG",
  organism = "rno",
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.05
)

compare_results_go <- compareCluster(
  geneClusters = up_input,
  fun = "enrichGO",
  OrgDb = org.Rn.eg.db,
  ont = "BP",
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.05
)

dotplot(compare_results_kegg)

compare_results_go_simplified <- simplify(compare_results_go, cutoff = 0.7, by = "p.adjust", select_fun = min)
dotplot(compare_results_go_simplified)
dotplot(compare_results_go_simplified, showCategory = 10, title = "Enriched Pathways" , 
        split=".sign") + facet_grid(.~.sign)

CR_KEGG <- pairwise_termsim(compare_results_kegg)
CR_GO <- pairwise_termsim(compare_results_go_simplified)
emapplot(CR_KEGG)

dotplot(compare_results_go_simplified)
dotplot(compare_results_go_simplified, showCategory = 5, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)


CR_KEGG_symbols <- setReadable(CR_KEGG, OrgDb = org.Rn.eg.db, keyType = "ENTREZID")
cnetplot(CR_KEGG_symbols)

#-------------------------------------------------------------------------------
kegg_results_2 <- enrichKEGG(
  dds2_down,
  organism = "rno",
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  universe = gene_universe
)

#go_results_2 <- enrichGO(
#  dds2_down,
#  OrgDb = org.Rn.eg.db,
#  ont = "BP",
#  minGSSize = 10,
#  maxGSSize = 500,
#  pvalueCutoff = 0.05, 
#)

significant_kegg_results_up <- kegg_results_2[kegg_results_2@result$pvalue <= 0.05, ]
significant_kegg_results_down <- kegg_results_2[kegg_results_2@result$pvalue <= 0.05, ]

SKRU <- significant_kegg_results_up %>%
  filter(Description %in% c("",
                            "Type II diabetes mellitus",
                            "Insulin secretion")) %>%
  mutate(Regulation = "Upregulation")

SKRD <- significant_kegg_results_down %>%
  filter(Description %in% c("AGE-RAGE signaling pathway in diabetic complications",
                            "Mitophagy - animal",
                            "PI3K-Akt signaling pathway",
                            "Protein processing in endoplasmic reticulum",
                            "")) %>%
  mutate(Regulation = "Downregulation")

combined_kegg_results <- bind_rows(SKRD, SKRU)

combined_kegg_results <- combined_kegg_results %>%
  mutate(Regulation = factor(Regulation, levels = c("Upregulation", "Downregulation")))

combined_kegg_results <- combined_kegg_results %>%
  mutate(GeneRatioNumeric = sapply(strsplit(GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2])))


ggplot(combined_kegg_results, aes(x = Regulation, y = Description)) +
  geom_point(aes(size = GeneRatioNumeric, color = pvalue)) +
  scale_color_gradient(low = "red", high = "blue") +  # Color based on p-value
  theme_minimal() +
  labs(y = "Pathway", color = "p-value", size = "Gene Ratio") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank())


go_results_2_s <- simplify(go_results_2, cutoff = 0.7, by = "p.adjust", select_fun = min)
go_results_s <- simplify(go_results, cutoff = 0.7, by = "p.adjust", select_fun = min)


