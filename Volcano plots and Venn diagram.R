library("readxl")
library("dplyr")

setwd("path/to/your/data")

dds1 <- read_excel("SIGNIF_Control_vs_Insulin.xlsx", col_names = TRUE)
dds2 <- read_excel("SIGNIF_Control_vs_Drp-1_CA.xlsx", col_names = TRUE)
dds3 <- read_excel("SIGNIF_Insulin_vs_Drp-1_CA_Insulin.xlsx", col_names = TRUE)
dds4 <- read_excel("SIGNIF_Drp-1_CA_vs_Drp-1_CA_Insulin.xlsx", col_names = TRUE)

f_dds1 <- dds1 %>%
  dplyr::filter(abs(log2FoldChange) > 0.26)

f_dds2 <- dds2 %>%
  dplyr::filter(abs(log2FoldChange) > 0.26)

f_dds3 <- dds3 %>%
  dplyr::filter(abs(log2FoldChange) > 0.26)

f_dds4 <- dds4 %>%
  dplyr::filter(abs(log2FoldChange) > 0.26)



#Volcano Plots------------------------------------------------------------------
keyvals <- ifelse(f_dds2$log2FoldChange > 0, '#009E73', '#56B4E9')
names(keyvals)[keyvals == '#009E73'] <- 'Upregulated'
names(keyvals)[keyvals == '#56B4E9'] <- 'Downregulated'

top_genes <- f_dds2 %>%
  dplyr::arrange(desc(abs(log2FoldChange))) %>%
  dplyr::slice_head(n = 15) %>%
  dplyr::pull(GeneSymbol)



f_dds2$label <- ifelse(f_dds2$GeneSymbol %in% top_genes & !grepl("^LOC", f_dds2$GeneSymbol), 
                       as.character(f_dds2$GeneSymbol), NA)

library("EnhancedVolcano")
EnhancedVolcano(f_dds2,
                lab = f_dds2$label,
                x = 'log2FoldChange',
                y = 'padj',
                title = NULL,
                subtitle = NULL,
                gridlines.major = FALSE,
                gridlines.minor = FALSE,
                colCustom = keyvals,
                legendPosition = "none",
                cutoffLineType = "blank",
                drawConnectors = TRUE,
                labSize = 4,
                caption = NULL,
                axisLabSize = 20
)

#Venn Diagram-------------------------------------------------------------------
library(VennDiagram)

dds1_up <- f_dds1 %>%
  filter(log2FoldChange > 0) %>%
  pull(GeneSymbol)

dds1_down <- f_dds1 %>%
  filter(log2FoldChange < 0) %>%
  pull(GeneSymbol)

dds4_up <- f_dds4 %>%
  filter(log2FoldChange > 0) %>%
  pull(GeneSymbol)

dds4_down <- f_dds4 %>%
  filter(log2FoldChange < 0) %>%
  pull(GeneSymbol)

dds3_up <- f_dds3 %>%
  filter(log2FoldChange > 0) %>%
  pull(GeneSymbol)

dds3_down <- f_dds3 %>%
  filter(log2FoldChange < 0) %>%
  pull(GeneSymbol)

dds2_up <- f_dds2 %>%
  filter(log2FoldChange > 0) %>%
  pull(GeneSymbol)

dds2_down <- f_dds2 %>%
  filter(log2FoldChange < 0) %>%
  pull(GeneSymbol)

circle_edge_colors <- c("#993366", "#004080")  # Darker versions of #CC79A7 and #0072B2

venn_up <- venn.diagram(
  x = list(Control = dds1_up, "Drp-1 CA" = dds4_up),
  filename = NULL,  # Set to NULL to plot directly in R
  fill = c("#CC79A7", "#0072B2"),
  alpha = 0.5,
  cex = 2.5,
  fontfamily = "serif",
  cat.cex = 0,
  cat.pos = 0,
  cat.dist = 0,
  euler.d = FALSE,
  scaled = FALSE,
  resolution = 20000,
  lwd = 2,                      
  col = circle_edge_colors
)
grid.draw(venn_up)

venn_down <- venn.diagram(
  x = list(Control = dds1_down, "Drp-1 CA" = dds4_down),
  filename = NULL,  # Set to NULL to plot directly in R
  fill = c("#CC79A7", "#0072B2"),
  alpha = 0.5,
  cex = 2.5,
  fontfamily = "serif",
  cat.cex = 0,
  cat.pos = 0,
  cat.dist = 0,
  euler.d = FALSE,
  scaled = FALSE,
  resolution = 20000,
  lwd = 3,                      
  col = circle_edge_colors
)
grid.draw(venn_down)

f_dds1 <- f_dds1 %>%
  mutate(GeneID = as.character(GeneID))


# Combined data with both conditions
combined_data <- bind_rows(
  f_dds1 %>%
    filter(GeneSymbol %in% dds1_up | GeneSymbol %in% dds1_down) %>%
    mutate(Condition = "Control"),
  f_dds4 %>%
    filter(GeneSymbol %in% dds4_up | GeneSymbol %in% dds4_down) %>%
    mutate(Condition = "Drp-1 CA")
)

# For unique genes in Drp-1 CA (up-regulated)
unique_to_drp1_up <- setdiff(dds4_up, dds1_up)
unique_to_drp1_up_data <- combined_data %>%
  filter(Condition == "Drp-1 CA" & GeneSymbol %in% unique_to_drp1_up) %>%
  arrange(desc(log2FoldChange))

# Get top 5 differentially expressed genes
top_5_drp1_up <- head(unique_to_drp1_up_data, 10)
print(top_5_drp1_up)

# For unique genes in Control (up-regulated)
unique_to_control_up <- setdiff(dds1_up, dds4_up)
unique_to_control_up_data <- combined_data %>%
  filter(Condition == "Control" & GeneSymbol %in% unique_to_control_up) %>%
  arrange(desc(log2FoldChange))

# Get top 5 differentially expressed genes
top_5_control_up <- head(unique_to_control_up_data, 10)
print(top_5_control_up)

#For unique genes in Drp-1 CA (down-regulated)
unique_to_drp1_down <- setdiff(dds4_down, dds1_down)
unique_to_drp1_down_data <- combined_data %>%
  filter(Condition == "Drp-1 CA" & GeneSymbol %in% unique_to_drp1_down) %>%
  arrange(desc(log2FoldChange))  # Sort by log2FoldChange

# Get top 5 differentially expressed genes
top_5_drp1_down <- head(unique_to_drp1_down_data, 10)
print(top_5_drp1_down)

#For unique genes in Control (down-regulated)
unique_to_control_down <- setdiff(dds1_down, dds4_down)
unique_to_control_down_data <- combined_data %>%
  filter(Condition == "Control" & GeneSymbol %in% unique_to_control_down) %>%
  arrange(desc(log2FoldChange)) 

# Get top 5 differentially expressed genes
top_5_control_down <- head(unique_to_control_down_data, 10)
print(top_5_control_down)

