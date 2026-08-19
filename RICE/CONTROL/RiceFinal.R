
############################################################
# Differential Expression Analysis of Rice RNA-seq
# Part 1: Data Import, DESeq2 Analysis and Result Export
############################################################

############################################################
# Set Working Directory
############################################################

setwd("/Users/animikhaghosh/Desktop/OL-RI/WD/PRJNA608550-Rice")

############################################################
# Load Required Library
############################################################

library(DESeq2)

############################################################
# Read Count Matrix
############################################################

countData <- read.delim(
  "counts-rice-copy.txt",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

############################################################
# Keep Only Rice Gene IDs
############################################################

countData <- countData[
  grep("^Os", countData$Geneid),
]

############################################################
# Set Gene IDs as Row Names
############################################################

rownames(countData) <- countData$Geneid

############################################################
# Keep Count Columns Only
############################################################

countData <- countData[, -(1:6)]

############################################################
# Sample Information
############################################################

condition <- factor(c(
  "Control",
  "Control",
  "Control",
  "Control",
  "Treatment",
  "Treatment",
  "Treatment",
  "Treatment"
))

colData <- data.frame(
  row.names = colnames(countData),
  condition = condition
)

############################################################
# Create DESeq2 Object
############################################################

dds <- DESeqDataSetFromMatrix(
  countData = round(as.matrix(countData)),
  colData = colData,
  design = ~ condition
)

############################################################
# Filter Lowly Expressed Genes
############################################################

keep <- rowSums(counts(dds)) >= 10

dds <- dds[keep, ]

############################################################
# Run DESeq2
############################################################

dds <- DESeq(dds)

############################################################
# Extract Results
############################################################

res <- results(dds)

############################################################
# Order by Adjusted p-value
############################################################

resOrdered <- res[
  order(res$padj),
]

############################################################
# Export Complete Results
############################################################

write.csv(
  as.data.frame(resOrdered),
  "DESeq2_all_results.csv"
)

############################################################
# Significant DEGs
############################################################

sig <- subset(
  as.data.frame(resOrdered),
  padj < 0.05 &
    abs(log2FoldChange) > 1
)

write.csv(
  sig,
  "DESeq2_significant_DEGs.csv",
  row.names = TRUE
)

############################################################
# Upregulated Genes
############################################################

upregulated <- subset(
  as.data.frame(resOrdered),
  padj < 0.05 &
    log2FoldChange > 1
)

write.csv(
  upregulated,
  "DESeq2_upregulated_genes.csv",
  row.names = TRUE
)

############################################################
# Downregulated Genes
############################################################

downregulated <- subset(
  as.data.frame(resOrdered),
  padj < 0.05 &
    log2FoldChange < -1
)

write.csv(
  downregulated,
  "DESeq2_downregulated_genes.csv",
  row.names = TRUE
)

############################################################
# Part 2: Annotation, Annotated Results, PCA and MA Plot
############################################################

############################################################
# Load Required Libraries
############################################################

library(DESeq2)
library(ggplot2)
library(ggrepel)
library(pheatmap)

############################################################
# Read Annotation File
############################################################

annot <- read.csv(
  "rice-annotation.txt",
  stringsAsFactors = FALSE
)

############################################################
# Rename Columns
############################################################

colnames(annot) <- c(
  "GeneID",
  "Description",
  "GeneSymbol",
  "Synonym"
)

############################################################
# Convert DESeq2 Results to Data Frame
############################################################

res_df <- as.data.frame(resOrdered)

res_df$GeneID <- rownames(res_df)

############################################################
# Merge Annotation
############################################################

library(dplyr)

annot_unique <- annot %>%
  distinct(GeneID, .keep_all = TRUE)

res_annotated <- merge(
  res_df,
  annot_unique,
  by = "GeneID",
  all.x = TRUE
)
############################################################
# Replace Missing Gene Symbols with Gene IDs
############################################################

res_annotated$GeneSymbol <- ifelse(
  is.na(res_annotated$GeneSymbol) |
    res_annotated$GeneSymbol == "",
  res_annotated$GeneID,
  res_annotated$GeneSymbol
)

############################################################
# Export Annotated Results
############################################################

write.csv(
  res_annotated,
  "DESeq2_all_results_annotated.csv",
  row.names = FALSE
)

############################################################
# Export Significant DEGs
############################################################

sig_annotated <- subset(
  res_annotated,
  padj < 0.05 &
    abs(log2FoldChange) > 1
)

write.csv(
  sig_annotated,
  "DESeq2_significant_DEGs_annotated.csv",
  row.names = FALSE
)

############################################################
# Export Upregulated Genes
############################################################

upregulated_annotated <- subset(
  res_annotated,
  padj < 0.05 &
    log2FoldChange > 1
)

write.csv(
  upregulated_annotated,
  "DESeq2_upregulated_genes_annotated.csv",
  row.names = FALSE
)

############################################################
# Export Downregulated Genes
############################################################

downregulated_annotated <- subset(
  res_annotated,
  padj < 0.05 &
    log2FoldChange < -1
)

write.csv(
  downregulated_annotated,
  "DESeq2_downregulated_genes_annotated.csv",
  row.names = FALSE
)

############################################################
# Variance Stabilizing Transformation
############################################################

vsd <- vst(
  dds,
  blind = FALSE
)

############################################################
# PCA Plot
############################################################

png(
  "PCA_plot.png",
  width = 3000,
  height = 2400,
  res = 300
)

plotPCA(
  vsd,
  intgroup = "condition"
)

dev.off()

############################################################
# MA Plot
############################################################

png(
  "MA_plot_Rice.png",
  width = 3000,
  height = 2400,
  res = 300
)

plotMA(
  resOrdered,
  ylim = c(-6, 6)
)

dev.off()

############################################################
# Part 3: Volcano Plot and Heatmap (Gene Symbols)
############################################################

############################################################
# Part 3: Volcano Plot (Gene Symbols)
############################################################

library(ggplot2)
library(ggrepel)

############################################################
# Prepare Volcano Data
############################################################

volcano_data <- subset(
  res_annotated,
  !is.na(padj)
)

############################################################
# Define Significance
############################################################

volcano_data$Significant <- "Not Significant"

volcano_data$Significant[
  volcano_data$padj < 0.05 &
    volcano_data$log2FoldChange > 1
] <- "Upregulated"

volcano_data$Significant[
  volcano_data$padj < 0.05 &
    volcano_data$log2FoldChange < -1
] <- "Downregulated"

############################################################
# Label for Plot
############################################################

volcano_data$Label <- ifelse(
  is.na(volcano_data$GeneSymbol) |
    volcano_data$GeneSymbol == "",
  volcano_data$GeneID,
  volcano_data$GeneSymbol
)

############################################################
# Top 10 Upregulated Genes
############################################################

top_up <- subset(
  volcano_data,
  padj < 0.05 &
    log2FoldChange > 1
)
top_up <- top_up[
  order(-top_up$log2FoldChange),
]




top_up <- head(top_up, 10)

############################################################
# Top 10 Downregulated Genes
############################################################

top_down <- subset(
  volcano_data,
  padj < 0.05 &
    log2FoldChange < -1
)

top_down <- top_down[
  order(top_down$log2FoldChange),
]



top_down <- head(top_down, 10)

############################################################
# Combine Top Genes
############################################################

top_genes <- rbind(
  top_up,
  top_down
)

############################################################
# Volcano Plot
############################################################

png(
  "Volcano_Plot_Rice_GeneSymbols.png",
  width = 3000,
  height = 2400,
  res = 300
)

ggplot(
  volcano_data,
  aes(
    x = log2FoldChange,
    y = -log10(padj),
    color = Significant
  )
) +
  
  geom_point(
    alpha = 0.7,
    size = 2
  ) +
  
  geom_vline(
    xintercept = c(-1, 1),
    linetype = "dashed",
    linewidth = 0.5
  ) +
  
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed",
    linewidth = 0.5
  ) +
  
  geom_text_repel(
    data = top_genes,
    aes(label = Label),
    size = 4,
    color = "black",
    segment.color = "black",
    box.padding = 0.5,
    point.padding = 0.3,
    max.overlaps = Inf
  ) +
  
  scale_color_manual(
    values = c(
      "Upregulated" = "red",
      "Downregulated" = "blue",
      "Not Significant" = "grey70"
    )
  ) +
  
  theme_bw(base_size = 14) +
  
  labs(
    title = "Rice Differential Expression",
    x = expression(log[2]~Fold~Change),
    y = expression(-log[10]~Adjusted~italic(P))
  )

dev.off()

##### Heatmap


############################################################
# Significant genes (FDR < 0.05)
############################################################

sig_heat <- subset(
  res_annotated,
  padj < 0.05 &
    abs(log2FoldChange) > 1
)

############################################################
# Top 50 Upregulated by Log2FC
############################################################

top_up <- sig_heat[
  order(-sig_heat$log2FoldChange),
]

top_up <- head(top_up, 50)

############################################################
# Top 50 Downregulated by Log2FC
############################################################

top_down <- sig_heat[
  order(sig_heat$log2FoldChange),
]

top_down <- head(top_down, 50)

############################################################
# Combine
############################################################

top100 <- c(
  top_up$GeneID,
  top_down$GeneID
)



############################################################
# Create Gene Symbol Mapping
############################################################

symbol_map <- annot$GeneSymbol
names(symbol_map) <- annot$GeneID

############################################################
# Gene Labels
############################################################

heatmap_labels <- symbol_map[top100]

heatmap_labels[
  is.na(heatmap_labels) |
    heatmap_labels == ""
] <- top100[
  is.na(heatmap_labels) |
    heatmap_labels == ""
]

############################################################
# Extract Expression Matrix
############################################################

mat <- assay(vsd)[top100, ]

############################################################
# Replace Row Names with Gene Symbols
############################################################

rownames(mat) <- make.unique(heatmap_labels)

############################################################
# Scale by Row (Z-score)
############################################################

mat_scaled <- t(
  scale(
    t(mat)
  )
)

############################################################
# Sample Annotation
############################################################

annotation_col <- data.frame(
  Condition = colData(vsd)$condition
)

rownames(annotation_col) <- colnames(mat_scaled)

############################################################
# Heatmap
############################################################

png(
  "Heatmap_Top50_Up_Top50_Down_GeneSymbols.png",
  width = 3200,
  height = 4200,
  res = 300
)

pheatmap(
  mat_scaled,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  annotation_col = annotation_col,
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 5,
  fontsize_col = 10,
  main = "Top 50 Upregulated + Top 50 Downregulated Genes"
)

dev.off()




############################################################
# PART 4: GO + KEGG ENRICHMENT
############################################################

library(gprofiler2)

sig_genes <- unique(sig_annotated$GeneID)

gost_res <- gost(
  query = sig_genes,
  organism = "osativa",
  correction_method = "fdr",
  significant = TRUE,
  sources = c(
    "GO:BP",
    "GO:MF",
    "GO:CC",
    "KEGG"
  )
)


############################################################
# PART 5: EXPORT ENRICHMENT RESULTS
############################################################

enrichment <- gost_res$result

# Remove list columns that cannot be written to CSV
enrichment_export <- enrichment[
  , !sapply(enrichment, is.list)
]

write.csv(
  enrichment_export,
  "Functional_Enrichment_GO_KEGG.csv",
  row.names = FALSE
)









############################################################
# PART 6: FUNCTIONAL ENRICHMENT PLOTS
############################################################

library(ggplot2)

plot_enrichment <- function(data, source_name, file_name, plot_title, fill_color){
  
  df <- subset(data, source == source_name)
  
  if(nrow(df) == 0){
    message("No ", source_name, " terms found.")
    return(NULL)
  }
  
  df <- df[order(df$p_value), ]
  df <- head(df, 20)
  
  p <- ggplot(
    df,
    aes(
      x = reorder(term_name, -log10(p_value)),
      y = -log10(p_value)
    )
  ) +
    geom_col(fill = fill_color) +
    coord_flip() +
    theme_minimal(base_size = 14) +
    labs(
      title = plot_title,
      x = "",
      y = expression(-log[10](FDR))
    )
  
  png(
    filename = file_name,
    width = 3000,
    height = 2400,
    res = 300
  )
  
  print(p)
  
  dev.off()
}

############################################################
# GO Biological Process
############################################################

plot_enrichment(
  enrichment_export,
  "GO:BP",
  "GO_BP_Barplot.png",
  "GO Biological Process",
  "steelblue"
)

############################################################
# GO Molecular Function
############################################################

plot_enrichment(
  enrichment_export,
  "GO:MF",
  "GO_MF_Barplot.png",
  "GO Molecular Function",
  "forestgreen"
)

############################################################
# GO Cellular Component
############################################################

plot_enrichment(
  enrichment_export,
  "GO:CC",
  "GO_CC_Barplot.png",
  "GO Cellular Component",
  "darkorange"
)

############################################################
# KEGG Pathways
############################################################

plot_enrichment(
  enrichment_export,
  "KEGG",
  "KEGG_Barplot.png",
  "KEGG Pathway Enrichment",
  "firebrick"
)

############################################################
# GO + KEGG ENRICHMENT : UPREGULATED GENES
############################################################

up_genes <- unique(upregulated_annotated$GeneID)

gost_up <- gost(
  query = up_genes,
  organism = "osativa",
  correction_method = "fdr",
  significant = TRUE,
  sources = c("GO:BP", "GO:MF", "GO:CC", "KEGG")
)

up_enrichment <- gost_up$result

up_enrichment_export <- up_enrichment[
  , !sapply(up_enrichment, is.list)
]

write.csv(
  up_enrichment_export,
  "Upregulated_GO_KEGG_Enrichment.csv",
  row.names = FALSE
)

############################################################
# GO + KEGG ENRICHMENT : DOWNREGULATED GENES
############################################################

down_genes <- unique(downregulated_annotated$GeneID)

gost_down <- gost(
  query = down_genes,
  organism = "osativa",
  correction_method = "fdr",
  significant = TRUE,
  sources = c("GO:BP", "GO:MF", "GO:CC", "KEGG")
)

down_enrichment <- gost_down$result

down_enrichment_export <- down_enrichment[
  , !sapply(down_enrichment, is.list)
]

write.csv(
  down_enrichment_export,
  "Downregulated_GO_KEGG_Enrichment.csv",
  row.names = FALSE
)

plot_enrichment(enrichment_export, "GO:BP", "GO_BP_All.png", "GO Biological Process (All DEGs)", "steelblue")
plot_enrichment(enrichment_export, "GO:MF", "GO_MF_All.png", "GO Molecular Function (All DEGs)", "forestgreen")
plot_enrichment(enrichment_export, "GO:CC", "GO_CC_All.png", "GO Cellular Component (All DEGs)", "darkorange")
plot_enrichment(enrichment_export, "KEGG", "KEGG_All.png", "KEGG Pathways (All DEGs)", "firebrick")


plot_enrichment(up_enrichment_export, "GO:BP", "GO_BP_Upregulated.png", "GO Biological Process (Upregulated)", "steelblue")
plot_enrichment(up_enrichment_export, "GO:MF", "GO_MF_Upregulated.png", "GO Molecular Function (Upregulated)", "forestgreen")
plot_enrichment(up_enrichment_export, "GO:CC", "GO_CC_Upregulated.png", "GO Cellular Component (Upregulated)", "darkorange")
plot_enrichment(up_enrichment_export, "KEGG", "KEGG_Upregulated.png", "KEGG Pathways (Upregulated)", "firebrick")



plot_enrichment(down_enrichment_export, "GO:BP", "GO_BP_Downregulated.png", "GO Biological Process (Downregulated)", "steelblue")
plot_enrichment(down_enrichment_export, "GO:MF", "GO_MF_Downregulated.png", "GO Molecular Function (Downregulated)", "forestgreen")
plot_enrichment(down_enrichment_export, "GO:CC", "GO_CC_Downregulated.png", "GO Cellular Component (Downregulated)", "darkorange")
plot_enrichment(down_enrichment_export, "KEGG", "KEGG_Downregulated.png", "KEGG Pathways (Downregulated)", "firebrick")









############################################################
# PART 7: STRING INPUT FILES
############################################################

############################################################
# KEGG Genes for STRING
############################################################

kegg <- subset(
  gost_res$result,
  source == "KEGG"
)

kegg_genes <- unique(
  unlist(kegg$intersection)
)

write.table(
  kegg_genes,
  "STRING_KEGG_Genes.txt",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

############################################################
# TOP 300 UPREGULATED GENES
############################################################

top300_up <- upregulated_annotated[
  order(
    upregulated_annotated$padj,
    -abs(upregulated_annotated$log2FoldChange)
  ),
]

top300_up <- head(top300_up, 300)

############################################################
# Use Gene Symbol if available, otherwise Gene ID
############################################################

top300_up$STRING_ID <- ifelse(
  is.na(top300_up$GeneSymbol) |
    top300_up$GeneSymbol == "",
  top300_up$GeneID,
  top300_up$GeneSymbol
)

############################################################
# Remove Duplicate STRING IDs
############################################################

top300_up <- top300_up[
  !duplicated(top300_up$STRING_ID),
]

############################################################
# Export STRING File
############################################################

write.table(
  top300_up$STRING_ID,
  "STRING_Top300_Upregulated.txt",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

############################################################
# TOP 300 DOWNREGULATED GENES
############################################################

top300_down <- downregulated_annotated[
  order(
    downregulated_annotated$padj,
    abs(downregulated_annotated$log2FoldChange)
  ),
]

top300_down <- head(top300_down, 300)

############################################################
# Use Gene Symbol if available, otherwise Gene ID
############################################################

top300_down$STRING_ID <- ifelse(
  is.na(top300_down$GeneSymbol) |
    top300_down$GeneSymbol == "",
  top300_down$GeneID,
  top300_down$GeneSymbol
)

############################################################
# Remove Duplicate STRING IDs
############################################################

top300_down <- top300_down[
  !duplicated(top300_down$STRING_ID),
]

############################################################
# Export STRING File
############################################################

write.table(
  top300_down$STRING_ID,
  "STRING_Top300_Downregulated.txt",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

############################################################
# Preview STRING Inputs
############################################################

head(
  top300_up[
    ,
    c("GeneID", "GeneSymbol", "STRING_ID")
  ],
  20
)

head(
  top300_down[
    ,
    c("GeneID", "GeneSymbol", "STRING_ID")
  ],
  20
)

############################################################
# Number of Genes with Official Symbols
############################################################

cat(
  "Upregulated genes with symbols:",
  sum(top300_up$GeneSymbol != top300_up$GeneID),
  "\n"
)

cat(
  "Downregulated genes with symbols:",
  sum(top300_down$GeneSymbol != top300_down$GeneID),
  "\n"
)


