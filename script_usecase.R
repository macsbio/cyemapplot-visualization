#install the renv ppackage to restore the R enviroment to install the required R packages
install.packages("renv")

#load the renv R package
library(renv)

#install the R packages
renv::restore()

#Load required libraries
library(clusterProfiler)
library(enrichplot)
library(RCy3)
library(org.Hs.eg.db)
library(tidyverse)
library(dplyr)
library(DOSE)
library(igraph)

# Set working directory  
# Note: Replace `"path/to/your/project"` with the actual path to the directory 
# where you downloaded and extracted the GitHub repository.  
setwd("path/to/your/project")  
#setting the working directory (auto)
#setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) if scripts run from the directory thsi can be used

# Source the cyemapplot function  
source("cyemapplot.R") 

##create output directory
outdir <- "Output-plots"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

#Note: First download and open Cytoscape application
# Test connection to Cytoscape 
cytoscapePing() 
#Expected output: "You are connected to Cytoscape!" 

# Check Cytoscape version 
cytoscapeVersionInfo() 

#import usecase deg data 
ibm.data<- read.delim("res_sIBM_AMP_main_updated.txt", sep = " ")
ibm.data$ID<- rownames(ibm.data)
ibm.data <- ibm.data %>%
  separate(ID, into = c("ensembl_id", "symbol"), sep = ";", remove = FALSE)
ibm.data$ensembl_id <- sub("\\..*$", "", ibm.data$ensembl_id)
ibm.data<- ibm.data[complete.cases(ibm.data),]

ibm.data<- subset(ibm.data, select = c("ensembl_id", "log2FoldChange", "pvalue", "padj"))
rownames(ibm.data)<- ibm.data$ensembl_id
colnames(ibm.data)<- c("ID", "log2FC", "pvalue", "padj")

# Convert ENSEMBL to ENTREZID for enrichment analysis  
mapped_ids <- bitr(ibm.data$ID, fromType = "ENSEMBL", toType = "ENTREZID",
                   OrgDb = org.Hs.eg.db)

# Merge converted IDs back to original data 
data_mapped <- merge(ibm.data, mapped_ids, by.x = "ID", by.y = "ENSEMBL")
data_mapped$ID<- NULL
data_mapped$ID <- data_mapped$ENTREZID 
data_mapped$ENTREZID <- NULL
data_mapped<-data_mapped[!duplicated(data_mapped$ID),]

#Define significance thresholds 
sig.deg <- subset(data_mapped, pvalue < 0.05 & abs(log2FC) > 1) 

# Extract gene IDs  
input_genes <- sig.deg$ID 

# Define background universe  
background_genes <- data_mapped$ID 


# Perform Disease Ontology enrichment 
ora_do <- enrichDO(gene = input_genes, 
                   ont = "HDO", 
                   organism = "hsa", 
                   pvalueCutoff = 0.05, 
                   pAdjustMethod = "BH",
                   universe = background_genes, 
                   maxGSSize = 500, 
                   qvalueCutoff = 0.05)


# Perform WikiPathways enrichment 
ora_wp <- enrichWP(gene = input_genes, 
                   organism = "Homo sapiens",
                   maxGSSize = 500, 
                   pAdjustMethod = "BH", 
                   universe = background_genes) 


# Create ranking metric, e.g. log2FC × -log10(p-value)
data_mapped$ranking <- data_mapped$log2FC * -log10(data_mapped$pvalue) 

# Prepare named numeric vector for GSEA 
geneList <- data_mapped$ranking  
names(geneList) <- data_mapped$ID 

# Sort in decreasing order 
geneList <- sort(geneList, decreasing = TRUE) 

#Perform GSEA on GO Biological Process
gsea_go_bp <- gseGO(geneList, 
                    OrgDb = "org.Hs.eg.db", 
                    keyType = "ENTREZID", 
                    ont = "BP", 
                    pvalueCutoff = 0.01, 
                    pAdjustMethod = "fdr")

# Dot plot for ORA results 
dotplot(ora_do, showCategory = 30) + ggtitle("Disease Ontology Enrichment") 


# Dot plot for GSEA results 
dotplot(gsea_go_bp, showCategory = 30) + ggtitle("GO Biological Process Enrichment") 

#pairwise term calculculation
do.ora.sim <- pairwise_termsim(ora_do, showCategory = 328, method = "JC")

# Bar plot for ORA results 
barplot(ora_do, showCategory = 30) + ggtitle("Disease Ontology Enrichment")

# Bar plot for GSEA results
barplot(gsea_go_bp, showCategory = 30) + ggtitle("GO Biological Process Enrichment")

# Calculate pairwise term similarity for ORA results 
ora_sim <- pairwise_termsim(ora_do, method = "JC", showCategory = nrow(ora_do)) 

# Calculate pairwise term similarity for GSEA results 
gsea_sim <- pairwise_termsim(gsea_go_bp, method = "JC", showCategory = nrow(gsea_go_bp)) 

# Tree plot visualization 
treeplot(ora_sim, 
         showCategory = 50, 
         cluster_method = 	"ward.D", 
         cladelab_offset = rel(6.5))


# Enrichment map visualization 
p <- emapplot(ora_sim, 
              showCategory = 100, 
              min_edge = 0.4,
              node_label = "all", 
              nWords = 4) 

# Adjust plot margins for better label visibility
p + coord_equal(clip = "off") +
  scale_x_continuous(expand = expansion(mult = 0.15)) + 
  scale_y_continuous(expand = expansion(mult = 0.15)) + 
  theme(plot.margin = margin(10, 60, 10, 10))


###cyemapplot visualisation
# Calculate pairwise term similarity for ORA results 
do.ora.sim <- pairwise_termsim(ora_do, method = "JC", showCategory = nrow(ora_do)) 

# Create basic Cytoscape visualization 
cyemapplot(do.ora.sim, 
           show_category = nrow(ora_do), # Number of terms to display 
           min_edge = 0.4,				# Minimum similarity threshold
           visualization = "basic",		# Basic style 
           ig_layout = igraph::layout_with_kk,	# Kamada-Kawai layout
           layout_scale = 800,			# Layout scaling factor 
           min_cluster_size = 8,			# Min cluster size filter 
           analysis_name = "IBM-DO-ORA") 

# Create pie chart Cytoscape visualization 
cyemapplot(wp.ora.sim,  
           show_category = nrow(ora_wp),		# Number of enriched pathways 
           min_edge = 0.2,				# Lower threshold for pathway overlap 
           visualization = "pie",			# Pie chart style 
           ig_layout = igraph::layout_with_kk, 
           layout_scale = 800, 
           min_cluster_size = 1,			# Keep all components (no size filter) 
           analysis_name = "IBM WP-ORA") 

# Calculate pairwise term similarity for GSEA results 
gsea.go.bp.sim <- pairwise_termsim(gsea_go_bp, method = "JC", showCategory = nrow(gsea_go_bp)) 

# Create DEG visualization in Cytoscape 
cyemapplot(gsea.go.bp.sim, 
           show_category = nrow(gsea_go_bp),	# All enriched GO terms 
           min_edge = 0.4,                	# Moderate similarity threshold 
           visualization = "deg",          	# DEG style with expression info
           degs_data = data_mapped,          	# Differential expression data 
           ig_layout = igraph::layout_with_kk, 
           layout_scale = 800, 
           min_cluster_size = 3,           	# Remove small components
           plot_clusters = TRUE,		# Create subnetworks for top clusters
           top_clusters = 5,       		# Top 10 clusters  
           analysis_name = "IBM GO-GSEA") 

