# cyemapplot: Advanced Network Visualization for Enrichment Analysis

## Overview

`cyemapplot` is an R function that extends the capabilities of the [clusterProfiler](https://bioconductor.org/packages/clusterProfiler/) package by exporting functional enrichment analysis results to [Cytoscape](https://cytoscape.org/) for advanced network visualization and customization. Built upon the `emapplot` function from the [enrichplot](https://bioconductor.org/packages/enrichplot/) package, `cyemapplot` enables researchers to create publication-quality network visualizations that reveal relationships between enriched terms through shared gene content.

## Key Features

- **Three Visualization Styles:**
  - **Basic**: Standard network showing enriched terms and similarity relationships
  - **Pie**: Nodes display proportion of input genes within each enriched term
  - **DEG**: Enhanced visualization with gene expression directionality and regulation patterns

- **Cytoscape Integration:** Leverages Cytoscape's powerful layout algorithms, styling options, and network analysis tools

- **Flexible Network Filtering:** Filter networks by cluster size and visualize top components separately

- **Compatible with ORA and GSEA:** Works with both over-representation analysis and gene set enrichment analysis results

- **Multiple Databases Supported:** Disease Ontology, WikiPathways, Gene Ontology, KEGG, Reactome, and custom gene sets

### Prerequisites

1. **R** (version 4.5.1 or higher): [Download R](https://www.r-project.org/)
2. **RStudio** (recommended): [Download RStudio](https://posit.co/download/rstudio-desktop/)
3. **Cytoscape** (version 3.10.0 or higher): [Download Cytoscape](https://cytoscape.org/)

### Install Required R Packages

This repository uses `renv` for package management. After cloning the repository:

```r
# Install renv if not already installed
install.packages("renv")

# Restore the R environment with all required packages
renv::restore()
```

**Required packages:**
- clusterProfiler
- enrichplot
- DOSE
- org.Hs.eg.db
- RCy3
- igraph
- tidyverse

### Clone the Repository

```bash
git clone https://github.com/macsbio/cyemapplot-visualization.git
cd cyemapplot-visualization
```

Or download the ZIP file:
1. Click the green "Code" button
2. Select "Download ZIP"
3. Extract to your preferred location

### cyemapplot function parameters overview
## Function Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `ea_sim` | enrichResult or gseaResult object with termsim matrix | Required |
| `analysis_name` | Name for Cytoscape network collection | "Enrichment" |
| `show_category` | Number of top enriched terms to display | 30 |
| `min_edge` | Minimum similarity score for edge creation | 0.2 |
| `degs_data` | Data frame with differential expression data (required for "deg" style) | NULL |
| `visualization` | Visualization style: "basic", "pie", or "deg" | "basic" |
| `ig_layout` | igraph layout function | layout_with_kk |
| `layout_scale` | Scaling factor for node positions | 500 |
| `min_cluster_size` | Minimum cluster size to retain | 1 |
| `plot_clusters` | Create separate subnetworks for top clusters | FALSE |
| `top_clusters` | Number of top clusters to plot separately | 5 |






