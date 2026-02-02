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

