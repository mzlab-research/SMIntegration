# SMIntegration : Spatial Multi-omics Integration Platform
[![Docker Pulls](https://img.shields.io/docker/pulls/yourdocker/smintegration)](https://hub.docker.com/r/yourdocker/smintegration)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

## Introduction
**SMIntegration** is an innovative open-source platform for integrated analysis of spatial transcriptomics and metabolomics data. It unifies spatial pattern recognition, differential comparison, network construction, and functional annotation in a single environment. Designed to address key challenges in spatial multi-omics correlation analysis, SMIntegration enables researchers to explore gene-metabolite co-regulation mechanisms through an intuitive web interface, revealing spatial heterogeneity in tissue development and disease progression.

## Key Features
- **Multi-modal Spatial Pattern Discovery**
  - Automated identification of spatial expression patterns for genes/metabolites
  - Moran's I correlation analysis between transcriptomic/metabolomic modules
- **Pixel-level Clustering & Cell Annotation**
  - 4 clustering algorithms: Seurat-LV, Seurat-LM, Seurat-SLM, UMAP-kmeans
  - SingleR-based automatic cell type annotation with reference datasets
  - Custom cell type annotation support
- **Flexible Differential Analysis**
  - ROI selection via: 
    - Interactive tissue imaging
    - Clustering results 
    - Cell annotation mapping
  - Differentially expressed genes (DEGs) and metabolites (DAMs) detection
- **Condition-Specific Networks**
  - Group-specific gene-metabolite correlation networks (|r| > 0.6, adj.p < 0.01)
  - Comparative network visualization across conditions
- **Interactive Visualization**
  - Gene-metabolite co-localization analysis
  - RGB overlay imaging for multi-feature visualization
  - Spatial distribution heatmaps for single molecules
- **Functional Enrichment**
  - KEGG pathway co-enrichment analysis
  - Fisher's exact test for DEG-DAM pathway associations

## Getting Started
### Web Platform (Recommended for small/medium datasets)
Access our cloud deployment directly:  
🔗 [https://metax.genomics.cn/app/SMIntegration](https://metax.genomics.cn/app/SMIntegration)

### Local Deployment (Docker)
For large datasets (>50GB memory requirement):
```bash
# Pull the latest Docker image
docker pull yourdocker/smintegration:latest

# Run with ShinyProxy
docker run -d -p 8080:8080 shinyproxy/smintegration-app
