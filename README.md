# SMIntegration : Spatial Multi-omics Integration Platform
## Introduction
**SMIntegration** is an innovative open-source platform for integrated analysis of spatial transcriptomics and metabolomics data. It integrates spatial pattern recognition, differential comparison, network construction, and functional annotation into a unified workflow. Designed to address key challenges in spatial multi-omics correlation analysis, SMIntegration enables researchers to explore gene-metabolite co-regulation mechanisms through an intuitive web interface, revealing spatial heterogeneity in tissue development and disease progression.

## Key Features
### 🔍 **Spatial Pattern Discovery**
  - Automated identification of spatial expression patterns for genes/metabolites
  - Moran's I correlation analysis between transcriptomic/metabolomic modules
### 🧩 **Pixel-level Clustering & Cell Annotation**
  - 4 clustering algorithms: Seurat-LV, Seurat-LM, Seurat-SLM, UMAP-kmeans
  - SingleR-based automatic cell type annotation with reference datasets
  - Custom cell type annotation support
### ⚖️ **Flexible Differential Analysis**
  - ROI selection via: 
    - Interactive tissue imaging
    - Clustering results 
    - Cell annotation mapping
  - Differentially expressed genes (DEGs) and metabolites (DAMs) detection
  - Group-specific gene-metabolite correlation networks
### 🧬 **Functional Enrichment**
  - KEGG pathway co-enrichment analysis
  - Fisher's exact test for DEG-DAM pathway associations
### 📊 **Interactive Visualization**
  - Gene-metabolite co-localization analysis
  - RGB overlay imaging for multi-feature visualization

## Getting Started

### 🌐 Online Access  
Access the live platform without installation:
🔗 [[[https://metax.genomics.cn/app/SMIntegration](https://metax.genomics.cn/app/smintegration)]

### 💻 Local Deployment (Docker)  
🔗 [[[Docker Installation Guide](https://docs.docker.com/get-started/get-docker/)]
```bash
# Pull the latest Docker image
docker pull mzlabresearch/smintegration:v-1.0

# Run with ShinyProxy
docker run -d -p 8787:3838  mzlabresearch/smintegration:v-1.0
```
Please use a web browser to access: http://localhost:8787

### 📥 Input Data Format
SMIntegration requires two **feature matrices** as input files:

**Key Requirements**:

1、Pre-aligned datasets
Spatial metabolomics + transcriptomics must share identical pixel coordinates，we recommend using [SpatialData](SpatialData.md) for registration.

2、Supported Formats:
SMIntegration requires two **feature matrices** in TXT or RDS format:

**TXT Format**:
  - Columns 1-4: Feature name (metabolite or gene), pixel coordinates (x, y), feature abundance (Intensity or MIDCount)
  - Each row represents one spatial pixel

**RDS Format (Seurat object)**:
  - Spatial location information stored in meta.data (x, y)
  - Feature matrix stored in data@assays$Spatial (Format: Features as rows, spatial pixels as columns, values represent feature abundance)
> 📖 Detailed formatting guide available in-app (Tutorial Panel → Data Preparation)

## Usage Workflow

1. **Upload Data**: Import feature matrix in **Overall Distribution Analysis** Panel
2. **Core Analysis**: 
 - Pattern Recognition: Identify spatial expression modules
 - Pixel Clustering: Group pixels using 4 algorithms
 - Cell Annotation: Automated (SingleR) or manual mapping
 - Differential Analysis: Compare regions (manual/automatic ROI selection)
 - Network Construction: Build condition-specific gene-metabolite networks
 - Functional Enrichment: Pathway mapping and co-enrichment analysis
 - Visualization：Dynamic exploration of spatial distributions and co-localization patterns

## Example Data

Test dataset: **Mouse brain coronal adjacent sections**
1. **Spatial Transcriptomics Data**: 
 - Acquisition & Processing: Stereo-seq, 0.05μm resolution, aggregated to 50μm resolution
 - Content: 14,530 pixels × 500 genes
2. **Spatial Metabolomics Data**: 
 - Acquisition & Processing: AFADESI (+) mode, 50μm resolution, spatially registered to identical pixel coordinates
 - Content: 14,530 pixels × 500 metabolites

Access:
Built-in dataset in SMIntegration (Tutorial Panel → Example Data)

Raw data: NGDC OMIX Repository ID: [OMIX011674](https://ngdc.cncb.ac.cn/omix/release/OMIX011674)

## Developer Guide

This section is intended for developers who wish to understand the software architecture or contribute to the SMIntegration codebase.

### 1. Software Architecture

SMIntegration is built as a modular Shiny application. The codebase is organized to separate the User Interface (UI), Server Logic, and Core Analytical Functions.

*   **`app.R`**: The main entry point. It initializes the application, loads required R libraries, and sources the UI and Server components.
*   **`ui/`**: Contains the visual layout definitions for each module.
    *   Files are named `U[X]_[Module].R` (e.g., `U1_Tutorial.R`, `U4_clustering.R`).
*   **`server/`**: Contains the backend reactive logic corresponding to each UI module.
    *   Files are named `S[X]_[Module].R` (e.g., `S1_Tutorial.R`, `S4_clustering.R`).
*   **`source/`**: Houses utility functions and core algorithms used across the application.
*   **`www/`**: Stores static assets such as images and CSS styles.

### 2. Validation Suite

To ensure reproducibility and reliability, SMIntegration includes a comprehensive dual-layer validation framework.

#### Analytical Pipeline Validation
*   **Entry Point**: `run_pipeline_validation.R`
*   **Description**: This script executes a series of unit tests located in the `validation_pipeline/` directory. It uses the `testthat` package to verify the correctness of core algorithms (preprocessing, clustering, differential analysis, etc.) using the included demo dataset.
*   **How to Run**:
    ```r
    Rscript run_pipeline_validation.R
    ```

#### Visualization Verification
*   **Entry Point**: `run_validation_figure.R`
*   **Description**: A regression test that automatically regenerates key manuscript figures from archived intermediate data (`.rds` files). This ensures that the visualization pipeline remains consistent.
*   **How to Run**:
    ```r
    Rscript run_validation_figure.R
    ```
    
## Community & Support

Developed by Haoke Deng (denghaoke\@genomics.cn)\
Last update: 2025-08-11
