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
docker pull crpi-pjcujkw4652f020i.cn-shenzhen.personal.cr.aliyuncs.com/mzlab/smintegration:v-1.0<img width="692" height="19" alt="image" src="https://github.com/user-attachments/assets/2a9a7ef3-c99f-42eb-9083-d710fe2ec132" />

# Run with ShinyProxy
docker run -d -p 8787:3838  crpi-pjcujkw4652f020i.cn-shenzhen.personal.cr.aliyuncs.com/mzlab/smintegration:v-1.0
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
 - Acquisition & Processing: AFAD-ESI (+) mode, 50μm resolution, spatially registered to identical pixel coordinates
 - Content: 14,530 pixels × 500 metabolites
Access:
Built-in dataset in SMIntegration (Tutorial Panel → Example Data)
Raw data: [NGDC OMIX Repository](https://ngdc.cncb.ac.cn/omix) ID: OMIX009541

## Community & Support

Developed by Haoke Deng (denghaoke\@genomics.cn)\
Last update: 2025-08-11
