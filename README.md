# Bioinformatics Projects 


This repository contains a variety of bioinformatics analysis workflows developed by MSc. Elif Duymaz Yılmaz. Each subdirectory houses a specialized analysis pipeline or code example.

## 📂 Project Structure

```
Bioinformatics/
├── scrnaseq_analysis-dockerImage/   # Docker + R-based single-cell RNA-Seq analysis
├── seurat.v5-todo/                  # Example script for clustering with Seurat v5
├── alternative-splicing-scrna-seq/  # Alternative splicing analysis pipeline (Python & R)
├── max-delbruck-center-march25/     # Spatial and multi-omics modules
└── README.md                        # This file
```

### 1. scrnaseq_analysis-dockerImage/
- `Dockerfile` → Base image definition for the analysis environment (R, Seurat, FastQC, Kallisto, etc.)
- `main.sh` → Pipeline control script (download, QC, alignment, quantification)
- `download.sh`, `quality_check.sh`, `kallisto.sh` → Step-by-step Bash scripts
- `Rconf.R`, `seurat.R` → Seurat-based data processing and visualization
- `DemoDockerImage` → Example container startup command

### 2. seurat.v5-todo/
- `GeneXClusterMatrix.R` → Example R script to generate gene expression–cluster matrices with Seurat v5

### 3. alternative-splicing-scrna-seq/
- `datadownload.sh` → Bash script for downloading raw data
- `scats.py`, `scats/` → Python implementation of the SCATS algorithm for single-cell splicing analysis
- `figures.R` → R script for plotting analysis results
- `README.md` → Project description and usage instructions

### 4. max-delbruck-center-march25/
- `spatial-omics-analysis/` → Modules for spatial omics data processing (under development)
- `multi-omics-data-integration-using-deeplearning/` → Deep learning–based multi-omics integration (CBioPortal and dimensionality reduction examples)

## 🚀 Usage

1. Navigate to the desired subdirectory:
   ```bash
   cd scrnaseq_analysis-dockerImage
   ```
2. Build and run the Docker environment:
   ```bash
   docker build -t scrna_pipeline .
   docker run --rm -v $(pwd)/data:/data scrna_pipeline bash main.sh
   ```
3. Alternatively, execute individual Bash or Nextflow scripts step by step.

## 🛠️ Requirements

- Docker or Singularity
- R (>=4.0) with required packages (`Seurat`, `scater`, `SCATS`)
- Python 3.7+ (for SCATS)
- Bash shell

---

*All rights reserved.*
