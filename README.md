# ML-DS: Myeloid Leukemia of Down Syndrome Analysis

This repository contains computational analysis code for single-cell RNA sequencing (scRNA-seq) data processing and analysis of Myeloid Leukemia of Down Syndrome (ML-DS). The code supports data processing, analyses, and figure generation for the associated publication.

## Overview

This project focuses on understanding the molecular mechanisms underlying ML-DS progression by analyzing:
- **Transient Abnormal Myelopoiesis (TAM)** - a pre-leukemic condition in Down syndrome newborns
- **ML-DS** - the acute leukemia that can develop from TAM
- **Fetal liver development** - normal hematopoietic development in diploid and trisomy 21 (Down syndrome) contexts
- **Megakaryocyte lineage trajectories** - developmental pathways from hematopoietic stem cells to mature megakaryocytes

## Key Analyses

### 1. Data Processing (`main_analyses/01_fetalLivers_processing/`, `main_analyses/02_MLDS_processing/`, `main_analyses/03_otherLeukaemia_processing/`)
- Quality control and preprocessing of scRNA-seq data
- Cell type annotation using reference-based label transfer
- Sample metadata integration and batch effect correction
- Processing of fetal liver references, ML-DS samples, and other leukemia datasets

### 2. Fetal Aneuploidy Analysis (`main_analyses/01.2_fetal_aneuploidy_analyses/`)
- Analysis of trisomy 21 effects on fetal hematopoietic development
- Comparison of diploid vs. trisomy 21 developmental trajectories
- Early molecular changes in Down syndrome hematopoiesis

### 3. Transcriptional Module Derivation (`main_analyses/04_derive_transcriptional_modules/`)
- **GATA1s module analysis** (`x4_GATA1s_module.R`) - characterization of truncated GATA1 regulatory network
- **ML-DS imprint in T21 fetal liver** (`x3_MLDS_imprint_in_fLiver_T21_v2.R`) - early molecular changes predisposing to leukemia
- **TAM vs ML-DS comparison** (`x5_good_vs_bad_TAM.R`, `x6_TAM_vs_MLDS.R`) - progression markers

### 4. Module Specificity Testing (`main_analyses/05_transcriptional_modules_specificity/`)
- Cross-dataset validation of transcriptional modules
- Module scoring across different cell types and conditions
- Bulk RNA-seq validation of single-cell derived signatures

### 5. ML-DS Relapse/Refactory Analysis (`main_analyses/06_MLDS_refractory_relapse/`)
- Analysis of treatment response and relapse patterns
- Identification of markers associated with clinical outcomes

### 6. Trajectory Analysis (`eda_scripts/`)
- **Palantir-based pseudotime analysis** (`2.1_fLiver_trajectory_2401.R`) of fetal liver hematopoietic development
- Construction of megakaryocyte/erythroid/mast cell developmental trajectories (`MLDS_trajectory_projection/`)
- Projection of TAM/ML-DS cells onto normal developmental trajectories

## Repository Structure

```
├── main_analyses/                              # Primary analysis pipeline
│   ├── 01_fetalLivers_processing/             # Fetal liver reference processing
│   ├── 01.2_fetal_aneuploidy_analyses/        # Trisomy 21 effects analysis
│   ├── 02_MLDS_processing/                    # ML-DS dataset processing
│   ├── 03_otherLeukaemia_processing/          # Other leukemia datasets
│   ├── 04_derive_transcriptional_modules/     # Module identification
│   ├── 05_transcriptional_modules_specificity/ # Module validation
│   ├── 06_MLDS_refractory_relapse/            # Clinical outcomes
│   ├── figures_generation/                    # Publication figures
│   └── utils/                                 # Shared utility functions
├── eda_scripts/                               # Exploratory analyses
│   ├── MLDS_trajectory_projection/            # Trajectory analysis
│   └── *.R                                    # Various exploratory scripts
└── LICENSE                                    # License file
```

## Key Dependencies

### R Packages
- **Seurat** - scRNA-seq analysis framework
- **DESeq2** / **edgeR** - differential expression analysis
- **UCell** - gene module scoring
- **ComplexHeatmap** - advanced visualization
- **tidyverse** - data manipulation and visualization

### Python Packages
- **Palantir** - trajectory inference and pseudotime analysis
- **scanpy** - single-cell analysis in Python
- **pandas** / **numpy** - data manipulation

### External Tools
- **CellTypist** - automated cell type annotation
- **SoupX** - ambient RNA removal
- **Harmony** - batch effect correction

## Data Types

The analysis incorporates multiple datasets:
- **ML-DS patient samples** - diagnostic and follow-up timepoints
- **TAM patient samples** - transient abnormal myelopoiesis
- **Fetal liver references** - normal hematopoietic development (diploid and T21)
- **Fetal adrenal** - additional developmental reference
- **Infant ALL** - pediatric leukemia comparison
- **Published atlases** - external validation datasets

## Usage

### Prerequisites
1. Install required R and Python packages
2. Configure file paths in scripts to match your data organization
3. Ensure access to reference genomes and annotation files

### Basic Workflow
1. **Data Processing**: Run scripts in `main_analyses/01_fetalLivers_processing/`, `main_analyses/02_MLDS_processing/`, and `main_analyses/03_otherLeukaemia_processing/`
2. **Aneuploidy Analysis**: Execute scripts in `main_analyses/01.2_fetal_aneuploidy_analyses/`
3. **Transcriptional Modules**: Run module derivation scripts in `main_analyses/04_derive_transcriptional_modules/`
4. **Module Validation**: Execute validation scripts in `main_analyses/05_transcriptional_modules_specificity/`
5. **Clinical Analysis**: Run outcome analysis scripts in `main_analyses/06_MLDS_refractory_relapse/`
6. **Trajectory Analysis**: Run exploratory trajectory analysis in `eda_scripts/MLDS_trajectory_projection/`
7. **Visualization**: Generate publication figures using `main_analyses/figures_generation/`

## Citation

If you use this code, please cite the associated publication:

*[Publication details to be added upon acceptance]*

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contact

For questions about the analysis or code implementation, please contact the authors through the associated publication or repository issues.

