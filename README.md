# Nematostella oral nerve ring analysis

Code and analysis files associated with:

**Morphological, molecular, and behavioral characterization of a ganglion-like oral nerve ring in the sea anemone Nematostella vectensis**

Ruohan Zhong, Chris W. Seidel, Anna M.L. Klompen, and Matthew C. Gibson

This repository contains custom analysis and figure-generation code used in the study, together with small input files required for the analyses.

## Repository structure

- `01_cell_body_quantification/`  
  Quantification, plotting, and statistical analysis of oral nerve ring cell-body size and density.

- `02_single_cell_RNAseq/`  
  Single-cell RNA-seq processing and downstream analyses used for Figs. 3, 4, 6a, S3, S5, and S6.

- `03_bulk_RNAseq/`  
  Bulk RNA-seq TPM processing and analysis used for Fig. 5c.

- `04_feeding_behavior/`  
  Plotting and statistical analysis of feeding behavior used for Figs. 6e-f and S8d-f.

- `05_phylogeny/`  
  Files and documentation associated with phylogenetic analyses.

## Data availability

Sequencing data are available through GEO:

- Single-cell RNA-seq: **GSE316974**
- Bulk RNA-seq: **GSE316975**

The processed single-cell Seurat object, `nv_oral_seurat.RDS`, is provided with the single-cell RNA-seq dataset and is not included in this repository.

Additional source data and imaging data associated with the study are available through the Stowers Original Data Repository (ODR).

## Reference genome

RNA-seq analyses used the *Nematostella vectensis* jaNemVect1.1 reference genome assembly:

**NCBI accession: GCF_932526225.1**

## Running the analyses

Scripts are organized by analysis type and figure. Small gene lists and figure input files are included in the corresponding analysis directories.

Large sequencing datasets and processed objects should be downloaded from the associated data repositories as needed.

## Citation

If you use this code, please cite the associated publication:

Zhong R, Seidel CW, Klompen AML, Gibson MC.  
*Morphological, molecular, and behavioral characterization of a ganglion-like oral nerve ring in the sea anemone Nematostella vectensis.*  
Nature Communications. [publication DOI to be added]

Archived code release:

[Zenodo DOI to be added]
