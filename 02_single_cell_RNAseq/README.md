# Single-cell RNA-seq

Code for processing and downstream analysis of the adult oral-region single-cell RNA-seq dataset.

Raw and processed data are available through GEO under accession **GSE316974**, including the final Seurat object `nv_oral_seurat.RDS`.

Reference genome: *Nematostella vectensis* jaNemVect1.1 (GCF_932526225.1).

`figure_inputs/` contains curated gene tables used directly by the figure scripts. For Fig. 4, neurotransmitter receptor and sensory gene lists are provided in `Fig4_neurotransmitter_receptor_gene_lists/` and `Fig4_sensory_gene_lists/`; results from these analyses were manually combined into the corresponding files in `figure_inputs/`.
