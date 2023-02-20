# Marker gene selection method comparison based on GO enrichment and Information Content

This repository includes some of the code used for the paper: **Comparing marker gene selection methods for single-cell RNA-seq
data in terms of GO enrichment**.

Files included:

- ``preprocessing.R``: Load and preprocessing of the ``SingleCellExperiment`` objects. Includes QC, normalization, feature selection, PCA, clustering. Packages needed: ``SingleCellExperiment``, ``scran``, ``scuttle``.
