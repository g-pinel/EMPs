# EMPs

This repository contains the script used to run this Shiny App:

This web app is related to our publication *Discovery of New Markers for Haemogenic Endothelium and Haematopoietic Progenitors in the Mouse Yolk Sac* ([pubmed link](https://pubmed.ncbi.nlm.nih.gov/41562864/)).

It allows exploration of our refined clusters of haematoendothelial cell types in embryonic day (E) 8.5 scRNA-seq data. The dashboard shown upon querying a gene includes:
- Feature plot: log-normalised expression level on a UMAP representation.
- Violin plot: log-normalised expression levels per cluster.
- Percentage plot: counts and percentage of cells per cluster that express the queried gene at any level.
- Subcellular location: a table with information regarding the subcellular location of the protein encoded by the queried gene (if protein-coding).

In addition, the app integrated differential-expression analysis tools, the following are supported:
- 1 cluster versus the rest
- One or more clusters versus another or more clusters (customize your comparisons).
