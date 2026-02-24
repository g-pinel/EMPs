# EMPs

This repository contains the script used to run the following Shiny App: <https://g-pinel.shinyapps.io/Haematoendothelial_E8_5/>

This web app is related to our publication *Discovery of New Markers for Haemogenic Endothelium and Haematopoietic Progenitors in the Mouse Yolk Sac* ([pubmed link](https://pubmed.ncbi.nlm.nih.gov/41562864/)).

It allows exploration of single-cell RNAseq data from embryonic day 8.5 mouse samples using our refined clusters of haematoendothelial cell types. The dashboard shown upon querying a gene includes:
- Feature plot: log-normalised expression level on a UMAP representation.
- Violin plot: log-normalised expression levels per cluster.
- Counts/Percentage plot: total counts and percentage of cells per cluster that express the queried gene at any level.
- Subcellular location: a table with information regarding the subcellular location of the protein encoded by the queried gene (if protein-coding).

In addition, the app integrated differential-expression analysis tools, the following are supported:
- 1 cluster versus the rest
- One or more clusters versus another or more clusters (customize your comparisons).

> [!NOTE]
> To accomodate the memory limit of 1Gb at shinyapps.io, erythrocytes and allantois endothelial cells were downscaled to 695 cells each, matching the third most abundant cell type in the dataset (embryo endothelium).
>
> However, this only affects the feature plot, violin plot, and custom cluster differential expression analyses, which generate data dynamically. 
