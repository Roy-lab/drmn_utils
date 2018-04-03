# drmn_utils
Ancillary utilities for DRMN. Input processing and interpreting results.


## find_transitioning_genesets_DRMN/findTransitionGenesetsDRMN 
This tool processes the output of DRMN to identify genes that transition states between cell types, and to find clusters of such genes with similar (or identical) trajectories. It generates a geneset file, a matrix of the trajectories for the transitioning genes, and several files that define heatmaps to plot with Heatmap.awk.
