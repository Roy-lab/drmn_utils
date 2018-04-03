# findTransitionGenesetsDRMN

Based on thefindTransitionGenesets utility, in https://github.com/Roy-lab/cmint

Debbie Chasman, April 2018


## Contents
This directory contains the following tools:

1. (C++ executable) ./findTransitionGenesetsDRMN is used to identify genes with non-constant trajectories, genesets with similar trajectories, and generate input files for making heatmap figures.
2. (Bash script) ./makeTransitionHeatmaps_DRMN.sh is a script that operates on the output of findTransitionGenesetsDRMN to make .svg figures for each set of genes and a bird's eye figure with the average pattern figure for all gene sets. 
3. (Bash script) ./find_and_view_transitionGenesets.sh is a wrapper around both 1. and 2. to make it easier to run in one go.

To plot the heatmaps, you will need the Heatmap.awk tool, written by Pouya Kheradpour. It is available from his website:
http://compbio.mit.edu/pouyak/software/Heatmap.awk

## findTransitionGenesetsDRMN
## --------------------------
This program will take the cell type-specific cluster assignments and cluster profiles and produce multiple files, one for each set of genes with a specific pattern. The size of the gene set is user defined to put a lower limit on the number of genes. 

In addition to per-geneset files defining heatmaps for heatmap.awk, it will produce three additional files:
-- all_clusterassign.txt, a matrix of cluster assignments for the genes with non-constant trajectories
-- all_assign.txt, a heatmap.awk-compatible file which has all the genes grouped by the specific patterns they exhibit (including gene sets including 1 gene). -- ordered_clusterset_means.txt, a heatmap.awk-compatible file that will have the mean cluster assignment profile and the mark profile for all gene sets
with more than the "minimum" number of genes. 

The results can be provided to makeTransitionHeatmaps_DRMN.sh for generating SVG images of the gene sets.

### USAGE
./findTransitionGenesetsDRMN drmn_result_dir_pattern celltypeorder OGID_file srccelltype threshold outputdir mingenesetsize

* drmn_result_dir_pattern: This is the location of the DRMN output directory, OR, a pattern that can be used to locate cell type-specific output directories, in the case where DRMN was run separately on each cell type (RMN). To specify a pattern, use the string #CELL# in place of the celltype name. The celltype name must match an entry in the celltypeorder file.

* celltypeorder: This is the same input as for DRMN and specifies the order of the cell types in the OGID file.

* OGID_file: This file is the same input as for DRMN and has the names of regions, one line per region.

* srccelltype: This is the root of the tree or the species with respect to which the cluster assignment files are created. 

* threshold: This is used in the hierarchical clustering to define clusters of regions with a specific pattern.

* outputdir: This directory stores the result of running findTransitionGenesetsDRMN

* mingenesetsize: This is the minimum number of genes a pattern should have to display its pattern in a separate clusterset file. The average pattern of this set is also outputted in the ordered_clusterset_means.txt file.

An example usage for a single DRMN output directory is:
./findTransitionGenesetsDRMN example_input_drmn example_data/celltype_order.txt example_data/ogids_exp.txt ips 0.05 example_out_drmn 10

An example usage for separate RMN output directories is:
./findTransitionGenesetsDRMN example_input_rmn/result_#CELL# example_data/celltype_order.txt example_data/ogids_exp.txt ips 0.05 example_out_rmn 10

Note 1:  when generating the average cluster profile for a cell type, there might be ties (we display the majority of the cluster assignment, but it is possible that for a cell type for those regions, there is no clear majority.). Then you will see a message: "Oops found a tie." This is a harmless message. The tie is resolved by using the max of the cluster ids.

## makeTransitionHeatmaps_DRMN.sh
## ------------------------------
This script will use the files in the output directory of findTransitionGenesetsDRMN to make heatmaps. 

It can also make a heatmap from all_assign.txt, although this line is commented out as the file can be very large. Please uncomment if you wish to make this figure.

### USAGE
./makeTransitionHeatmaps_DRMN.sh outputdir numberofclusters minExpr maxExpr heatmap_script_path
* outputdir: this is the output location directory of findTransistionGenesetsDRMN. All figures will be stored in a subdirectory called figs in outputdir.
* nuberofclusters: this is the number of clusters (expression states) that were specified when running DRMN.
* minExpr : A lower bound on the expression values. Can be positive or negative. Will be given a blue colormap.
* maxExpr : An upper bound on the expression values. Will be given a red colormap.
* heatmap_script_path: This is the location of the heatmap.awk script. We provide a version as part of the DRMN pacakge. This Heatmap.awk script is written by Pouya Kheradpour and was downloaded from http://compbio.mit.edu/pouyak/software/Heatmap.awk

Example usage:
./makeTransitionHeatmaps_DRMN.sh example_out_drmn 3 0 14 Heatmap.awk


## find_and_view_transitionGenesets.sh
## ------------------------------------
This is a wrapper script around findTransitionGenesetsDRMN and makeTransitionHeatmaps_DRMN.sh.
By default, it reports genesets of size 3+, using a hierarchical clustering threshold of 0.05. You may edit these in the script.
The output directory will be appended with "_min3_0.05" (or whatever the settings are for the min geneset size and threshold).

### USAGE
./find_and_view_transitionGenesets.sh drmn_results_dir output_dir_prefix orderfile ogidsfile srcCell k minExpressionValue maxExpressionValue path_to_Heatmap.awk

Example usages:
./find_and_view_transitionGenesets.sh example_input_drmn example_output_drmn example_data/celltype_order.txt example_data/ogids_exp.txt ips 3 0 14 Heatmap.awk
./find_and_view_transitionGenesets.sh example_input_rmn/results_#CELL# example_output_rmn example_data/celltype_order.txt example_data/ogids_exp.txt ips 3 0 14 Heatmap.awk


Acknowledgements
-----------------
We thank Pouya Kheradpour for writing the Heatmap.awk script.


