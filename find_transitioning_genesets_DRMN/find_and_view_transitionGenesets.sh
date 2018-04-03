#!/bin/bash
# Runs findTransitionGenesetsDRMN to identify 
set -u
if [[ $# != 9 ]]; then
	echo "USAGE: find_and_view_transitionGenesets.sh drmn_results_dir output_dir_prefix orderfile ogidsfile srcCell k minExpressionValue maxExpressionValue path_to_Heatmap.awk"
	echo "drmn_results_dir may be a single directory, or a string pattern for cell type-specific directories, containing #CELL# where each cell type specific result is located."
	echo "We will be appending the geneset size and threshold to the output dir prefix."
	exit 2
fi

results=$1	# DRMN results directory
out=$2
order=$3	# celltype order (for columns)
ogids=$4	# ogids file
srccell=$5	# source celltype for names
k=$6
minExp=$7   # minimum expression for color scale
maxExp=$8 	# maximum expression for color scale
HMAWK=$9	# path to heatmap.awk

threshold=0.05
mingeneset=3

echo $results
#if [[ ! -d ${results} ]]; then
#	echo "No directory: $results"
#	exit
#fi

findTransitionGenesets=./findTransitionGenesetsDRMN
makeHeatmap=./makeTransitionHeatmaps_DRMN.sh

# Run findTransitionGenesets
out=${out}_min${mingeneset}_${threshold}
mkdir ${out}
echo $out
${findTransitionGenesets} ${results} ${order} ${ogids} ${srccell} $threshold ${out} $mingeneset

# If results exist, plot.
if [[ -e ${out}/all_genesets.txt ]]; then
	ngenes=$(cat ${out}/all_clusterassign.txt | wc -l)
	nclust=$(cat ${out}/all_genesets.txt | wc -l)
	echo "Found $ngenes total genes that change state."
	echo "Found $nclust genesets with at least $mingeneset genes."
	${makeHeatmap} ${out} $k $minExp $maxExp $HMAWK
fi

