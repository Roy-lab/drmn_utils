#!/bin/bash
# Example script to test transitioning genesets for enrichment
set -u

home=$1 # Should be output for from find_transition_genesets
cutoff=0.05

enrich=/mnt/dv/wid/projects2/Roy-common/programs/programs/enrichanalyzer_Nongraph_Qval/enrichAnalyzer

cisbp=/mnt/dv/wid/projects2/Roy-common/data/data_new/mouse/cisbp_motifs/mammals_mus_musculus_common_regnet.txt 

outpref=${home}/enrichment/cisbp_mus_musculus
mkdir -p ${home}/enrichment

$enrich ${home}/all_genesets.txt <(cut -f1 ${home}/all_clusterassign.txt) $cisbp $cutoff ${outpref} persg








