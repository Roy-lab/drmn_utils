#!/bin/bash
# Runs a bunch of examples.
set -u 

finder=./findTransitionGenesetsDRMN

#### test for when we have fewer celltypes than in order file ####
drmnRes=example_input_missingcell
outNew=example_output_missingcell

${finder} ${drmnRes} example_data/celltype_order.txt example_data/ogids_exp.txt ips 0.05 ${outNew} 10

#### Test on single result directory from DRMN. ####
drmnRes=example_input_drmn
outNew=example_output_drmn

${finder} ${drmnRes} example_data/celltype_order.txt example_data/ogids_exp.txt ips 0.05 ${outNew} 10

#### Test on separate results dirs from RMN. ####
inSep="example_input_rmn/result_#CELL#"
outSep=example_output_rmn

${finder} ${inSep} example_data/celltype_order.txt example_data/ogids_exp.txt ips 0.05 ${outSep} 10

#### Run test of find_and_view
./find_and_view_transitionGenesets.sh example_input_drmn example_output_drmn example_data/celltype_order.txt example_data/ogids_exp.txt ips 3 0 14 /mnt/dv/wid/projects2/Roy-common/programs/scripts/figscripts/Heatmap.awk

./find_and_view_transitionGenesets.sh example_input_rmn/result_#CELL# example_output_rmn example_data/celltype_order.txt example_data/ogids_exp.txt ips 3 0 14 /mnt/dv/wid/projects2/Roy-common/programs/scripts/figscripts/Heatmap.awk
 
