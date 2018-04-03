#!/bin/bash
# Test to make sure we can now read from cell type-specific directories, 
# and that we get the same results as the original version did on DRMN.
set -u 

oldFinder=/mnt/dv/wid/projects2/Roy-common/programs/programs/cmint_all/utils/find_transitioning_genesets_local/findTransitionGenesets
newFinder=/mnt/dv/wid/projects2/Roy-common/programs/programs/cmint_all/utils/find_transitioning_genesets_DRMN/findTransitionGenesetsDRMN

#### Test on single result directory (original mode) ####
#### This should match the results from the oldFinder, with the exception that we add
#### a new output file: all_clusterassign.txt.
drmnRes=example_input_drmn

outOld=example_output_drmn_orig
outNew=example_output_drmn_update

mkdir -p $outOld
${oldFinder} ${drmnRes} example_data/celltype_order.txt example_data/ogids_exp.txt ips 0.05 ${outOld} 10

${newFinder} ${drmnRes} example_data/celltype_order.txt example_data/ogids_exp.txt ips 0.05 ${outNew} 10


diff -r --brief $outOld $outNew
if [ $? -ne 0 ]; then
	echo "Single results directory test failed."
fi

#### Test on separate results dirs ####
inCombined=example_input_rmn_singledir
inSep="example_input_rmn/result_#CELL#"
echo $inSep

outCombined=example_output_rmn_singledir
outSep=example_output_rmn

${newFinder} ${inCombined} example_data/celltype_order.txt example_data/ogids_exp.txt ips 0.05 ${outCombined} 10
${newFinder} ${inSep} example_data/celltype_order.txt example_data/ogids_exp.txt ips 0.05 ${outSep} 10

diff -r --brief $outCombined $outSep
if [ $? -ne 0 ]; then
	echo "Separate results directory test failed: "
	exit 1
fi
