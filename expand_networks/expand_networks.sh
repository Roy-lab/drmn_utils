#!/bin/bash
# Given the regulator-module networks and the module assignments from DRMN,
# joins them to create a regulator-gene network per cell type
# This is a placeholder to tide us over until we have the main DRMN code spit this out.
# It will output the networks in place as "cell_expanded_net.txt"
# USAGE: expand_networks.sh drmn_results_dir celltype_order.txt
set -u

if [[ $# != 2 ]]; then
	echo "USAGE: expand_networks.sh drmn_results_dir celltype_order.txt"
	exit 2
fi

resdir=$1
order=$2

if [[ ! -d $resdir ]]; then
	echo "Cannot find DRMN result dir $resdir"
	exit 2
fi
if [[ ! -e $order ]]; then
	echo "Cannot find celltype order file $order"
	exit 2
fi

# loop over cell types
while read cell
do
	innet=${resdir}/${cell}_net.txt
	inmod=${resdir}/${cell}_clusterassign.txt 
	outnet=${resdir}/${cell}_expanded_net.txt

	hasfiles=0
	for f in $innet $inmod;
	do
		if [[ -e $f ]]; then
			hasfiles=$((hasfiles+1))
		#else
		#	echo "No file $f"
		fi
	done
	# no files for this celltype 
	if [[ $hasfiles < 2 ]]; then
		continue 
	fi

	# to join correctly, we need to have the module assignments
	# sorted by module (col 2) (first sort command)
	# and the network file sorted by module integer ID (originally was column 1
	# after the join, the format is: moduleID target regulator weight.
	# The network format is regulator, target, weight.
	# All weights for the same module are the same.
	join -t$'\t' -1 2 -2 1 <(sort -k2,2n $inmod ) <(cat $innet | sed 's/Module//' | sort -k1,1n) | awk '{print $3"\t"$2"\t"$4}' > $outnet
	echo $outnet
done < $order

