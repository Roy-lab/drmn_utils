#!/bin/bash
# Runs example of scoring up motif-gene edges from accessibility
# DC Oct 2017

set -exu

if [[ ! -d example_output ]]; then
	mkdir example_output
fi

# Make a window this wide around the motif instances
# for printing profiles
window=20

# normalize to 1 million reads
normFactor=1000000


# for test, just use two small chromosomes
# based on /mnt/ws/sysbio/roygroup/shared/data_new/human/genomes/hg19_v2.fa.fai
sizes=example_input/hg19_chr20_chr21.fa.fai

signal=example_input/test_counts_motif.txt 

motifnet=example_input/motif_sample.txt
mybase=motif_sample

# do for one regions version
for regions in $motifnet
do
	echo $regions
	outfile=example_output/${mybase} # output filename to make
	cmd="./aggregateSignalMotifNet $regions $sizes $signal ${outfile} -p${window} -n ${normFactor}"
	echo $outfile
	eval "${cmd} > ${outfile}.log"
	#echo "valgrind --track-origins=yes --leak-check=full $cmd"
	#echo "gdb $cmd"

	# plot signal
	plotcmd="plot_signal('${outfile}_signalProfile.txt','$outfile','$mybase','pdf')"
	echo $plotcmd
	matlab -nodisplay -nodesktop < <(echo $plotcmd) >> ${outfile}.log
	echo "Done with $regions"
done
