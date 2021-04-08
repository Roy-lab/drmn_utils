# aggregateSignalMotifNet
This tool is used to aggregate counts for edges in a motif->gene network.

*Given:*
* Regions: motif instances assigned to TSSes for genes, eg by using a window around the TSS. We can have many individual motifs and many TSSes. Each TSS maps to one gene, although each gene may have many TSSs.
* Signal file: counts file output by, eg, bedtools genomecov
*Do:* 
* Score motif-gene edges by aggregated signal under motif instances for the gene.
* Optional: print count profiles around motif instances.


### How scoring works
The score for a motif->TSS edge is the sum of counts under the motif instances for that TSS. 
The score for a motif->gene edge is the maximum motif->TSS value for all TSSes for that gene.

### Usage
```
aggregateSignalToSum tss chromosomesizes signal outfilePrefix [-p/--printProfile <windowSize>, -n/--normalize <scalar>]
If -p/--printProfile option, will print per-base-pair signal for motif instances. Use with caution as it will create huge files. 
  This requires a parameter for defining a uniformly sized window around each motif center: center +/- window.
	You need to set the window size at least as wide as half of your largest motif.
If -n/--normalize <scalar>, will normalize to a sequencing depth of <scalar>*1 million reads.
```
### Input formatting

The <tss> argument will be a set of motifs mapped to genes in the following format:
```
(from ../example_files/mouse_chr1_example.txt)
chr1	CONVERT	TSS_ENSMUST00000160944_TSS1	3042328	3042337	0.000006	+	.	TGANNYRGCA;ENSMUST00000160944
chr1	CONVERT	TSS_ENSMUST00000160944_TSS1	3042487	3042500	0.000007	+	.	LM186;ENSMUST00000160944
chr1	CONVERT	TSS_ENSMUST00000160944_TSS1	3042531	3042544	0.000000	+	.	LM32;ENSMUST00000160944
```
The chromosomesizes input is a standard genome .fai file:
```
(from ../example_files/mm9.fa.fai)
chr1	197195432	6	50	51
chr2	181748087	201139353	50	51
chr3	159599783	386522408	50	51
```	
The signal file is a .bedGraph format file of read coverage (as from bedtools genomecov) in the following format:
```
(from ../example_files/mouse_ESC_chr1_example.counts)
chr1	3000463	3000514	1
chr1	3000621	3000672	1
chr1	3000938	3000989	1)
```	
Finally the prefix argument is a prefix intended to organize the output file names. 	
	
### Example

To run an example try:

```
	./aggregateSignalMotifNet ../example_files/mouse_chr1_example.txt ../example_files/mm9.fa.fai ../example_files/mouse_ESC_chr1_example.counts ../:wexample_files/mouse_ESC_chr1_example_Q_motif -n 1
```
### Outputs:

*_aggregate.txt*
*aggregated_values.txt*

### Preparing input file of regions mapped to genes

A python utility script matchMotifToGenePerTF2.py is included for organizing the first input to the aggregageSignalMotifNet progra). The usage is as follows:
```
python matchMotifToGenePerTF2.py <motif regions> <tss list> <upstream window> <downstream window> <output>
```
Here the motif instances are anticipated in a format as follows:
chr1   3004298 3004313 +   8.830266

The list of gene TSS sites are indicated to be listed as follows:
ENSMUST00000166088 chr1   75032585    75032585

The produced list of mappings can then be used as input for signal aggregation in the format for aggregateSignalMotifNet above. The specific usage for this working example is:

```
python matchMotifToGenePerTF2.py ../example_files/mouse_chr1_motif_example.txt ../example_fil    es/Mus_musculus.NCBIM37.67.TSS.txt 2500 2500 ../example_files/mouse_chr1_example.txt
```

The mapping is done according to a defined upstream and downstream window relative to a gene TSS site. This window can be asymetric. A motif site within the defined window will be mapped to that respective gene. Note genes and motifs are in general multi-way mapped because a given motif site may fall within the specified TSS window of multiple genes,.




