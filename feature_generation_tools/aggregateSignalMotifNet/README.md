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

### Example
To run the example, call ./run_example.sh at your bash prompt.






