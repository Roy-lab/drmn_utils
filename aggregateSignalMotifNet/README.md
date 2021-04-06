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
aggregateSignalMotifNet tss chromosomesizes signal outfilePrefix [-p/--printProfile <windowSize>, -n/--normalize <scalar>]
If -p/--printProfile option, will print per-base-pair signal for motif instances. Use with caution as it will create huge files. 
  This requires a parameter for defining a uniformly sized window around each motif center: center +/- window.
	You need to set the window size at least as wide as half of your largest motif.
If -n/--normalize <scalar>, will normalize to a sequencing depth of <scalar>*1 million reads.
```

### Example
To run the example, call ./run_example.sh at your bash prompt.

### Preparing input file of regions mapped to genes

A python utility script is included here, matchMotifToGenePerTF2.py, for organizing an inputinformation of motif regions mapped to genes (the first input to this program). The usage is as follows:
```
python matchMotifToGenePerTF2.py <motif regions> <tss list> <upstream window> <downstream window> <output>
```

Here the motif instances are anticipated in a format as follows:
chr1   3004298 3004313 +   8.830266

The list of gene TSS sites are indicated to be listed as follows:
ENSMUST00000166088 chr10   75032585    75032585

The produced list of mappings can then be used as input for signal aggregation. 

The mapping is done according to a defined upstream and downstream window relative to a gene TSS site, which can be asymetric, and a motif site within that region will be mapped to that gene if this criterion is met. Genes and motifs are multi-way mapped because a given motif site may fall within this specified window for multiple genes, and multiple motif sites can in general be mapped to a single gene.





