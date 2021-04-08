# Summary of Q-motif and promoter accessibility feature generation. 

  This directory contains code and scripts for generating and normalizing Q-motif and promoter accessibility feature data for input to DRMN. The steps shown here (outlined in Example.sh) will walk through the several steps necessary to aggregate and normalize ATACseq accessibility signals for use in DRMN. The files in example_files are meant to represent the full feature preparation pipeline for the chr1 subset of the Chronis et al. data set. The goal is for the interseted user to be able to run this example feature-generation analysis as implemented with Example.sh.

This directory contains code for three C++ utilities (aggregateSignalMotifNet, aggregateSignalRegion, and mergeData). To run the main example analysis in Example.sh you will need to first compile each of these. The first two (aggregateSignalMotifNet, and aggregateSignalRegion) include make files, and mergeData can be compiled with a simple command (g++ Framework.C -g -o mergeData, also noted in the corresponding readme). Running the complete analysis presented in Example.sh analysis requires a working matlab and python installation. With these requirements met the Example.sh analysis on the included data completes in the time scale of a few minutes, 

### Prepare coverage data from an aligned (ATACseq) library (.bam file). 

  The bedtools genomecov tool is an expedient choice to generate coverage data from an aligned data set, which in this work are generally ATACseq data tracing genome acessibility or chromatin mark ChIPseq data. Here we use the "-bg" option to generate genome-wide coverage data in the bedGraph format. When assessing paried end data the "-pc" option has been applied to prepare coverage based on full fragment lengths of each pair. The format for such data will be as follows:
  
```
chr1	3000857	3000906	1
chr1	3001108	3001159	1
chr1	3001583	3001634	1
```

Such data has been compiled for the chr1 coverage data from the Chronis et al. data set in the example_files directory:

```
example_files/mouse_ESC_chr1_example.counts
example_files/mouse_MEF_48hr_chr1_example.counts
example_files/mouse_MEF_chr1_example.counts
example_files/mouse_pips1_chr1_example.counts
example_files/mouse_pips2_chr1_example.counts
```
What follows is an example tutorial on aggregating and organizing such data into Q-motif and promoter signal feature data inputs for DRMN. 

### Aggregate (score) motif instances for Q-motif features

Given a set of motif instances (such as those listed in example_files/mouse_chr1_motif_example.txt), we first map those instances to gene TSS sites and aggregate the data for those regions. This is done in two steps first, using the aggregateSignalMotifNet/matchMotifToGenePerTF2.py script to apply that mapping, and the aggregateSignalMotifNet program to generate aggregated data for these motif instances for each data set. 

The usage and formatting infromation for these two tools is presented in the master Example.sh script, and in particular in the README for the aggregateSignalMotifNet tool. Note when using the -n1 argument for aggregateSignalMotifNet the output data values will be divided by the mean per-bp, genome-wide coverage of signal, but will not be further normalized.

From this step a set of data files with unique identifiers representing motif and gene pairings is obtained:
```
example_files/mouse_ESC_chr1_example_Q_motif_aggregated_values.txt <==
AAANWWTGC_ENSMUST00000000834	2.50564
AAANWWTGC_ENSMUST00000001284	2.76939
AAANWWTGC_ENSMUST00000003219	1.18688
AAANWWTGC_ENSMUST00000006467	3.16502
AAANWWTGC_ENSMUST00000010049	2.63752
```

Where each row has a signifier for <motif>_<gene> pairings, and an aggregated signal value. In this working example such inputs can be found in the following files:

```
example_files/mouse_ESC_chr1_example_Q_motif_aggregated_values.txt
example_files/mouse_MEF_48hr_chr1_example_Q_motif_aggregated_values.txt
example_files/mouse_MEF_chr1_example_Q_motif_aggregated_values.txt
example_files/mouse_pips1_chr1_example_Q_motif_aggregated_values.txt
example_files/mouse_pips2_chr1_example_Q_motif_aggregated_values.txt
```

### Merge and normalize Q-motif feature data

The aggregated signal values (already likely normalized by the mean per-bp coverage in each respective data set) are then merged, log transformed and quantile normalized. Such steps can indeed be applied flexibly depending on user preferences but here we've included the mergeData tool and the quantile_normalize.m script for one working option, with usage shown in Example.sh. 

Here the relevant files of interest (indicated in Example.sh) are:
```
example_files/mouse_chr1_example_merged.txt - motif singal data merged across conditions
example_files/mouse_chr1_example_merged_uniform.txt - merged rows with no missing data
example_files/mouse_chr1_example_merged_normalized.txt - example data matrix with final log-transformed and quantile normalized values
```
### Aggregate (score) and normalize promoter signal data

Similarly, as the set of Q-motif data aggregate and normalized based on an initial set of regions, we'll first define a set of promoter regions (here in example_files/promoters.bed) using the aggregateSignalRegion tool, and likewise merge and normalize these data separately, but in analogy to what has been shown for Q-motifs. The appropriate usage and formats are again exemplified in Example.sh (lines 29-45), the data in example_files, and the readme for the aggregateSignalRegion tool.

```
example_files/mouse_chr1_example_promoter_merged.txt - promoter signal data merged across conditions
example_files/mouse_chr1_example_promoter_merged_uniform.txt - merged rows with no missing data
example_files/mouse_chr1_example_promoter_merged_normalized.txt - example data matrix with final log-transformed and quantile normalized values
```

### Prepare final feature input data files for DRMN 

The input feature files for DRMN include both the Q-motif and promoter signal data as well as the number of regulatory features and genes included. The example feature files generated with Example.sh are:
```
example_files/mouse_ESC_chr1_example_features.txt
example_files/mouse_MEF_48hr_chr1_example_features.txt
example_files/mouse_MEF_chr1_example_features.txt
example_files/mouse_pips1_chr1_example_features.txt
example_files/mouse_pips2_chr1_example_features.txt
```

There is one file for each condition with the following format:

```
795	5158
AAANWWTGC	ENSMUST00000010049_ESC	1.683187
AAANWWTGC	ENSMUST00000027139_ESC	2.504774
AAANWWTGC	ENSMUST00000027263_ESC	2.405853
...
ATAC	ENSMUST00000177472_ESC	2.805337
ATAC	ENSMUST00000177501_ESC	1.175440
ATAC	ENSMUST00000177532_ESC	1.175743
```
The top two numbers are the number of regulatory features (the number of motifs represetned plus one for the ATAC promoter signal features) and the feature names designated here are the motif names and "ATAC" for the promoter signals. Here the gene names are also given a unique identifier for the condition the respective data points represent, here the "ESC" condtion for this working example data set. 

Such feature files are then utilized as direct input to DRMN along corresponding expression and cluster assignment data. 



