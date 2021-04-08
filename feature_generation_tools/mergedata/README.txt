# mergeData: will vertically join data matrices on the same row IDs.

### Compilation: g++ Framework.C -g -o mergeData

### Usage:
```
mergeData infileset outputexpr
```

### Example(s)

Can be issued from above directory as in Example.sh

```
mergedata/mergeData example_files/mouse_chr1_example_list.txt example_files/mouse_chr1_example_merged.txt
mergedata/mergeData example_files/mouse_chr1_example_promoter_list.txt example_files/mouse_chr1_example_promoter_merged.txt
```

Where the input and output is correspondingly:

```
(from ../example_files/mouse_chr1_example_list.txt)
example_files/mouse_ESC_chr1_example_Q_motif_aggregated_values.txt
example_files/mouse_MEF_48hr_chr1_example_Q_motif_aggregated_values.txt
example_files/mouse_MEF_chr1_example_Q_motif_aggregated_values.txt
example_files/mouse_pips1_chr1_example_Q_motif_aggregated_values.txt
example_files/mouse_pips2_chr1_example_Q_motif_aggregated_values.txt
```

and

```
(from ../example_files/mouse_chr1_example_merged.txt)
AAANWWTGC_ENSMUST00000006570	<nodata>	3.08405	2.92899	<nodata>	<nodata>
AAANWWTGC_ENSMUST00000010049	2.63752	40.9739	42.633	26.5423	54.8784
AAANWWTGC_ENSMUST00000012331	<nodata>	3.96521	2.92899	<nodata>	<nodata>
AAANWWTGC_ENSMUST00000015987	2.24189	3.96521	1.62721	<nodata>	<nodata>
GTGGGTGK_ENSMUST00000111815	6.97294	3.96521	5.85797	1.10593	5.87983
GTGGGTGK_ENSMUST00000111836	15.5778	4.46086	10.9837	47.5549	33.319
GTGGGTGK_ENSMUST00000111887	16.1713	11.8956	4.7596	29.86	49.9785
GTGGGTGK_ENSMUST00000112232	16.468	2.47826	4.7596	82.9445	23.5193
GTGGGTGK_ENSMUST00000112538	39.1671	101.113	84.5745	143.771	97.9971
```

To select elements with uniform data representation you can apply "grep -F -v nodata" to the merged data file, as exemplified in Example.sh.

### Summary

infileset: a text file with one filename per row
outputexpr : a new filename for the merged result

Example input file 1:
	Gene	Sample1	Sample2
	gene1	0.5	0.7
	gene2	0.0	0.2

Example input file 2:
	Gene	Sample3
	gene1	0.0
	gene2	1.1

Example input file list:
	file1.txt 
	file2.txt

Example output file:
	Gene	Sample1	Sample2	Sample3
	gene1	0.5	0.7	0.0
	gene2	0.0	0.2	1.1
	
### Final notes

NB1: In order to get headers in the output file, you must have "Gene" as the entry in the very first row/column for every file.
NB2: If a gene/region label is unobserved in one file, it will receive the value "<nodata>".

