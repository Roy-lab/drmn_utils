mergeData will vertically join data matrices on the same row IDs.

Compilation:g++ Framework.C -g -o mergeData

Usage:
>> ./mergeData
mergeData infileset outputexpr

infileset: a text file with one filename per row to define the set of input files to merge columnwise.
outputexpr : a new filename for the merged data file result.

Example input file 1:
	Gene	Sample1	Sample2
	gene1	0.5	0.7
	gene2	0.0	0.2

Example input file 2:
	Gene	Sample3
	gene1	0.0
	gene2	1.1


Example output file:
	Gene	Sample1	Sample2	Sample3
	gene1	0.5	0.7	0.0
	gene2	0.0	0.2	1.1

NB1: In order to get headers in the output file, you must have "Gene" as the entry in the very first row/column for every file.
NB2: If a gene is unobserved in one file, it will receive the value "<nodata>".

