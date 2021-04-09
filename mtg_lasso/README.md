# Preparing inputs for MTG-LASSO

The python script createRegMat.py is used to create the inputs:

```
Usage: python createRegMat.py --config config.txt --order order.txt --clusters all_genesets.txt --outprefix data/out
```

where config.txt is config file input of DRMN, order.txt is order of cell lines/time points input of DRMN, and data/out is the prefix of the output files.
all\_genesets.txt is the list of genes in transitioning gene sets.
The format is, first column is gene name, second column is the gene set number
```
gene1	1
gene2	1
gene3	2
gene4	2
``` 

The matlab script runMTGLASSO.m is the wrapper script that will uses SLEP package to perform MTG-LASSO for one transitioning gene set.
SLEP is available at http://yelabs.net/software/SLEP/

The same generated input files can be used to apply linear regression (predicting the expression that is the first column of the data matrix as a function of regulatory features that are the remaining columns of the data matrix).
