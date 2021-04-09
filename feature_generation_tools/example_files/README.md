# Summary

This directory contains all the input, intermmediate and output files for the Q-motif and promoter accessibility feature generation pipeline. The formatting of these files has been generally described in the readme files of the main directory and those of the aggregateSignalMotifNet and aggregateSignalRegion tools. The included chart below summarizes the key contents in this directory, and how they were generated in the implementation presented in ../Example.sh.

![Summary flow chart for example feature processing analysis](https://github.com/Roy-lab/drmn_utils/blob/master/feature_generation_tools/example_files/FileFlowChartSummary.pdf?raw=true)

# Decompressing files

Several files have been uploaded in a .tar.xz format to make this directory more tractable to utilize/share. These are in particular:

```
mouse_chr1_example_merged.tar.xz
mouse_chr1_example.tar.xz
mouse_ESC_chr1_example_Q_motif_aggregate.tar.xz
mouse_ESC_chr1_example.tar.xz
mouse_MEF_48hr_chr1_example.tar.xz
mouse_MEF_chr1_example.tar.xz
mouse_pips1_chr1_example.tar.xz
mouse_pips2_chr1_example.tar.xz
```
Each of these were generated in the following way:

```
tar cf - filename | xz -4e > filename.tar.xz
```

and they can be readilly compressed with:

```
tar -xf archive.tar.xz
```
This has been tested with tar (GNU tar) 1.26. 
