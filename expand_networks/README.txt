expand_networks.sh

A simple tool to combine the module assignments with the regulatory programs to generate a regulator-gene network.
This network is weighted by regression weight. The weight for one regulator to all genes in a module is the same.

USAGE: ./expand_networks.sh drmn_results_dir celltype_order.txt

results_dir contains the DRMN output.
celltype_order.txt lists the cell types.

Example:
./expand_networks.sh example celltype_order.txt
