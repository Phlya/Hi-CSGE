# Hi-CSGE
A pipeline fastq->heatmaps for an SGE cluster based on `hiclib`. In the end you can get 'cooler' files and files ready to be converted to .hic format for visualisation in Juicebox.

The idea is to split input fastq files into chunks and process them in parallel, then once we have fragment files, merge them for each sample and apply filters all together, and then create final heatmaps from the combined file.

Tested on an SGE cluster (eddie3, University of Edinburgh).
