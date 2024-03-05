#!/bin/bash

# Run cellranger version 7.0.1

for sample in D19-4295 D19-4296 D19-4297 D19-4298 D19-4299 D19-4300 D19-4301 D19-4302 D19-5859 D19-5860 D19-5861 D19-5862 D19-5863 D19-5864 D19-5865 D19-5866 rep-D19-4295 rep-D19-4296 rep-D19-5859 rep-D19-5860 ; do
	cellranger count \
		--id="$sample" \
		--sample="$sample" \
		--fastqs=/home/leandro/Documents/IC/BipolarDisorder/snRNAseq/PsychENCODE/ \
		--transcriptome=/home/leandro/Documents/IC/BipolarDisorder/snRNAseq/PsychENCODE/human_ref/refdata-gex-GRCh38-2020-A
		
done
