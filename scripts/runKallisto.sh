#!/bin/bash

cd /Volumes/target_nbl_ngs/PPTC-PDX-genomics/bcm-processed/rnaseq_fastq/fastq

for sample in ; do kallisto quant \
-i /Volumes/target_nbl_ngs/Kendsersky/b7h3/rnaseq_data/reference/Homo_sapiens.GRCh38.cdna.all.index \
-o /mnt/isilon/maris_lab/target_nbl_ngs/Kendsersky/b7h3/rnaseq_data/aln/${sample}.aligned \
-t 12 \
"${sample}_R1.fastq.gz" \
"${sample}_R2.fastq.gz" \
&> /Volumes/target_nbl_ngs/Kendsersky/b7h3/rnaseq_data/aln/logs/${sample%%.*}.aln.log \
; done