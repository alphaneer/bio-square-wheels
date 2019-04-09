#!/bin/bash
samplename=$1
raw_pe1=$2
raw_pe2=$3
ref=$4 #bwa indexed
cpu=$5

~/bin/awk -v pe1=$raw_pe1 -v pe2=$raw_pe2 -v samplename=$samplename -v ref=$ref -v cpu=$cpu '\
BEGIN{
  outdir=ENVIRON["PWD"]"/"samplename;
  ID=pe1;gsub(".*/","",ID);
  cleanreads=sprintf("fastp -w %s -i %s -I %s --stdout -l 50 -j %s/%s.fasp.json -h %s/%s.fasp.html",cpu/2,pe1,pe2,outdir,ID,outdir,ID);
  cmd=sprintf("bwa mem  -t %d -R %s -p %s  <(%s) 2>%s/logs/%s.map.bam.log| tee >(samtools  stats - > %s/%s.bam.stats)  >(samtools  flagstat - > %s/%s.bam.flagstat) |samtools view -Sb |bedtools bamtobed -i stdin |sortm --parallel=%d --temporary-directory=./ -k1,1 -k2,2n |pigz > %s/%s.bed.gz ",cpu,"\"@RG\\tID:"ID"\\tLB:L"ID"\\tPL:ILLUMINA\\tSM:"samplename"\"",ref,cleanreads,outdir,ID,outdir,ID,outdir,ID,cpu,outdir,ID);
  print "mkdir -p "outdir"/logs;"cmd;
}'
