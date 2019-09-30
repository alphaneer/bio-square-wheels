#!/bin/bash
REF=$(realpath -s $1)
READS=$(realpath -s $2) #fq.fofn:pe1 pe2
CPU=$3
ROUND=$4
BUSCO_CLS=$5

[ -z $BUSCO_CLS ]  && BUSCO_CLS=embryophyta_odb9

for i in $(seq 1 $ROUND)
do
  OUTDIR="$(readlink -f ngs_round${i})"
  mkdir -p $OUTDIR
  ln -sf $REF* $OUTDIR/
  REF=${OUTDIR}/$(basename $REF)  
  
  [ ! -f ${REF}.fai ] && samtools faidx ${REF}

  #bwa  
  [ ! -f ${REF}.bwt ] && bwa index ${REF}

  bamlist=""
  while read pe
  do
    s=$(basename $pe)
    cmd="bwa mem -t $CPU $REF $pe | samtools sort -@ $[$CPU/2+1] -m 2g -o ${OUTDIR}/${s}.bam && samtools index -@ $CPU ${OUTDIR}/${s}.bam"
    echo $cmd
    bamlist="$bamlist ${OUTDIR}/${s}.bam"
  done <<<"$(cat $READS)"  > ${OUTDIR}/bwa.sh

  tmp_var_outdir=${OUTDIR}/tmp_var_outdir
  mkdir -p $tmp_var_outdir
  
  #samtools view -H $bam | grep "\@SQ" | sed 's/^.*SN://g' | cut -f 1
  
  parallel --dry-run --will-cite --colsep '\t' bcftools mpileup -Ou -f $REF $bamlist -r {1} '|' bcftools call -mv -o $tmp_var_outdir/{1}.vcf :::: ${REF}.fai |tac > ${OUTDIR}/snp.sh
  
  cnsfa="${OUTDIR}/$(basename $REF).cns.fa"  
  cmd="cd $OUTDIR; \
~/bin/python3 $SWPATH/busco/scripts/run_BUSCO.py -o $(basename $cnsfa) -i $cnsfa -l \
$SWPATH/busco/lineage/${BUSCO_CLS}/ -m genome -c $CPU"
  echo $cmd > ${OUTDIR}/busco.sh
  REF=$cnsfa
done
