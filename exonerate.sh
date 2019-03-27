#!/bin/bash
#run gth without 
genome=$(realpath $1)
protein=$(realpath $2)
cpu=$3
guidefile=$(realpath $4) #by DIAMOND

#exonerate gt tabix samtools
EXONERATEHOME=/dirpath/to/exonerate
export PATH=/path/to/gt:$PATH

pwd=$(pwd)

[[ -f ${genome}.fai ]] || samtools faidx $genome
replace_str=".addstop" ;protein=${protein%%replace_str}
[[ -f ${protein}.addstop ]] || gt seqtransform -addstopaminos $protein > ${protein}.addstop
protein=${protein}.addstop
[[ -f ${protein}.fai ]] || samtools faidx $protein

mkdir -p $pwd/errs
mkdir -p $pwd/aln
> command.list

for x in $(tabix -l $guidefile)
do  
  #all the operations take place in the folder created in /tmp  directory  to avoid lock problems in nfs storage
  tmpdir="mktemp -dp /tmp/";
  cmd="tmpdir=\$($tmpdir);cd \$tmpdir;trap 'rm -rf "\$tmpdir"' EXIT; cd \$tmpdir";
  cmd="$cmd; samtools faidx $protein '$x' > pro.fa ; samtools faidx $genome \$(tabix $guidefile '$x' |cut -f5|awk '{s=s\" \"\$1} END{print s}') > genome.fa ";
  cmd="$cmd && $EXONERATEHOME/bin/exonerate --model protein2genome --fsmmemory 1024 --showalignment no --percent 80  --showsugar no --showvulgar no --showtargetgff yes -q pro.fa -t  genome.fa > '$pwd/aln/exon.result.${x}' 2> '$pwd/errs/exon.stderr.${x}' "  
  cmd="$cmd && rm -rf \$tmpdir"
  cmd="test -f '$pwd/aln/gth.result.${x}' ||($cmd)"
  echo $cmd >> command.list 
done

cmd="parallel -j $cpu < command.list && cat $pwd/aln/gth.result\* > $pwd/exon.aln.concat && cat $pwd/errs/gth.stderr*> $pwd/exon.stderr.concat"
eval $cmd
