#!/bin/bash
genome=$(realpath $1)
protein=$(realpath $2)
cpu=$3
guidefile=$(realpath $4) #by DIAMOND
GTHHOME=/path/to/GenomeThreader/
export PATH=/path/to/gt/bin:$PATH
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
  cmd="tmpdir=\$($tmpdir);trap 'rm -rf "\$tmpdir"' EXIT; cd \$tmpdir";
  cmd="$cmd; samtools faidx $protein '$x' > pro.fa ; samtools faidx $genome \$(tabix $guidefile '$x' |cut -f5|awk '{s=s\" \"\$1} END{print s}') > genome.fa ";
    
  cmd="$cmd && $GTHHOME/bin/gth -genomic genome.fa -protein pro.fa -gff3out -skipalignmentout -gcmincoverage 80 -prseedlength 20 -prminmatchlen 20 -prhdist 2 -o '$pwd/aln/gth.result.${x}' 2> '$pwd/errs/gth.stderr.${x}'"
  cmd="$cmd && rm -rf \$tmpdir"
  cmd="test -f '$pwd/aln/gth.result.${x}' ||($cmd)"
  echo $cmd >> command.list 
done

cmd="parallel -j $cpu < command.list && cat $pwd/aln/gth.result\* > $pwd/gth.aln.concat && cat $pwd/errs/gth.stderr*> $pwd/gth.stderr.concat"
eval $cmd
