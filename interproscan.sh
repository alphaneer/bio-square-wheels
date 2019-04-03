#!/bin/bash
#run interproscan by spliting the genome into contig blocks and run it one after another
protein=$(realpath -s $1)
ncontigs=$2 # the number of contigs running at one time

[[ -f ${protein}.fai ]] || samtools faidx $protein

#samtools
INTERPROSCANHOME=$SWPATH/interproscan-5.30-69

pwd=$(pwd)

[[ -f ${protein}.fai ]] || samtools faidx $protein

mkdir -p $pwd/errs
mkdir -p $pwd/res

cut -f1 ${protein}.fai |xargs -n $ncontigs|while IFS= read -r x
do
  tmpdir="mktemp -dp $pwd";
  cmd="tmpdir=\$($tmpdir);trap 'rm -rf "\$tmpdir"' EXIT; cd \$tmpdir";
  cmd="$cmd; samtools faidx $protein $x > pro.fa"
  #make the output file name short to avoid system error info 'Too long parameter'
  result=$(awk -v chrs="$x" 'BEGIN{n=split(chrs,tA," ");print tA[1]".."tA[n]}')
  cmd="$cmd; $INTERPROSCANHOME/interproscan.sh -goterms -pa -i pro.fa -o $pwd/res/res.${result} -f tsv 2> $pwd/errs/stderr.${result}"
  cmd="$cmd && rm -rf \$tmpdir;"
  cmd="test -f '$pwd/res/res.${result}' ||($cmd)"
  echo $cmd
done
