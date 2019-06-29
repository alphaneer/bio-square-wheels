#!/bin/bash

readonly help_msg="\
Usage: $(basename $0) [OPTION]
`basename $0` is used to run blasr efficiently. It first run minimap2 and then blasr
eg `basename $0` -t 3 -q qry.fa  -r ref.fa -R minimap.paf -x PAF/MAF -o out;  
the outdir are as below when finished
outdir
|-- minimap.paf (must exist to skip the time-consuming generating process,let it an empty file when minimap.paf.gz provided)
|-- minimap.paf.gz
 -- minimap.paf.gz.tbi
"
readonly last_header="\
# LAST version 876
#
# a=7 b=1 A=7 B=1 e=33 d=24 x=32 y=9 z=32 D=1e+06 E=161.501
# R=01 u=2 s=2 S=0 M=0 T=0 m=10 l=1 n=10 k=1 w=1000 t=0.910239 j=3 Q=0
# ../ref/human_cR01db
# Reference sequences=594 normal letters=3095950843
# lambda=1.09602 K=0.335388
#
#    A  C  G  T
# A  1 -1 -1 -1
# C -1  1 -1 -1
# G -1 -1  1 -1
# T -1 -1 -1  1
#
# Coordinates are 0-based.  For - strand matches, coordinates
# in the reverse complement of the 2nd sequence are used.
#
# name start alnSize strand seqSize alignment
#
# m=0.01 s=33
#
"
readonly get_region=$(cat <<'COMMAND'
awk '{\
    if(! ($6 in min_map)) {min_map[$6]=$8;max_map[$6]=$9;next;} if($9>max_map[$6]) max_map[$6]=$9;if($8<min_map[$6]) min_map[$6]=$8;
  } END{
    s="";for(i in min_map) s=s" "i":"(min_map[i]>5000?min_map[i]-5000:0)"-"max_map[i]+5000;
    print s;
  }'
COMMAND
)

parse_args () {
    [[ "$#" == 0 ]] && show_msg "missing file argument" && exit -1;
    local OPTIND=1
    while getopts "t:q:r:R:o:x:" opt "$@"; do
        case "$opt" in
            q) qry=$(realpath -s "$OPTARG")
               ;;
            r) ref=$(realpath -s "$OPTARG")
               ;;
            t) ncpu="$OPTARG"
               ;;
            R) minimap_out=$(realpath -s  "$OPTARG")
               ;;
            o) outdir=$(realpath "$OPTARG")
               ;;
            x) type="$OPTARG"
               ;;
            h) show_msg
               exit 0
               ;;
            *) show_msg "Try -h for more information."
               exit 1
               ;;
        esac
    done
    shift $((OPTIND-1))
}

show_msg () {
    [[ "$1" ]] && echo "`basename $0`: $1" >&2
    echo "$help_msg" >&2
}

parse_args "$@"

pwd=$(pwd)
[[ $outdir == "" ]] && outdir="." 
mkdir -p $outdir

[[ ! -e ${ref}.fai ]] && samtools faidx $ref
[[ ! -e ${qry}.fai ]] && samtools faidx $qry

if [[ $minimap_out == "" ]]
then
   echo $last_header > ${outdir}/minimap.raw.maf
   minimap2 -t $ncpu -k 15 -w 5 -c --cs=long -r2k $ref $qry 2>/dev/null|k8 $SWPATH/minimap2/misc/paftools.js view -f maf - |tail -n +2 >> ${outdir}/minimap.raw.maf
   last-split ${outdir}/minimap.raw.maf > ${outdir}/minimap.maf 
   samtools view -@ $ncpu -ht ${ref}.fai  <(maf-convert sam ${outdir}/minimap.maf)|k8 $SWPATH/minimap2/misc/paftools.js sam2paf - > ${outdir}/minimap.paf 
   minimap_out=${outdir}/minimap.paf
   type="PAF"
elif [[ ! -e ${outdir}/minimap.paf  ]];then
   if [[ $type == "PAF" ]];then
      ln -sf $minimap_out ${outdir}/minimap.paf
      minimap_out=${outdir}/minimap.paf
   elif [[ $type == "MAF" ]];then
      samtools view -@ $ncpu -ht ${ref}.fai  <(maf-convert sam $minimap_out)|k8 $SWPATH/minimap2/misc/paftools.js sam2paf - > ${outdir}/minimap.paf
      minimap_out=${outdir}/minimap.paf
   fi
else  
   minimap_out=${outdir}/minimap.paf
fi

echo $qry >&2
echo $ref >&2
echo $minimap_out >&2

[[ ! -e ${outdir}/minimap.paf.gz ]] && cat ${outdir}/minimap.paf|sort -k1,1 -k3,3n -k4,4n |bgzip -@ $ncpu > ${outdir}/minimap.paf.gz
[[ ! -e ${outdir}/minimap.paf.gz.tbi ]] && tabix -f -s 1 -b 3 -e 4 ${outdir}/minimap.paf.gz


for q in $(tabix -l ${outdir}/minimap.paf.gz)
do
  
 cmd="region=\$(tabix ${outdir}/minimap.paf.gz $q |$get_region)" 
 cmd="cd $outdir; $cmd && samtools faidx $qry $q > ${q}.qry.fa && samtools faidx $ref \$region > ${q}.ref.fa \
&& blasr ${q}.qry.fa ${q}.ref.fa  --nproc 1 --noSplitSubreads --bestn 1 -m 1 --out ${q}.out 1>/dev/null 2>/dev/null && rm -rf ${q}.*.fa"
  echo $cmd 
done
