#!/bin/bash

readonly help_msg="\
Usage: $(basename $0) [OPTION]
`basename $0` is used to run allhic
eg `basename $0` -r ref.fa -b bam -t 20 -k 12 -m 25 -e AAGCTT  
-b FILE 	bam file can be generated with command:
          	   bwa mem -SP5M -t 20 ref.fa r1.fq.gz r2.fq.gz |samtools view -Shb - -o out.bam
                unreliable alignments in bam are filtered out, skip the bam filter by command:
		   ln -sf out.bam out.REduced.bam
-r FILE		reference fasta file
-e STRING	enzyme_sites HindIII: AAGCTT; MboI: GATC
-k NUMBER       number of groups (user defined K value)
-m INT		minimum number of restriction sites (default, 25)
-t INT		max threads
"
parse_args () {
    [[ "$#" == 0 ]] && show_msg "missing file argument" && exit -1;
    local OPTIND=1
    while getopts "r:b:t:k:m:e:h" opt "$@"; do
        case "$opt" in
            r) ref=$(realpath -s "$OPTARG")
               ;;
            b) bam=$(readlink -f "$OPTARG")
               ;;
	    t) ncpu="$OPTARG"
               ;;
	    k) chr_cnt="$OPTARG"
               ;;
	    m) min_restriction_sites="$OPTARG"
	       ;;
            e) enz="$OPTARG"
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

ncpu=${ncpu:-20}

min_restriction_sites=${min_restriction_sites:-25}

enz=${enz:-AAGCTT}

#env
export PATH=$SWPATH/ALLHiC/scripts/:$SWPATH/ALLHiC/bin/:$PATH

s=$(basename $bam .bam)
ln -s $ref .
ref=$(basename $ref)

BED_RE_file=${ref}.near_${enz}.500.bed
cmd="`which make_bed_around_RE_site.pl` $ref $enz 500"
cmd="[ ! -s $BED_RE_file ] && $cmd "
eval $cmd


cmd="samtools view $bam -L $BED_RE_file -q 30 -h -F 0x0c -@ $ncpu |awk '{if( (match($0,/NM:i:([0-9]+)/,tA) && tA[1]>5) || (match($0,/XA:Z:.*/))   ) next;print;   }'| samtools view -Sbh -@ $ncpu > ${s}.REduced.bam"
cmd="[ ! -s  ${s}.REduced.bam ] && $cmd "
eval $cmd

cmd="samtools flagstat -@ $ncpu ${s}.REduced.bam > ${s}.REduced.bam.flagstat & "
cmd="[ ! -s ${s}.REduced.bam.flagstat ] && $cmd"
eval $cmd

cleanbam=${s}.REduced.bam

if ! (ls ${s}.REduced.counts_${enz}*g*.txt && ls ${s}.REduced.clm)
then
  `which ALLHiC_partition` -b $cleanbam -r $ref -e $enz -k $chr_cnt -m $min_restriction_sites
  `which allhic` extract $cleanbam $ref --RE $enz
fi

parallel_cnt=$ncpu
[ $parallel_cnt -gt $ch_cnt ] && parallel_cnt=$chr_cnt
for x in $(ls ${s}.REduced.counts_${enz}*g*.txt)
do
   cmd="`which allhic` optimize $x ${s}.REduced.clm"
   echo $cmd
done | parallel -j $parallel_cnt

`which ALLHiC_build` $ref

seq_len.pl groups.asm.fasta |grep 'SP5M.rmdup.REduced.counts_AAGCTT' > chrn.list
`which ALLHiC_plot` $cleanbam groups.agp chrn.list 500k pdf
